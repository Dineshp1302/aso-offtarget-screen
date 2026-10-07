#!/usr/bin/env python3
"""
Re-score ASO binding sites against ALT alleles.
 
The premise
-----------
Standard off-target screens align candidate oligos against THE reference
genome. But there is no such thing as the human genome -- there are eight
billion of them. A site that carries 2 mismatches in GRCh38 may carry only 1
in a person homozygous for a common variant at the right position. That ASO
is acceptably specific for most of the population and not for them.
 
This script asks, for every variant overlapping every predicted binding site:
does the ALT allele make the duplex BETTER or WORSE?
 
  variant_gained    : fewer mismatches on ALT than on REF
                      -> an off-target that exists only in carriers
  variant_protected : more mismatches on ALT than on REF
                      -> safer in carriers than the reference implies
  position_shift    : same count, but a mismatch moves into or out of the
                      RNase H central window -> cleavage competence changes
  neutral           : no change in count or central status
 
For sites classed on_target the same arithmetic means something different:
a variant that ADDS mismatches to the intended site reduces efficacy in
carriers. That is reported as ontarget_efficacy_loss.
 
Validation built in
-------------------
Before scoring anything, the script recomputes the REF mismatch count from
sequence and compares it against the NM tag the aligner reported. If these
disagree the coordinates, strand convention, or FASTA build are wrong, and
the script says so loudly rather than emitting plausible nonsense.
 
Scope
-----
SNVs only. Indels shift local coordinates and need haplotype reconstruction;
they are counted and reported, not scored. See README limitations.
 
Usage:
    python variant_rescore.py \
        --hits results/hits_table.tsv \
        --vcf results/offtarget_snvs.vcf.gz \
        --fasta ref/chr18.fa \
        --out results/variant_rescored.tsv \
        --summary results/variant_summary.tsv
"""
 
import argparse
import sys
from collections import Counter
 
import pandas as pd
import pysam
 
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
 
 
def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]
 
 
def compare(expected, observed, is_reverse, length):
    """Mismatch count and ASO-relative positions between two equal-length seqs.
 
    `expected` is the plus-strand sequence the ASO should pair with; for a
    minus-strand hit that is revcomp(ASO). Positions are converted back to
    ASO 5'->3' coordinates, 1-based.
    """
    n = 0
    positions = []
    for i, (e, o) in enumerate(zip(expected, observed)):
        if e != o:
            n += 1
            aso_pos = (length - 1 - i) if is_reverse else i
            positions.append(aso_pos + 1)
    return n, sorted(positions)
 
 
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hits", required=True, help="hits_table.tsv")
    p.add_argument("--vcf", required=True, help="bgzipped, indexed, SNVs only")
    p.add_argument("--fasta", required=True, help="indexed reference FASTA")
    p.add_argument("--central-start", type=int, default=8)
    p.add_argument("--central-end", type=int, default=12)
    p.add_argument("--af-field", default="AF")
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    args = p.parse_args()
 
    central = set(range(args.central_start, args.central_end + 1))
 
    hits = pd.read_csv(args.hits, sep="\t")
    fa = pysam.FastaFile(args.fasta)
    vcf = pysam.VariantFile(args.vcf)
 
    rows = []
    n_validated = n_mismatched_nm = n_indels_skipped = 0
 
    for h in hits.itertuples(index=False):
        aso = h.aso_seq
        if not isinstance(aso, str) or not aso:
            continue
        L = len(aso)
        is_rev = (h.strand == "-")
 
        # Plus-strand sequence the ASO is expected to pair with.
        expected = revcomp(aso) if is_rev else aso
 
        try:
            ref_seq = fa.fetch(h.chrom, int(h.start), int(h.end)).upper()
        except (KeyError, ValueError):
            continue
        if len(ref_seq) != L:
            continue
 
        # ---- validation: does sequence agree with the aligner? ----
        ref_n, ref_pos = compare(expected, ref_seq, is_rev, L)
        if ref_n == h.n_mismatch:
            n_validated += 1
        else:
            n_mismatched_nm += 1
            if n_mismatched_nm <= 3:
                print(f"WARN {h.aso_id} {h.chrom}:{h.start} strand {h.strand}: "
                      f"NM tag says {h.n_mismatch}, sequence says {ref_n}",
                      file=sys.stderr)
 
        # ---- variants overlapping this site ----
        try:
            recs = list(vcf.fetch(h.chrom, int(h.start), int(h.end)))
        except ValueError:
            recs = []
 
        for rec in recs:
            if len(rec.ref) != 1 or any(len(a) != 1 for a in rec.alts or []):
                n_indels_skipped += 1
                continue
 
            offset = rec.pos - 1 - int(h.start)   # VCF is 1-based
            if not (0 <= offset < L):
                continue
 
            # REF allele must agree with the FASTA, else the VCF and the
            # genome build do not match and nothing downstream is valid.
            if ref_seq[offset] != rec.ref.upper():
                print(f"WARN REF mismatch at {rec.chrom}:{rec.pos} "
                      f"VCF={rec.ref} FASTA={ref_seq[offset]}", file=sys.stderr)
                continue
 
            af = None
            if args.af_field in rec.info:
                v = rec.info[args.af_field]
                af = float(v[0]) if isinstance(v, (tuple, list)) else float(v)
 
            for alt in rec.alts or []:
                alt_seq = ref_seq[:offset] + alt.upper() + ref_seq[offset + 1:]
                alt_n, alt_pos = compare(expected, alt_seq, is_rev, L)
 
                ref_central = sum(1 for x in ref_pos if x in central)
                alt_central = sum(1 for x in alt_pos if x in central)
 
                if alt_n < ref_n:
                    cls = ("ontarget_efficacy_gain"
                           if h.site_class == "on_target" else "variant_gained")
                elif alt_n > ref_n:
                    cls = ("ontarget_efficacy_loss"
                           if h.site_class == "on_target" else "variant_protected")
                elif ref_central != alt_central:
                    cls = "position_shift"
                else:
                    cls = "neutral"
 
                # ASO-relative position of the variant itself
                var_aso_pos = (L - 1 - offset) if is_rev else offset
 
                rows.append({
                    "aso_id": h.aso_id,
                    "chrom": h.chrom,
                    "site_start": h.start,
                    "site_end": h.end,
                    "strand": h.strand,
                    "site_class": h.site_class,
                    "var_pos": rec.pos,
                    "ref_allele": rec.ref,
                    "alt_allele": alt,
                    "af": af,
                    "var_aso_pos": var_aso_pos + 1,
                    "var_in_central": (var_aso_pos + 1) in central,
                    "n_mm_ref": ref_n,
                    "n_mm_alt": alt_n,
                    "delta_mm": alt_n - ref_n,
                    "mm_pos_ref": ",".join(map(str, ref_pos)) or ".",
                    "mm_pos_alt": ",".join(map(str, alt_pos)) or ".",
                    "central_ref": ref_central,
                    "central_alt": alt_central,
                    "classification": cls,
                })
 
    fa.close()
    vcf.close()
 
    if not rows:
        sys.exit("ERROR: no variant-site overlaps scored. Check that the VCF "
                 "and hits table use the same contig names and coordinates.")
 
    df = pd.DataFrame(rows).sort_values(
        ["delta_mm", "n_mm_alt", "af"], ascending=[True, True, False]
    )
    df.to_csv(args.out, sep="\t", index=False)
 
    # ---- per-candidate summary ----
    g = df.groupby("aso_id")
    summary = pd.DataFrame({
        "n_variant_sites": g.size(),
        "n_gained": g.apply(
            lambda x: (x.classification == "variant_gained").sum()),
        "n_protected": g.apply(
            lambda x: (x.classification == "variant_protected").sum()),
        "n_efficacy_loss": g.apply(
            lambda x: (x.classification == "ontarget_efficacy_loss").sum()),
        "n_position_shift": g.apply(
            lambda x: (x.classification == "position_shift").sum()),
        "worst_alt_mm": g["n_mm_alt"].min(),
        "max_af_gained": g.apply(
            lambda x: x.loc[x.classification == "variant_gained", "af"].max()
            if (x.classification == "variant_gained").any() else None,
            ),
    }).reset_index()
    summary.to_csv(args.summary, sep="\t", index=False)
 
    # ---- report ----
    print(f"\nValidation: {n_validated} sites agree with NM tag, "
          f"{n_mismatched_nm} disagree", file=sys.stderr)
    if n_mismatched_nm > n_validated * 0.01:
        print("  *** WARNING: >1% disagreement. Check genome build, strand "
              "convention, and coordinate offsets before trusting results.",
              file=sys.stderr)
    if n_indels_skipped:
        print(f"  {n_indels_skipped} indel records skipped (SNVs only)",
              file=sys.stderr)
 
    print(f"\nScored {len(df)} variant-site pairs across "
          f"{df.aso_id.nunique()} candidates", file=sys.stderr)
    for k, v in Counter(df.classification).most_common():
        print(f"  {k:28s} {v}", file=sys.stderr)
 
    gained = df[df.classification == "variant_gained"]
    if len(gained):
        print(f"\nTop variant-gained off-targets "
              f"(reference-invisible risk):", file=sys.stderr)
        cols = ["aso_id", "chrom", "var_pos", "ref_allele", "alt_allele",
                "af", "n_mm_ref", "n_mm_alt", "central_alt"]
        print(gained.head(10)[cols].to_string(index=False), file=sys.stderr)
 
        worst = gained[gained.n_mm_alt == 0]
        if len(worst):
            print(f"\n  *** {len(worst)} variants create a PERFECT "
                  f"off-target match ***", file=sys.stderr)
 
    print(f"\nWrote {args.out} and {args.summary}", file=sys.stderr)
 
 
if __name__ == "__main__":
    main()
 
