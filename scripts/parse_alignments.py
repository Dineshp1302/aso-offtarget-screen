"""
Parse Bowtie2 alignments into a tidy per-hit table.
 
For each alignment of each ASO candidate this records:
  - genomic location and strand
  - mismatch count (NM tag)
  - mismatch positions EXPRESSED IN ASO COORDINATES (1-based, 5'->3')
  - whether any mismatch falls in the RNase H central window
  - on-target vs off-target classification
 
Why ASO-relative mismatch positions matter
------------------------------------------
RNase H1 cleaves the RNA strand of a DNA/RNA heteroduplex. The enzyme needs a
contiguous, well-formed duplex across its footprint, which sits over the
central region of the oligo. A mismatch at position 1 or 20 barely perturbs
cleavage -- the duplex still forms and the terminus frays anyway. A mismatch at
position 10 disrupts the footprint and can abolish cleavage entirely.
 
So a 2-mismatch off-target with both mismatches central is relatively safe,
while a 2-mismatch off-target with both at the termini is dangerous. Mismatch
COUNT alone cannot distinguish these; mismatch POSITION can. This is the main
reason to parse alignments properly instead of just counting NM tags.
 
The strand correction is the subtle part. When an alignment is on the minus
strand, the SEQ field in the SAM record holds the reverse complement of the
query as submitted. pysam reports query positions against that stored,
flipped sequence. To recover the position in the original ASO (5'->3' as
synthesised) we must flip back:  aso_pos = read_len - 1 - stored_pos.
Getting this wrong mirrors every mismatch position around the oligo centre,
which is nearly invisible in aggregate statistics but wrong in every detail.
 
Requires the MD tag, which Bowtie2 emits by default.
 
Usage:
    python parse_alignments.py \
        --bam results/hits.bam \
        --candidates results/candidates.tsv \
        --target-region chr18:31591877-31598821 \
        --out-tsv results/hits_table.tsv \
        --out-bed results/offtarget_sites.bed
"""
 
import argparse
import sys
from collections import Counter
 
import pandas as pd
import pysam
 
 
def parse_region(region):
    """'chr18:31591877-31598821' -> ('chr18', 31591877, 31598821)"""
    try:
        chrom, span = region.split(":")
        start, end = span.replace(",", "").split("-")
        return chrom, int(start), int(end)
    except ValueError:
        sys.exit(f"ERROR: could not parse --target-region '{region}'. "
                 "Expected format chr18:31591877-31598821")
 
 
def mismatch_positions_aso(read):
    """Mismatch positions in ASO coordinates, 1-based from the 5' end.
 
    get_aligned_pairs(with_seq=True) yields (query_pos, ref_pos, ref_base)
    where ref_base is lowercase at mismatched positions. query_pos indexes the
    SEQ field as stored, which is reverse-complemented for minus-strand hits.
    """
    length = read.query_length or len(read.query_sequence or "")
    out = []
 
    for qpos, rpos, rbase in read.get_aligned_pairs(with_seq=True):
        if qpos is None or rpos is None or rbase is None:
            continue              # indel position; we forbid gaps anyway
        if rbase.islower():       # lowercase ref base == mismatch
            aso_pos = (length - 1 - qpos) if read.is_reverse else qpos
            out.append(aso_pos + 1)
 
    return sorted(out)
 
 
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--bam", required=True)
    p.add_argument("--candidates", required=True,
                   help="candidates.tsv from tile_candidates.py")
    p.add_argument("--target-region", required=True,
                   help="intended locus, e.g. chr18:31591877-31598821")
    p.add_argument("--central-start", type=int, default=8,
                   help="1-based start of RNase H central window")
    p.add_argument("--central-end", type=int, default=12)
    p.add_argument("--max-mismatch", type=int, default=3)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-bed", required=True,
                   help="BED of OFF-target sites only, for variant intersection")
    p.add_argument("--bed-flank", type=int, default=0,
                   help="extend off-target intervals by N bp on each side")
    args = p.parse_args()
 
    tchrom, tstart, tend = parse_region(args.target_region)
    central = set(range(args.central_start, args.central_end + 1))
 
    cand = pd.read_csv(args.candidates, sep="\t")
    cand_pass = cand[cand["pass"] == "PASS"] if "pass" in cand.columns else cand
    aso_seqs = dict(zip(cand_pass["aso_id"], cand_pass["aso_seq"]))
    tx_start = dict(zip(cand_pass["aso_id"], cand_pass["tx_start"]))
 
    rows = []
    skipped_gapped = 0
    skipped_overmm = 0
 
    bam = pysam.AlignmentFile(args.bam, "rb")
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
 
        # Gaps are chemically implausible for a 20-mer heteroduplex and we
        # penalised them heavily in bowtie2; drop any that slipped through.
        cig = read.cigartuples or []
        if any(op in (1, 2, 3) for op, _ in cig):
            skipped_gapped += 1
            continue
 
        nm = read.get_tag("NM") if read.has_tag("NM") else None
        if nm is None or nm > args.max_mismatch:
            skipped_overmm += 1
            continue
 
        try:
            mm_pos = mismatch_positions_aso(read)
        except ValueError:
            # MD tag missing
            sys.exit("ERROR: alignments lack MD tags. Re-run bowtie2 (it emits "
                     "MD by default) or add them with "
                     "`samtools calmd -b hits.bam ref.fa`.")
 
        chrom = read.reference_name
        start = read.reference_start                 # 0-based
        end = start + read.query_length
        strand = "-" if read.is_reverse else "+"
 
        on_target = (chrom == tchrom and start < tend and end > tstart)
        n_central = sum(1 for x in mm_pos if x in central)
 
        rows.append({
            "aso_id": read.query_name,
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "n_mismatch": nm,
            "mismatch_pos": ",".join(map(str, mm_pos)) or ".",
            "n_central_mismatch": n_central,
            "has_central_mismatch": n_central > 0,
            "site_class": "on_target" if on_target else "off_target",
            "aln_score": read.get_tag("AS") if read.has_tag("AS") else None,
            "aso_seq": aso_seqs.get(read.query_name, ""),
            "tx_start": tx_start.get(read.query_name, pd.NA),
        })
    bam.close()
 
    if not rows:
        sys.exit("ERROR: no usable alignments found. Check the BAM is sorted, "
                 "indexed, and non-empty.")
 
    df = pd.DataFrame(rows).sort_values(
        ["aso_id", "n_mismatch", "chrom", "start"]
    )
    df.to_csv(args.out_tsv, sep="\t", index=False)
 
    # BED of off-target sites only -- this is what bcftools/bedtools consume.
    off = df[df.site_class == "off_target"].copy()
    if args.bed_flank:
        off["start"] = (off["start"] - args.bed_flank).clip(lower=0)
        off["end"] = off["end"] + args.bed_flank
 
    bed = off[["chrom", "start", "end", "aso_id", "n_mismatch", "strand"]]
    bed = bed.sort_values(["chrom", "start"])
    bed.to_csv(args.out_bed, sep="\t", index=False, header=False)
 
    # ---- summary ----
    n_aso = df.aso_id.nunique()
    on = df[df.site_class == "on_target"]
    print(f"Parsed {len(df)} alignments for {n_aso} candidates", file=sys.stderr)
    if skipped_gapped:
        print(f"  skipped {skipped_gapped} gapped alignments", file=sys.stderr)
    if skipped_overmm:
        print(f"  skipped {skipped_overmm} alignments over "
              f"{args.max_mismatch} mismatches", file=sys.stderr)
    print(f"  on-target  : {len(on)}", file=sys.stderr)
    print(f"  off-target : {len(off)}", file=sys.stderr)
 
    print("\n  off-target hits by mismatch count:", file=sys.stderr)
    for k, v in sorted(Counter(off.n_mismatch).items()):
        tag = "  <-- HIGH RISK" if k <= 1 else ""
        print(f"    {k} mismatch: {v}{tag}", file=sys.stderr)
 
    perfect = set(on[on.n_mismatch == 0].aso_id)
    missing = sorted(set(df.aso_id) - perfect)
    print(f"\n  CONTROL: {len(perfect)}/{n_aso} candidates have a perfect "
          f"on-target hit", file=sys.stderr)
    if missing:
        print(f"  {len(missing)} without one (expect splice-junction spanners): "
              f"{', '.join(missing[:5])}"
              f"{' ...' if len(missing) > 5 else ''}", file=sys.stderr)
 
    print(f"\nWrote {args.out_tsv} and {args.out_bed}", file=sys.stderr)
 
 
if __name__ == "__main__":
    main()
 
