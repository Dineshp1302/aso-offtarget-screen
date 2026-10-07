#!/usr/bin/env python3
"""
Tile ASO candidates across a target transcript and apply design filters.
 
An antisense oligonucleotide is the reverse complement of its target site.
We tile TARGET SITES across the transcript, then emit the ASO (revcomp) for
each one. Both sequences are recorded so the strand convention is explicit
and auditable downstream -- getting this backwards is the single most common
silent error in ASO pipelines.
 
Usage:
    python tile_candidates.py \
        --transcripts ref/gencode.v44.transcripts.fa.gz \
        --transcript-id ENST00000237014 \
        --out-fa results/candidates.fa \
        --out-tsv results/candidates.tsv
"""
 
import argparse
import gzip
import sys
from pathlib import Path
 
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
 
 
def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]
 
 
def open_maybe_gzip(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)
 
 
def load_transcript(fasta_path, transcript_id):
    """Pull one transcript from a GENCODE transcripts FASTA.
 
    GENCODE headers are pipe-delimited:
    >ENST00000237014.8|ENSG00000118271.12|...|TTR-201|TTR|...
    We match on the ENST accession ignoring the version suffix.
    """
    want = transcript_id.split(".")[0]
    keeping, header, chunks = False, None, []
 
    with open_maybe_gzip(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if keeping:
                    break
                fields = line[1:].strip().split("|")
                if fields and fields[0].split(".")[0] == want:
                    keeping, header = True, line[1:].strip()
            elif keeping:
                chunks.append(line.strip())
 
    if not keeping:
        sys.exit(f"ERROR: transcript {transcript_id} not found in {fasta_path}")
    return header, "".join(chunks).upper()
 
 
def gc_fraction(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)
 
 
def max_homopolymer(seq):
    best = run = 1
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best
 
 
def max_base_fraction(seq):
    return max(seq.count(b) for b in "ACGT") / len(seq)
 
 
def evaluate(site, cfg):
    """Return (passes, [reasons for failure]). Filters apply to the TARGET SITE."""
    reasons = []
    if "N" in site:
        reasons.append("ambiguous_base")
    gc = gc_fraction(site)
    if gc < cfg["gc_min"]:
        reasons.append(f"gc_low({gc:.2f})")
    if gc > cfg["gc_max"]:
        reasons.append(f"gc_high({gc:.2f})")
    hp = max_homopolymer(site)
    if hp > cfg["max_homopolymer"]:
        reasons.append(f"homopolymer({hp})")
    # G-quadruplex-prone: four or more G-runs
    if site.count("GGG") >= 2 and hp >= 3:
        reasons.append("g_rich")
    if cfg["reject_cpg"] and "CG" in site:
        reasons.append("cpg")
    if max_base_fraction(site) > cfg["max_low_complexity"]:
        reasons.append("low_complexity")
    return (len(reasons) == 0), reasons
 
 
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--transcripts", required=True)
    p.add_argument("--transcript-id", required=True)
    p.add_argument("--gene", default="TGT")
    p.add_argument("--length", type=int, default=20)
    p.add_argument("--step", type=int, default=1)
    p.add_argument("--gc-min", type=float, default=0.40)
    p.add_argument("--gc-max", type=float, default=0.60)
    p.add_argument("--max-homopolymer", type=int, default=3)
    p.add_argument("--keep-cpg", action="store_true",
                   help="do not reject CpG-containing candidates")
    p.add_argument("--max-low-complexity", type=float, default=0.70)
    p.add_argument("--out-fa", required=True)
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--emit-failed", action="store_true",
                   help="also write rejected candidates to the TSV (fa stays pass-only)")
    args = p.parse_args()
 
    cfg = {
        "gc_min": args.gc_min,
        "gc_max": args.gc_max,
        "max_homopolymer": args.max_homopolymer,
        "reject_cpg": not args.keep_cpg,
        "max_low_complexity": args.max_low_complexity,
    }
 
    header, tx = load_transcript(args.transcripts, args.transcript_id)
    print(f"Loaded {args.transcript_id}: {len(tx)} nt", file=sys.stderr)
    print(f"  header: {header[:90]}...", file=sys.stderr)
 
    Path(args.out_fa).parent.mkdir(parents=True, exist_ok=True)
 
    n_total = n_pass = 0
    with open(args.out_fa, "w") as fa, open(args.out_tsv, "w") as tsv:
        tsv.write("aso_id\ttx_start\ttx_end\ttarget_site\taso_seq\t"
                  "gc\thomopolymer\tpass\treasons\n")
 
        for i in range(0, len(tx) - args.length + 1, args.step):
            site = tx[i:i + args.length]
            n_total += 1
            ok, reasons = evaluate(site, cfg)
            aso = revcomp(site)
            aso_id = f"ASO_{args.gene}_{i + 1:05d}"
 
            if ok:
                n_pass += 1
                # FASTA holds the ASO sequence. Bowtie2 reports hits on the
                # minus strand where the ASO matches a plus-strand target.
                fa.write(f">{aso_id}\n{aso}\n")
 
            if ok or args.emit_failed:
                tsv.write(
                    f"{aso_id}\t{i + 1}\t{i + args.length}\t{site}\t{aso}\t"
                    f"{gc_fraction(site):.3f}\t{max_homopolymer(site)}\t"
                    f"{'PASS' if ok else 'FAIL'}\t{';'.join(reasons) or '.'}\n"
                )
 
    print(f"Tiled {n_total} windows -> {n_pass} passed filters "
          f"({100 * n_pass / max(n_total, 1):.1f}%)", file=sys.stderr)
    print(f"Wrote {args.out_fa} and {args.out_tsv}", file=sys.stderr)
 
 
if __name__ == "__main__":
    main()
 