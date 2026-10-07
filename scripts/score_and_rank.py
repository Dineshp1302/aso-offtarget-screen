#!/usr/bin/env python3
"""
Composite risk score and final candidate ranking.
 
Design principle: every component is kept as its own column. The final score
is a transparent weighted sum, not a black box. A reviewer who disagrees with
a weight can recompute from the components without rerunning the pipeline.
 
Score components (higher = worse)
---------------------------------
1. Reference off-target burden, weighted by mismatch count and by whether
   the mismatch sits in the RNase H central window. A terminal mismatch
   leaves the duplex cleavage-competent and is therefore worse than a
   central one.
 
2. Gene context. An off-target in a protein-coding gene matters more than one
   in intergenic sequence. Genes on a user-supplied watchlist (tumour
   suppressors, essential genes) are escalated further.
 
3. Variant-driven risk. Sites that gain mismatch quality on the ALT allele,
   scaled by allele frequency: a liability at AF 0.5 affects half the
   population, one at AF 0.001 affects very few.
 
4. On-target efficacy loss. Variants inside the intended binding site reduce
   potency in carriers.
 
Candidates that span splice junctions are excluded from ranking and reported
separately: genomic alignment in end-to-end mode cannot evaluate them.
 
Usage:
    python score_and_rank.py \
        --hits results/hits_table.tsv \
        --variants results/variant_rescored.tsv \
        --candidates results/candidates.tsv \
        --annotated results/offtarget_annotated.tsv \
        --exclude results/junction_spanning_candidates.txt \
        --out results/candidates_ranked.tsv
"""
 
import argparse
import sys
 
import numpy as np
import pandas as pd
 
# Default watchlist: tumour suppressors and DNA repair genes where
# knockdown carries obvious risk. Extend via --watchlist.
DEFAULT_WATCHLIST = {
    "RBBP8", "TP53", "BRCA1", "BRCA2", "ATM", "ATR", "CHEK2", "PTEN",
    "RB1", "APC", "NF1", "NF2", "VHL", "STK11", "MLH1", "MSH2", "MSH6",
    "PMS2", "PALB2", "BARD1", "BRIP1", "RAD51C", "RAD51D", "CDKN2A",
}
 
 
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--hits", required=True)
    p.add_argument("--variants", required=False,
                   help="variant_rescored.tsv; omit to score reference only")
    p.add_argument("--candidates", required=True)
    p.add_argument("--annotated", required=False,
                   help="TSV: aso_id, chrom, start, gene_name, gene_type")
    p.add_argument("--exclude", required=False,
                   help="file of aso_ids to exclude, one per line")
    p.add_argument("--watchlist", required=False,
                   help="file of gene symbols to escalate, one per line")
    p.add_argument("--out", required=True)
 
    # weights
    p.add_argument("--w-mm0", type=float, default=100.0)
    p.add_argument("--w-mm1", type=float, default=30.0)
    p.add_argument("--w-mm2", type=float, default=8.0)
    p.add_argument("--w-mm3", type=float, default=1.0)
    p.add_argument("--central-discount", type=float, default=0.3,
                   help="multiplier applied when a mismatch is central")
    p.add_argument("--w-coding", type=float, default=2.0,
                   help="multiplier for off-targets in protein-coding genes")
    p.add_argument("--w-watchlist", type=float, default=5.0)
    p.add_argument("--w-variant-gain", type=float, default=40.0,
                   help="scaled by allele frequency")
    p.add_argument("--w-efficacy-loss", type=float, default=25.0)
    args = p.parse_args()
 
    mm_weight = {0: args.w_mm0, 1: args.w_mm1,
                 2: args.w_mm2, 3: args.w_mm3}
 
    watchlist = set(DEFAULT_WATCHLIST)
    if args.watchlist:
        with open(args.watchlist) as fh:
            watchlist |= {l.strip().upper() for l in fh if l.strip()}
 
    cand = pd.read_csv(args.candidates, sep="\t")
    cand = cand[cand["pass"] == "PASS"] if "pass" in cand.columns else cand
    hits = pd.read_csv(args.hits, sep="\t")
 
    excluded = set()
    if args.exclude:
        with open(args.exclude) as fh:
            excluded = {l.strip() for l in fh if l.strip()}
 
    # ---- 1. reference off-target burden ----
    off = hits[hits.site_class == "off_target"].copy()
    off["base_w"] = off.n_mismatch.map(mm_weight).fillna(0.0)
    # a central mismatch degrades cleavage, so the hit is less dangerous
    off["pos_factor"] = np.where(off.has_central_mismatch.astype(str) == "True",
                                 args.central_discount, 1.0)
    off["hit_score"] = off.base_w * off.pos_factor
 
    # ---- 2. gene context ----
    if args.annotated:
        ann = pd.read_csv(args.annotated, sep="\t")
        off = off.merge(ann, on=["aso_id", "chrom", "start"], how="left")
        gname = off.get("gene_name", pd.Series([np.nan] * len(off))).fillna("")
        gtype = off.get("gene_type", pd.Series([np.nan] * len(off))).fillna("")
        ctx = np.ones(len(off))
        ctx = np.where(gtype.str.contains("protein_coding"), args.w_coding, ctx)
        ctx = np.where(gname.str.upper().isin(watchlist), args.w_watchlist, ctx)
        off["ctx_factor"] = ctx
    else:
        off["ctx_factor"] = 1.0
        off["gene_name"] = ""
    off["hit_score"] = off.hit_score * off.ctx_factor
 
    agg = off.groupby("aso_id").agg(
        n_offtarget=("aso_id", "size"),
        n_mm0=("n_mismatch", lambda s: (s == 0).sum()),
        n_mm1=("n_mismatch", lambda s: (s == 1).sum()),
        n_mm2=("n_mismatch", lambda s: (s == 2).sum()),
        n_mm3=("n_mismatch", lambda s: (s == 3).sum()),
        offtarget_score=("hit_score", "sum"),
    ).reset_index()
 
    flagged = (off[off.gene_name.str.upper().isin(watchlist)]
               .groupby("aso_id").gene_name.apply(
                   lambda s: ";".join(sorted(set(s)))).rename("watchlist_genes")
               .reset_index()) if "gene_name" in off.columns else None
 
    # ---- 3 & 4. variant-driven risk ----
    if args.variants:
        var = pd.read_csv(args.variants, sep="\t")
        var["af"] = pd.to_numeric(var.af, errors="coerce").fillna(0.0)
 
        gained = var[var.classification == "variant_gained"].copy()
        # severity rises as ALT mismatch count falls
        gained["sev"] = gained.n_mm_alt.map(
            {0: 4.0, 1: 2.0, 2: 1.0, 3: 0.4}).fillna(0.2)
        gained["gain_score"] = args.w_variant_gain * gained.af * gained.sev
 
        loss = var[var.classification == "ontarget_efficacy_loss"].copy()
        loss["loss_score"] = args.w_efficacy_loss * loss.af
 
        vagg = gained.groupby("aso_id").agg(
            n_variant_gained=("aso_id", "size"),
            max_af_gained=("af", "max"),
            best_alt_mm=("n_mm_alt", "min"),
            variant_gain_score=("gain_score", "sum"),
        ).reset_index()
 
        lagg = loss.groupby("aso_id").agg(
            n_efficacy_loss=("aso_id", "size"),
            efficacy_loss_score=("loss_score", "sum"),
        ).reset_index()
    else:
        vagg = pd.DataFrame(columns=["aso_id", "n_variant_gained",
                                     "max_af_gained", "best_alt_mm",
                                     "variant_gain_score"])
        lagg = pd.DataFrame(columns=["aso_id", "n_efficacy_loss",
                                     "efficacy_loss_score"])
 
    # ---- assemble ----
    out = cand[["aso_id", "tx_start", "tx_end", "target_site",
                "aso_seq", "gc"]].copy()
    out = out.merge(agg, on="aso_id", how="left")
    if flagged is not None:
        out = out.merge(flagged, on="aso_id", how="left")
    out = out.merge(vagg, on="aso_id", how="left")
    out = out.merge(lagg, on="aso_id", how="left")
 
    num_cols = ["n_offtarget", "n_mm0", "n_mm1", "n_mm2", "n_mm3",
                "offtarget_score", "n_variant_gained", "variant_gain_score",
                "n_efficacy_loss", "efficacy_loss_score"]
    for c in num_cols:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0)
 
    out["total_risk_score"] = (out.offtarget_score
                               + out.variant_gain_score
                               + out.efficacy_loss_score)
 
    out["excluded"] = out.aso_id.isin(excluded)
    out["exclude_reason"] = np.where(out.excluded,
                                     "splice_junction_spanning", "")
 
    rankable = out[~out.excluded].copy()
    rankable = rankable.sort_values("total_risk_score")
    rankable["rank"] = range(1, len(rankable) + 1)
    out = pd.concat([rankable, out[out.excluded]], ignore_index=True)
 
    out.to_csv(args.out, sep="\t", index=False)
 
    # ---- report ----
    print(f"Scored {len(out)} candidates "
          f"({out.excluded.sum()} excluded as junction-spanning)",
          file=sys.stderr)
    print(f"\nTop 10 lowest-risk candidates:", file=sys.stderr)
    cols = ["rank", "aso_id", "tx_start", "n_offtarget", "n_mm1", "n_mm2",
            "n_variant_gained", "total_risk_score"]
    cols = [c for c in cols if c in rankable.columns]
    print(rankable.head(10)[cols].to_string(index=False), file=sys.stderr)
 
    print(f"\nTop 5 highest-risk candidates:", file=sys.stderr)
    print(rankable.tail(5)[cols].to_string(index=False), file=sys.stderr)
 
    if "watchlist_genes" in out.columns:
        wl = out[out.watchlist_genes.notna() & (out.watchlist_genes != "")]
        if len(wl):
            print(f"\n{len(wl)} candidates hit watchlist genes:",
                  file=sys.stderr)
            print(wl[["aso_id", "watchlist_genes", "total_risk_score"]]
                  .to_string(index=False), file=sys.stderr)
 
    print(f"\nWrote {args.out}", file=sys.stderr)
 
 
if __name__ == "__main__":
    main()
