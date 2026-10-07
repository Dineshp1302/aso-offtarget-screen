#!/usr/bin/env python3
"""
Generate figures for the ASO off-target screen.
 
Five figures, each answering one question:
 
  fig1  Design filters   -- why did 57% of tiled windows fail?
  fig2  Off-target burden -- how do candidates differ, and at what severity?
  fig3  Mismatch position -- does position (not just count) separate hits?
  fig4  Variant effect    -- how often do variants create vs destroy off-targets?
  fig5  Final ranking     -- composite score across the transcript
 
Usage:
    python make_figures.py \
        --candidates results/candidates.tsv \
        --hits results/hits_table.tsv \
        --variants results/variant_rescored.tsv \
        --ranked results/candidates_ranked.tsv \
        --outdir results/figures
"""
 
import argparse
import os
from collections import Counter
 
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
 
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
})
 
C = {"mm0": "#8B0000", "mm1": "#D7301F", "mm2": "#FC8D59",
     "mm3": "#FDCC8A", "gain": "#B2182B", "prot": "#2166AC",
     "neut": "#999999", "ok": "#2C7FB8"}
 
 
def fig1_filters(cand, outdir):
    """What the design filters removed, and why."""
    if "pass" not in cand.columns:
        return
    reasons = Counter()
    for r in cand.loc[cand["pass"] == "FAIL", "reasons"].dropna():
        for item in str(r).split(";"):
            key = item.split("(")[0].strip()
            if key and key != ".":
                reasons[key] += 1
 
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.4),
                                 gridspec_kw={"width_ratios": [1, 1.6]})
 
    n_pass = (cand["pass"] == "PASS").sum()
    n_fail = (cand["pass"] == "FAIL").sum()
    a1.bar(["pass", "fail"], [n_pass, n_fail],
           color=[C["ok"], "#CCCCCC"], width=0.6)
    for i, v in enumerate([n_pass, n_fail]):
        a1.text(i, v, f"{v}\n({100*v/(n_pass+n_fail):.0f}%)",
                ha="center", va="bottom", fontsize=8)
    a1.set_ylabel("tiled windows")
    a1.set_title("Design filter outcome", fontsize=10, loc="left")
    a1.set_ylim(0, max(n_pass, n_fail) * 1.25)
 
    if reasons:
        items = reasons.most_common()
        labels = [k for k, _ in items]
        vals = [v for _, v in items]
        y = np.arange(len(labels))
        a2.barh(y, vals, color="#999999", height=0.65)
        a2.set_yticks(y)
        a2.set_yticklabels(labels)
        a2.invert_yaxis()
        a2.set_xlabel("windows rejected (reasons overlap)")
        a2.set_title("Rejection reasons", fontsize=10, loc="left")
        for i, v in enumerate(vals):
            a2.text(v, i, f" {v}", va="center", fontsize=8)
 
    fig.suptitle("Figure 1. ASO candidate design filtering",
                 fontsize=11, x=0.01, ha="left", weight="bold")
    fig.savefig(os.path.join(outdir, "fig1_design_filters.png"))
    plt.close(fig)
 
 
def fig2_burden(hits, outdir, top_n=40):
    """Off-target burden per candidate, stacked by mismatch count."""
    off = hits[hits.site_class == "off_target"]
    if off.empty:
        return
    tab = (off.groupby(["aso_id", "n_mismatch"]).size()
           .unstack(fill_value=0))
    for c in (0, 1, 2, 3):
        if c not in tab.columns:
            tab[c] = 0
    tab = tab[[0, 1, 2, 3]]
    # order by severity: perfect hits first, then 1mm, then total
    tab = tab.sort_values([0, 1, 2, 3], ascending=False)
    show = tab.head(top_n)
 
    fig, ax = plt.subplots(figsize=(10, 3.8))
    bottom = np.zeros(len(show))
    for c, col in zip([3, 2, 1, 0], [C["mm3"], C["mm2"], C["mm1"], C["mm0"]]):
        ax.bar(range(len(show)), show[c], bottom=bottom, color=col,
               label=f"{c} mismatch", width=0.85)
        bottom += show[c].values
 
    ax.set_xticks(range(len(show)))
    ax.set_xticklabels([s.replace("ASO_TTR_", "") for s in show.index],
                       rotation=90, fontsize=6)
    ax.set_xlabel("candidate (transcript position)")
    ax.set_ylabel("off-target sites")
    ax.legend(frameon=False, fontsize=8, ncol=4, loc="upper right")
    ax.set_title(f"Figure 2. Off-target burden, {top_n} worst candidates",
                 fontsize=11, loc="left", weight="bold")
    fig.savefig(os.path.join(outdir, "fig2_offtarget_burden.png"))
    plt.close(fig)
 
 
def fig3_positions(hits, outdir, central=(8, 12)):
    """Mismatch position distribution and its effect on risk tier."""
    off = hits[hits.site_class == "off_target"].copy()
    if off.empty:
        return
 
    pos_by_mm = {}
    for mm in (1, 2):
        sub = off[off.n_mismatch == mm]
        counts = Counter()
        for s in sub.mismatch_pos.dropna():
            if s == ".":
                continue
            for p in str(s).split(","):
                counts[int(p)] += 1
        pos_by_mm[mm] = counts
 
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.4))
 
    width = 0.4
    xs = np.arange(1, 21)
    for i, (mm, col) in enumerate(zip((1, 2), [C["mm1"], C["mm2"]])):
        vals = [pos_by_mm[mm].get(p, 0) for p in xs]
        a1.bar(xs + (i - 0.5) * width, vals, width=width,
               color=col, label=f"{mm} mismatch")
    a1.axvspan(central[0] - 0.5, central[1] + 0.5, color="#2166AC", alpha=0.12)
    a1.text(np.mean(central), a1.get_ylim()[1] * 0.92, "RNase H\ncentral window",
            ha="center", fontsize=7, color="#2166AC")
    a1.set_xticks(xs)
    a1.set_xticklabels(xs, fontsize=6)
    a1.set_xlabel("mismatch position in ASO (5'→3')")
    a1.set_ylabel("off-target hits")
    a1.legend(frameon=False, fontsize=8)
    a1.set_title("Where mismatches fall", fontsize=10, loc="left")
 
    grp = (off.groupby(["n_mismatch", "has_central_mismatch"])
           .size().unstack(fill_value=0))
    idx = np.arange(len(grp))
    for i, (col_key, col, lbl) in enumerate([
            (False, "#D7301F", "terminal only (cleavage competent)"),
            (True, "#2166AC", "central (cleavage impaired)")]):
        if col_key in grp.columns:
            a2.bar(idx + (i - 0.5) * 0.38, grp[col_key], width=0.38,
                   color=col, label=lbl)
    a2.set_xticks(idx)
    a2.set_xticklabels([f"{int(i)} mm" for i in grp.index])
    a2.set_yscale("log")
    a2.set_ylabel("off-target hits (log)")
    a2.legend(frameon=False, fontsize=7)
    a2.set_title("Cleavage competence by tier", fontsize=10, loc="left")
 
    fig.suptitle("Figure 3. Mismatch position, not just count, determines risk",
                 fontsize=11, x=0.01, ha="left", weight="bold")
    fig.savefig(os.path.join(outdir, "fig3_mismatch_positions.png"))
    plt.close(fig)
 
 
def fig4_variants(var, outdir):
    """Variant effect on binding-site quality."""
    if var is None or var.empty:
        return
    var = var.copy()
    var["af"] = pd.to_numeric(var.af, errors="coerce").fillna(0)
 
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 3.4),
                                 gridspec_kw={"width_ratios": [1, 1.4]})
 
    counts = var.classification.value_counts()
    cols = {"variant_protected": C["prot"], "variant_gained": C["gain"],
            "neutral": C["neut"], "position_shift": "#F4A582",
            "ontarget_efficacy_loss": "#762A83",
            "ontarget_efficacy_gain": "#5AAE61"}
    a1.barh(range(len(counts)), counts.values,
            color=[cols.get(k, "#999") for k in counts.index], height=0.65)
    a1.set_yticks(range(len(counts)))
    a1.set_yticklabels([k.replace("_", " ") for k in counts.index], fontsize=8)
    a1.invert_yaxis()
    a1.set_xlabel("variant–site pairs")
    a1.set_title("Variant effect on duplex quality", fontsize=10, loc="left")
    for i, v in enumerate(counts.values):
        a1.text(v, i, f" {v}", va="center", fontsize=8)
 
    gained = var[var.classification == "variant_gained"]
    other = var[var.classification != "variant_gained"]
    a2.scatter(other.af, other.n_mm_alt + np.random.uniform(-.12, .12, len(other)),
               s=10, c=C["neut"], alpha=0.35, label="other", edgecolors="none")
    a2.scatter(gained.af, gained.n_mm_alt + np.random.uniform(-.12, .12, len(gained)),
               s=28, c=C["gain"], alpha=0.85, label="variant gained",
               edgecolors="white", linewidths=0.4)
    a2.set_xscale("log")
    a2.set_xlabel("allele frequency (gnomAD)")
    a2.set_ylabel("mismatches on ALT allele")
    a2.invert_yaxis()
    a2.legend(frameon=False, fontsize=8, loc="lower left")
    a2.set_title("Risk vs population frequency", fontsize=10, loc="left")
 
    if len(gained):
        top = gained.nlargest(1, "af").iloc[0]
        a2.annotate(f"{top.aso_id.replace('ASO_TTR_','')}\nAF={top.af:.2f}",
                    xy=(top.af, top.n_mm_alt), xytext=(12, -18),
                    textcoords="offset points", fontsize=7,
                    arrowprops=dict(arrowstyle="->", lw=0.6, color="#555"))
 
    fig.suptitle("Figure 4. Variant-aware re-scoring of binding sites",
                 fontsize=11, x=0.01, ha="left", weight="bold")
    fig.savefig(os.path.join(outdir, "fig4_variant_effect.png"))
    plt.close(fig)
 
 
def fig5_ranking(ranked, outdir):
    """Composite risk across the transcript."""
    if ranked is None or ranked.empty:
        return
    r = ranked[~ranked.excluded.astype(bool)].copy() \
        if "excluded" in ranked.columns else ranked.copy()
    r["total_risk_score"] = pd.to_numeric(r.total_risk_score, errors="coerce")
    r = r.dropna(subset=["total_risk_score", "tx_start"])
 
    fig, ax = plt.subplots(figsize=(10, 3.4))
 
    has_wl = (r.watchlist_genes.notna() & (r.watchlist_genes != "")) \
        if "watchlist_genes" in r.columns else pd.Series(False, index=r.index)
 
    ax.scatter(r.loc[~has_wl, "tx_start"], r.loc[~has_wl, "total_risk_score"],
               s=18, c=C["ok"], alpha=0.7, edgecolors="none", label="candidate")
    if has_wl.any():
        ax.scatter(r.loc[has_wl, "tx_start"], r.loc[has_wl, "total_risk_score"],
                   s=42, c=C["gain"], alpha=0.9, edgecolors="white",
                   linewidths=0.5, label="hits watchlist gene")
 
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("position in TTR transcript (nt)")
    ax.set_ylabel("composite risk score (log)")
    ax.legend(frameon=False, fontsize=8)
 
    best = r.nsmallest(3, "total_risk_score")
    for _, b in best.iterrows():
        ax.annotate(b.aso_id.replace("ASO_TTR_", ""),
                    xy=(b.tx_start, b.total_risk_score), xytext=(0, 10),
                    textcoords="offset points", fontsize=7, ha="center",
                    color="#2C7FB8")
    worst = r.nlargest(3, "total_risk_score")
    for _, b in worst.iterrows():
        ax.annotate(b.aso_id.replace("ASO_TTR_", ""),
                    xy=(b.tx_start, b.total_risk_score), xytext=(0, 8),
                    textcoords="offset points", fontsize=7, ha="center",
                    color=C["gain"])
 
    ax.set_title("Figure 5. Composite risk score across the target transcript",
                 fontsize=11, loc="left", weight="bold")
    fig.savefig(os.path.join(outdir, "fig5_ranking.png"))
    plt.close(fig)
 
 
def main():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--hits", required=True)
    p.add_argument("--variants")
    p.add_argument("--ranked")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()
 
    os.makedirs(args.outdir, exist_ok=True)
 
    cand = pd.read_csv(args.candidates, sep="\t")
    hits = pd.read_csv(args.hits, sep="\t")
    var = pd.read_csv(args.variants, sep="\t") if args.variants else None
    ranked = pd.read_csv(args.ranked, sep="\t") if args.ranked else None
 
    fig1_filters(cand, args.outdir)
    fig2_burden(hits, args.outdir)
    fig3_positions(hits, args.outdir)
    fig4_variants(var, args.outdir)
    fig5_ranking(ranked, args.outdir)
 
    print(f"Figures written to {args.outdir}:")
    for f in sorted(os.listdir(args.outdir)):
        if f.endswith(".png"):
            print(f"  {f}")
 
 
if __name__ == "__main__":
    main()
 
