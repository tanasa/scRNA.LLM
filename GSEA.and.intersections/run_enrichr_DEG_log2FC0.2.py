#!/usr/bin/env python3
"""
Enrichr (gseapy) + barplots in ONE script.
  - DEG mode (default): tables/DEG_log2FC0.2/genes_only/*.txt
      -> tables/DEG_log2FC0.2/enrichr_results/*_enrichr.csv + barplots/
  - genes_in_common mode: genes_in_common/*_genes_in_common.txt
      -> genes_in_common/enrichr_results/*_enrichr.csv + barplots/

Usage:
  python3 run_enrichr_DEG_log2FC0.2.py              # DEG log2FC0.2
  python3 run_enrichr_DEG_log2FC0.2.py --genes-in-common
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ENRICHR_LIBRARIES = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
    "WikiPathway_2023_Human",
    "GO_Biological_Process_2023",
]

ADJ_PVAL = 0.2
TOP_N = 25
LIBRARIES = [
    ("WikiPathway_2023_Human", "WikiPathways"),
    ("MSigDB_Hallmark_2020", "Hallmark"),
    ("GO_Biological_Process_2023", "GO_BP"),
    ("Reactome_2022", "Reactome"),
    ("KEGG_2021_Human", "KEGG"),
]


def plot_enrichr_csv(csv_path: Path, barplots_dir: Path) -> None:
    """Barplots for one combined Enrichr CSV (Gene_set + Term + Adjusted P-value)."""
    barplots_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(csv_path)
    pcol = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
    if pcol not in df.columns or "Gene_set" not in df.columns or "Term" not in df.columns:
        print("  Skip plots:", csv_path.name, "(missing columns)")
        return
    base_name = csv_path.stem.replace("_enrichr", "")
    for gene_set_value, short_name in LIBRARIES:
        sub = df[(df["Gene_set"] == gene_set_value) & (df[pcol] < ADJ_PVAL)].copy()
        if sub.empty:
            continue
        sub = sub.sort_values(pcol).head(TOP_N)
        pvals = np.clip(sub[pcol].values, 1e-300, 1.0)
        y = -np.log10(pvals)
        terms = sub["Term"].tolist()
        terms_short = [t[:60] + "..." if len(t) > 60 else t for t in terms]
        n = len(terms_short)
        fig, ax = plt.subplots(figsize=(8, max(5, n * 0.35)))
        y_pos = np.arange(n)[::-1]
        ax.barh(y_pos, y, color="steelblue", edgecolor="navy", linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(terms_short, fontsize=8)
        ax.set_xlabel("-log10(Adjusted P-value)")
        ax.set_title(f"{base_name} — {short_name}\n(Adj P < {ADJ_PVAL}, top {n} terms)")
        ax.axvline(-np.log10(ADJ_PVAL), color="gray", linestyle="--", linewidth=1, alpha=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        plt.tight_layout()
        plt.savefig(barplots_dir / f"{base_name}_{short_name}_barplot.png", dpi=150, bbox_inches="tight")
        plt.close()


def run_enrichr_and_plot(txt_path: Path, out_dir: Path) -> None:
    try:
        import gseapy as gp
    except ImportError:
        print("ERROR: pip install gseapy")
        sys.exit(1)
    genes = [g.strip() for g in txt_path.read_text().strip().splitlines() if g.strip()]
    if not genes:
        print(txt_path.name, ": empty, skip")
        return
    print(txt_path.name, f"({len(genes)} genes)...", end=" ", flush=True)
    all_dfs = []
    for lib in ENRICHR_LIBRARIES:
        enr = gp.enrichr(genes, gene_sets=lib, organism="Human", outdir=None, no_plot=True, verbose=False)
        res = enr.results.copy()
        res["Gene_set"] = lib
        all_dfs.append(res)
    combined = pd.concat(all_dfs, ignore_index=True)
    out_csv = out_dir / f"{txt_path.stem}_enrichr.csv"
    combined.to_csv(out_csv, index=False)
    print("->", out_csv.name)
    # Barplots next to tables
    barplots_dir = out_dir / "barplots"
    plot_enrichr_csv(out_csv, barplots_dir)


def main():
    genes_in_common = "--genes-in-common" in sys.argv
    base = Path(__file__).resolve().parent.parent

    if genes_in_common:
        genes_dir = base / "genes_in_common"
        out_dir = genes_dir / "enrichr_results"
        pattern = "*_genes_in_common.txt"
    else:
        genes_dir = base / "tables" / "DEG_log2FC0.2" / "genes_only"
        out_dir = base / "tables" / "DEG_log2FC0.2" / "enrichr_results"

    if not genes_dir.is_dir():
        print("ERROR: Not found:", genes_dir)
        sys.exit(1)
    out_dir.mkdir(parents=True, exist_ok=True)

    if genes_in_common:
        txt_files = sorted(genes_dir.glob(pattern))
    else:
        txt_files = sorted(genes_dir.glob("*.txt"))

    if not txt_files:
        print("No input files in", genes_dir)
        sys.exit(1)

    print("Mode:", "genes_in_common" if genes_in_common else "DEG_log2FC0.2")
    print("Input:", genes_dir)
    print("Output:", out_dir)
    print("Barplots:", out_dir / "barplots", "\n")

    for path in txt_files:
        run_enrichr_and_plot(path, out_dir)

    print("\nDone. Tables + barplots in:", out_dir)


if __name__ == "__main__":
    main()
