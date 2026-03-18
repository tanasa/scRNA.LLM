#!/usr/bin/env python3
# =============================================================================
# FGSEA (prerank GSEA) using DESeq2 CSV results + GMT files in databases_GSEA
# =============================================================================
# - Reads DESeq2 CSVs (Gene_symbol, log2FoldChange, padj, ...) from project root
# - Ranks genes by log2FoldChange
# - Loads GMTs from databases_GSEA (Hallmark, KEGG, Reactome, WikiPathways, GO_BP)
# - Runs gseapy.prerank per cell type x pathway set; saves tables and plots
# Requires: pandas, gseapy, matplotlib
# =============================================================================

import os
import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import gseapy as gp
except ImportError:
    print("ERROR: gseapy not installed. Install with: pip install gseapy")
    sys.exit(1)

# ----- Paths (run from scripts_FGSEA or project root) -----
if Path("databases_GSEA").is_dir():
    base_dir = Path.cwd()
else:
    base_dir = Path.cwd().parent
csv_dir = base_dir
gmt_dir = base_dir / "databases_GSEA"
out_dir = base_dir / "GSEA_results_DESeq2_python"
out_dir.mkdir(exist_ok=True)
(out_dir / "tables").mkdir(exist_ok=True)
(out_dir / "plots").mkdir(exist_ok=True)

# ----- DESeq2 CSV files -----
csv_files = sorted(csv_dir.glob("*.csv"))
csv_files = [f for f in csv_files if f.name != "summary_DEG_counts.csv"]
if not csv_files:
    print("ERROR: No CSV files found in", csv_dir)
    sys.exit(1)
print("DESeq2 CSV files to process:", len(csv_files))
print("GMT directory:", gmt_dir, "\n")

# ----- GMT files in databases_GSEA -----
gmt_sets = {
    "Hallmark":     gmt_dir / "h.all.v2026.1.Hs.symbols.gmt",
    "KEGG":         gmt_dir / "c2.cp.kegg_medicus.v2026.1.Hs.symbols.gmt",
    "Reactome":     gmt_dir / "c2.cp.reactome.v2026.1.Hs.symbols.gmt",
    "WikiPathways": gmt_dir / "c2.cp.wikipathways.v2026.1.Hs.symbols.gmt",
    "GO_BP":        gmt_dir / "c5.go.bp.v2026.1.Hs.symbols.gmt",
}
gmt_sets = {k: v for k, v in gmt_sets.items() if v.is_file()}
if not gmt_sets:
    print("ERROR: No GMT files found in", gmt_dir)
    sys.exit(1)
print("GMT sets to use:", ", ".join(gmt_sets.keys()), "\n")

# ----- Process each DESeq2 CSV -----
for csv_path in csv_files:
    cell_type = csv_path.stem
    print("======== ", cell_type, " ========", sep="")

    res = pd.read_csv(csv_path)
    if "Gene_symbol" not in res.columns:
        print("  Skip: no Gene_symbol column.\n")
        continue
    if "log2FoldChange" not in res.columns:
        print("  Skip: no log2FoldChange column.\n")
        continue

    # Rank by log2FoldChange; drop NA, dedupe by mean; uppercase symbols for GMT match
    res2 = res[["Gene_symbol", "log2FoldChange"]].dropna()
    res2["Gene_symbol"] = res2["Gene_symbol"].astype(str).str.upper()
    rnk = (
        res2.groupby("Gene_symbol", as_index=True)["log2FoldChange"]
        .mean()
        .sort_values(ascending=False)
    )
    print("  Ranked genes:", len(rnk))

    for set_name, gmt_file in gmt_sets.items():
        print("  Running prerank:", set_name, "... ", end="", flush=True)
        try:
            pr = gp.prerank(
                rnk=rnk,
                gene_sets=str(gmt_file),
                outdir=None,
                permutation_num=1000,
                no_plot=True,
                seed=42,
                min_size=10,
                max_size=500,
            )
            # gseapy .res2d: pathway names are in "Term". A column sometimes named "pathway"
            # can wrongly contain the method label "prerank" for every row — use Term for labels.
            df = pr.res2d.copy()
            if "Term" in df.columns:
                df["pathway"] = df["Term"].astype(str)
                df = df.drop(columns=["Term"], errors="ignore")
            elif "Name" in df.columns:
                df = df.rename(columns={"Name": "pathway"})
            # p-value columns (gseapy may use NOM p-val / FDR q-val)
            if "NOM p-val" in df.columns and "pval" not in df.columns:
                df = df.rename(columns={"NOM p-val": "pval"})
            if "FDR q-val" in df.columns and "padj" not in df.columns:
                df = df.rename(columns={"FDR q-val": "padj"})
            df = df.sort_values("NES", ascending=False)

            out_table = out_dir / "tables" / f"{cell_type}_{set_name}_fgsea.csv"
            df.to_csv(out_table, index=False)

            # Barplot (top pathways by NES)
            n_plot = min(30, len(df))
            plot_df = df.head(n_plot)
            padj_col = "padj" if "padj" in plot_df.columns else "pval"
            fig, ax = plt.subplots(figsize=(10, max(6, n_plot * 0.2)))
            y_pos = range(len(plot_df))
            nes = plot_df["NES"].values
            colors = ["#e41a1c" if p < 0.05 else "#377eb8" for p in plot_df[padj_col]]
            ax.barh(y_pos, nes, color=colors)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(plot_df["pathway"].tolist(), fontsize=8)
            ax.set_xlabel("NES")
            ax.set_ylabel("")
            ax.set_title(f"{cell_type} - {set_name}")
            ax.invert_yaxis()
            plt.tight_layout()
            plt.savefig(
                out_dir / "plots" / f"{cell_type}_{set_name}_barplot.png",
                dpi=150,
                bbox_inches="tight",
            )
            plt.close()
            print("OK")
        except Exception as e:
            print("ERROR:", e)

    print()

print("Done. Results in:", out_dir)
print("  tables/: *_fgsea.csv")
print("  plots/:  *_barplot.png")
