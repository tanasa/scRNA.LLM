#!/usr/bin/env python
# coding: utf-8

# In[5]:


print(
    
"""
Script to process h5ad file for Geneformer compatibility
Updated for: SCANVI_INTEGRATED_FILTERED.h5ad
Location: /mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026
"""

)

# reuse the old script : process_scanvi_h5ad_for_geneformer.py


# In[6]:


import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================

print("=" * 70)
print("Setting up environment...")
print("=" * 70)

# Set up single-cell environment
env_path = "/home/btanasa/miniconda3/envs/torch-cu"
os.environ["CONDA_PREFIX"] = env_path
os.environ["CONDA_DEFAULT_ENV"] = "torch-cu"
os.environ["RETICULATE_PYTHON"] = f"{env_path}/bin/python"
print("✅ Environment set to single-cell")
print(f"CONDA_PREFIX: {os.environ.get('CONDA_PREFIX')}")
print(f"CONDA_DEFAULT_ENV: {os.environ.get('CONDA_DEFAULT_ENV')}")
print(f"Python executable: {sys.executable}")


# In[7]:


# =============================================================================
# WORKING DIRECTORY AND FILE PATHS
# =============================================================================

# Set working directory
workdir = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026"
os.chdir(workdir)
print(f"\n📁 Working directory: {os.getcwd()}")

# List files in directory
print("\n📄 Files in directory:")
for file in os.listdir(workdir):
    if os.path.isfile(os.path.join(workdir, file)):
        print(f"  - {file}")

# h5ad file
h5ad_file = "SCANVI_INTEGRATED_FILTERED.h5ad"
h5ad_path = os.path.join(workdir, h5ad_file)

if not os.path.exists(h5ad_path):
    raise FileNotFoundError(f"❌ h5ad file not found: {h5ad_path}")

print(f"\n✅ Found h5ad file: {h5ad_file}")

# Output prefix
prefix = Path(h5ad_file).stem
print(f"📝 Output prefix: {prefix}")


# In[8]:


# =============================================================================
# LOAD DATA
# =============================================================================

print("\n" + "=" * 70)
print("Loading h5ad file...")
print("=" * 70)

adata = sc.read_h5ad(h5ad_path)
print(f"✅ Successfully loaded dataset: {adata.shape}")
print(f"   Cells: {adata.n_obs:,}")
print(f"   Genes: {adata.n_vars:,}")


# In[9]:


# =============================================================================
# INSPECT DATA STRUCTURE
# =============================================================================

print("\n" + "=" * 70)
print("Data Structure Inspection")
print("=" * 70)

# Cells (metadata)
print("\n📋 Cell metadata (.obs):")
print(f"   Columns: {len(adata.obs.columns)}")
print(f"   Columns: {list(adata.obs.columns)}")
# print("\n   First few rows:")
# print(adata.obs.head())


# In[10]:


import pandas as pd

# Show all columns
pd.set_option('display.max_columns', None)

# Also useful:
pd.set_option('display.width', None)  # Don't wrap to multiple lines
pd.set_option('display.max_colwidth', None)  # Show full column content

# Now use head()
print(adata.obs.head())


# In[11]:


# =============================================================================
# IDENTIFY CATEGORICAL COLUMNS
# =============================================================================

print("=" * 80)
print("CATEGORICAL COLUMNS - UNIQUE VALUES")
print("=" * 80)

# Get categorical columns (dtype = 'object' or 'category')
categorical_cols = []
for col in adata.obs.columns:
    if adata.obs[col].dtype == 'object' or pd.api.types.is_categorical_dtype(adata.obs[col]):
        categorical_cols.append(col)

print(f"\nFound {len(categorical_cols)} categorical columns\n")
print(categorical_cols)


# In[12]:


# =============================================================================
# SHOW UNIQUE VALUES FOR EACH CATEGORICAL COLUMN
# =============================================================================

for col in categorical_cols:
    print("=" * 80)
    print(f"Column: {col}")
    print("=" * 80)
    
    # Get unique values
    unique_vals = adata.obs[col].unique()
    n_unique = len(unique_vals)
    
    # Get value counts
    value_counts = adata.obs[col].value_counts(dropna=False)
    
    print(f"Number of unique values: {n_unique}")
    print(f"Total cells: {len(adata.obs)}")
    print(f"\nValue counts:")
    print(value_counts)
    print()


# In[13]:


# =============================================================================
# KEY COLUMNS DETAILED VIEW
# =============================================================================

print("\n" + "=" * 80)
print("KEY ANNOTATION COLUMNS - DETAILED VIEW")
print("=" * 80)

key_cols = [
    "condition",
    "class",
    "donor",
    "orig_ident",
    "labels_scanvi",
    "pre_ingest_annots",
    "post_filter_annots",
    "predictions_scanvi",
    "celltypes",
    "fine_annotation",
    "individual",
    "batch",
    "sample",
]

unique_vals = {c: adata.obs[c].unique() for c in key_cols if c in adata.obs}
for k, v in unique_vals.items():
    print(f"{k}: {v}")


for col in key_cols:
    if col in adata.obs.columns:
        print("\n" + "=" * 80)
        print(f"\n{col}:")
        print(f"  Unique values: {adata.obs[col].unique()}")
        print(f"  Value counts:")
        
        # counts = adata.obs[col].value_counts()
        counts = adata.obs[col].value_counts(dropna=False)

        for val, count in counts.items():
            pct = (count / len(adata.obs)) * 100
            print(f"    {val}: {count} ({pct:.2f}%)")


# In[14]:


# =============================================================================
# CROSSTAB: KEY RELATIONSHIPS
# =============================================================================

print("\n" + "=" * 80)
print("KEY RELATIONSHIPS")
print("=" * 80)

# Condition × Cell Type
if 'condition' in adata.obs.columns and 'post_filter_annots' in adata.obs.columns:
    print("\n1. Condition × Cell Type (post_filter_annots):")
    crosstab = pd.crosstab(adata.obs['condition'], adata.obs['post_filter_annots'], margins=True)
    print(crosstab)
print("=" * 80)

# Sample × Condition
if 'sample' in adata.obs.columns and 'condition' in adata.obs.columns:
    print("\n2. Sample × Condition:")
    crosstab = pd.crosstab(adata.obs['sample'], adata.obs['condition'], margins=True)
    print(crosstab)
print("=" * 80)

# Donor × Condition
if 'donor' in adata.obs.columns and 'condition' in adata.obs.columns:
    print("\n3. Donor × Condition:")
    crosstab = pd.crosstab(adata.obs['donor'], adata.obs['condition'], margins=True)
    print(crosstab)
print("=" * 80)

# Class × post_filter_annots
if 'class' in adata.obs.columns and 'post_filter_annots' in adata.obs.columns:
    print("\n4. Class × Cell Type (post_filter_annots):")
    crosstab = pd.crosstab(adata.obs['class'], adata.obs['post_filter_annots'], margins=True)
    print(crosstab)
print("=" * 80)

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE!")
print("=" * 80)


# In[ ]:





# In[15]:


qc_vars = [
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "total_counts_mt",
    "log1p_total_counts_mt",
    "pct_counts_mt",
    "total_counts_ribo",
    "log1p_total_counts_ribo",
    "pct_counts_ribo",
    "total_counts_hb",
    "log1p_total_counts_hb",
    "pct_counts_hb",
    "n_genes",
    "doublet_score",
    "predicted_doublet",
]

# to USE :

# orig.ident
# donor

# condition
# class

# predictions_scanvi
# post_filter_annots


# In[ ]:


# Layers and embeddings
print("\n📊 Available layers:", list(adata.layers.keys()))
print("📊 Available embeddings (.obsm):", list(adata.obsm.keys()))

# To watch : 

# orig.ident
# donor

# condition
# class

# predictions_scanvi
# post_filter_annots


# In[36]:


import scanpy as sc

embeddings = ["X_pca", "X_scANVI", "X_scVI", "X_scvi", "X_umap"]
bases = [e[2:] if e.startswith("X_") else e for e in embeddings]
print(bases)


# In[43]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

cols = [
    "condition",
    "class",
    "donor",
    "orig_ident",
    "labels_scanvi",
    "pre_ingest_annots",
    "post_filter_annots",
    "predictions_scanvi",
    "celltypes",
    "fine_annotation",
    "individual",
]

if "X_umap" not in adata.obsm:
    raise ValueError("UMAP not found in adata.obsm['X_umap'].")

for col in cols:
    if col not in adata.obs.columns:
        continue

    # Only makes sense for categorical data
    is_cat = (
        pd.api.types.is_categorical_dtype(adata.obs[col])
        or adata.obs[col].dtype == object
    )

    sc.pl.umap(
        adata,
        color=col,
        legend_loc="on data",      # labels on plot
        legend_fontsize=6,
        legend_fontweight="normal",
        frameon=False,
        title=f"umap: {col}",
        show=False,                # 👈 important
    )

    ax = plt.gca()

    # Re-create legend manually
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1.05, 0.5),
            frameon=False,
            fontsize=8,
        )

    plt.show()


# In[ ]:





# In[16]:


import pandas as pd
import matplotlib.pyplot as plt

cols = [
    "condition",
    "class",
    "donor",
    "orig_ident",
    "labels_scanvi",
    "pre_ingest_annots",
    "post_filter_annots",
    "predictions_scanvi",
    "celltypes",
    "fine_annotation",
    "individual",
]

top_n = 300  # change if you want more/less categories shown

for col in cols:
    if col not in adata.obs.columns:
        continue

    print("=" * 80)
    print(col)
    print("=" * 80)

    counts = adata.obs[col].value_counts(dropna=False)

    # Print counts (including NaN)
    for val, n in counts.items():
        label = "NaN" if pd.isna(val) else val
        print(f"{label}: {n}")
    print(f"\nTotal cells: {counts.sum()}\n")

    # --- Bar plot (top N categories) ---
    plot_counts = counts.copy()
    plot_counts.index = [
        "NaN" if pd.isna(v) else str(v) for v in plot_counts.index
    ]

    if len(plot_counts) > top_n:
        top = plot_counts.iloc[:top_n]
        other = plot_counts.iloc[top_n:].sum()
        plot_counts = pd.concat([top, pd.Series({"Other": other})])

    fig_height = max(3, 0.25 * len(plot_counts) + 1)  # scale with #categories
    fig, ax = plt.subplots(figsize=(10, fig_height))
    plot_counts.sort_values().plot(kind="barh", ax=ax)  # horizontal for readability

    ax.set_title(f"Cells per category: {col}", fontweight="normal")
    ax.set_xlabel("Number of cells", fontweight="normal")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.show()


# In[52]:


# Let's work with these LABELS :

cols = [
    "condition",
    "class",
    
    "post_filter_annots",
    "predictions_scanvi",
]

# To watch : 

# orig.ident
# donor

# condition
# class

# predictions_scanvi
# post_filter_annots


# In[17]:


print('''

Making the scanpy object compatible with GENEFORMER :

resuming the script : 

process_scanvi_h5ad_for_geneformer.py

''')


# In[19]:


import os
print(os.getcwd())


# In[20]:


# =============================================================================
# INSPECT DATA STRUCTURE
# =============================================================================

print("\n" + "=" * 70)
print("Data Structure Inspection")
print("=" * 70)

# Cells (metadata)
print("\n📋 Cell metadata (.obs):")
print(f"   Columns: {len(adata.obs.columns)}")
print(f"   Columns: {list(adata.obs.columns)}")
print("\n   First few rows:")
print(adata.obs.head())

# Genes / features (metadata)
print("\n🧬 Gene metadata (.var):")
print(f"   Columns: {len(adata.var.columns)}")
print(f"   Columns: {list(adata.var.columns)}")
print("\n   First few rows:")
print(adata.var.head())

# Layers and embeddings
print("\n📊 Available layers:", list(adata.layers.keys()))
print("📊 Available embeddings (.obsm):", list(adata.obsm.keys()))

# Save metadata to CSV
print("\n💾 Saving metadata to CSV...")
adata.obs.to_csv(f"{prefix}.cells_metadata.csv", index=True)
adata.var.to_csv(f"{prefix}.genes_metadata.csv", index=True)
print(f"   ✅ Saved: {prefix}.cells_metadata.csv")
print(f"   ✅ Saved: {prefix}.genes_metadata.csv")


# In[24]:


print("\n" + "=" * 70)
print("Data Structure Inspection")
print("=" * 70)

# Overall AnnData
print("\n📦 AnnData object:")
print(f"   Cells × Genes: {adata.shape}")

# Expression matrix
print("\n🧮 Expression matrix (.X):")
print(f"   Shape: {adata.X.shape}")
print(f"   Type: {type(adata.X)}")

# Cell metadata
print("\n📋 Cell metadata (.obs):")
print(f"   Rows (cells): {adata.obs.shape[0]}")
print(f"   Columns: {adata.obs.shape[1]}")
print(f"   Column names: {list(adata.obs.columns)}")
print("\n   First few rows:")
print(adata.obs.head())

# Gene metadata
print("\n🧬 Gene metadata (.var):")
print(f"   Rows (genes): {adata.var.shape[0]}")
print(f"   Columns: {adata.var.shape[1]}")
print(f"   Column names: {list(adata.var.columns)}")
print("\n   First few rows:")
print(adata.var.head())


# In[ ]:





# In[21]:


# =============================================================================
# GENE NAME PROCESSING
# =============================================================================

print("\n" + "=" * 70)
print("Gene Name Processing")
print("=" * 70)

# Add gene_name column
adata.var["gene_name"] = adata.var_names.astype(str)
print(f"✅ Added 'gene_name' column")

# Check gene name formats
s = adata.var["gene_name"]
adata.var["is_ENST"] = s.str.match(r"^ENST\d+(\.\d+)?$", na=False)
adata.var["is_ENSG"] = s.str.match(r"^ENSG\d+(\.\d+)?$", na=False)
adata.var["is_ENSP"] = s.str.match(r"^ENSP\d+(\.\d+)?$", na=False)

print(f"\n📊 Gene name format summary:")
print(f"   ENST transcripts: {int(adata.var['is_ENST'].sum())} / {adata.n_vars}")
print(f"   ENSG genes: {int(adata.var['is_ENSG'].sum())} / {adata.n_vars}")
print(f"   ENSP proteins: {int(adata.var['is_ENSP'].sum())} / {adata.n_vars}")

# Count genes with version suffix
mask = s.str.contains(r"\.\d+$", na=False)
n_versioned = int(mask.sum())
print(f"\n📊 Genes with version suffix (.number): {n_versioned} / {adata.n_vars}")

# Check specific genes
test_genes = ["TET1", "ORAI1"]
s_upper = s.str.upper()
print(f"\n🔍 Checking specific genes:")
for g in test_genes:
    hits = s[s_upper == g]
    print(f"   {g}: {'FOUND' if len(hits) else 'NOT FOUND'} (n={len(hits)})")


# In[22]:


# =============================================================================
# REMOVE VERSION SUFFIXES
# =============================================================================

print("\n" + "=" * 70)
print("Removing Version Suffixes")
print("=" * 70)

# Create gene_no_version column
adata.var["gene_no_version"] = s.str.replace(r"\.\d+$", "", regex=True)
print(f"✅ Created 'gene_no_version' column")

# Check for duplicates
genes_no_version = adata.var["gene_no_version"].astype(str)
n_dup_rows = int(genes_no_version.duplicated().sum())
vc = genes_no_version.value_counts()
n_dup_names = int((vc > 1).sum())

print(f"\n📊 Duplicate analysis:")
print(f"   Duplicated rows (beyond first): {n_dup_rows}")
print(f"   Unique gene names with duplicates: {n_dup_names}")

if n_dup_rows > 0:
    print(f"\n   Top duplicated names:")
    print(vc[vc > 1].head(5))


# In[23]:


# =============================================================================
# REMOVE DUPLICATES (KEEP MOST EXPRESSED)
# =============================================================================

print("\n" + "=" * 70)
print("Removing Duplicate Genes (keeping most expressed)")
print("=" * 70)

# Calculate gene expression sums
gene_sum = np.asarray(adata.X.sum(axis=0)).ravel()
df = pd.DataFrame({
    "gene": adata.var["gene_no_version"].astype(str).values,
    "s": gene_sum
})

# Keep most expressed per gene
keep_idx = df.groupby("gene")["s"].idxmax().sort_values().to_numpy()
adata = adata[:, keep_idx].copy()
adata.var_names = df["gene"].values[keep_idx]

# Also ensure uniqueness
genes = adata.var["gene_no_version"].astype(str).values
n_before = adata.n_vars
_, keep_idx = np.unique(genes, return_index=True)
keep_idx = np.sort(keep_idx)

adata2 = adata[:, keep_idx].copy()
adata2.var_names = genes[keep_idx]

print(f"✅ Kept {adata2.n_vars:,} of {n_before:,} genes")
print(f"   Removed {n_before - adata2.n_vars:,} duplicates")


# In[ ]:





# In[25]:


# =============================================================================
# MAP TO ENSEMBL IDs
# =============================================================================

print("\n" + "=" * 70)
print("Mapping Gene Names to Ensembl IDs")
print("=" * 70)

# Load gene mapping file
map_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer-pickle/gene_name_id_dict_gc95M.txt"

if not os.path.exists(map_file):
    raise FileNotFoundError(f"❌ Mapping file not found: {map_file}")

print(f"📖 Loading gene mapping file: {map_file}")
gene_map = pd.read_csv(
    map_file, sep=r"\s+", header=None, names=["gene", "ensembl_id"], dtype=str
)
print(f"   ✅ Loaded {len(gene_map):,} gene mappings")

# Create lookup table (case-insensitive)
gm = gene_map.copy()
gm["gene_upper"] = gm["gene"].str.upper()
lut = gm.set_index("gene_upper")["ensembl_id"]

# Map gene names to Ensembl IDs
s = adata2.var["gene_no_version"].astype(str)
adata2.var["ensembl_id"] = s.str.upper().map(lut)
n_match = int(adata2.var["ensembl_id"].notna().sum())

print(f"\n📊 Mapping results:")
print(f"   Matched: {n_match:,} / {adata2.n_vars:,} genes ({n_match/adata2.n_vars*100:.1f}%)")
print(f"   Unmatched: {adata2.n_vars - n_match:,} genes")

# Show sample matches
hits = pd.DataFrame({
    "gene_no_version": s,
    "ensembl_id": adata2.var["ensembl_id"]
}).dropna()
print(f"\n   Sample matches:")
print(hits.head())


# In[26]:


# =============================================================================
# FILTER TO MATCHED GENES ONLY
# =============================================================================

print("\n" + "=" * 70)
print("Filtering to Matched Genes Only")
print("=" * 70)

mask = adata2.var["ensembl_id"].notna()
adata3 = adata2[:, mask].copy()

print(f"✅ Filtered dataset: {adata3.n_obs:,} cells × {adata3.n_vars:,} genes")
print(f"   Removed {adata2.n_vars - adata3.n_vars:,} unmatched genes")


# In[27]:


# =============================================================================
# PREPARE FOR GENEFORMER
# =============================================================================

print("\n" + "=" * 70)
print("Preparing for Geneformer Compatibility")
print("=" * 70)

# Add n_counts column (required by Geneformer)
if "n_counts" not in adata3.obs.columns:
    if "total_counts" in adata3.obs.columns:
        adata3.obs["n_counts"] = adata3.obs["total_counts"]
        adata3.obs["counts"] = adata3.obs["total_counts"]
        print("✅ Added 'n_counts' and 'counts' columns from 'total_counts'")
    else:
        print("⚠️  Warning: 'total_counts' not found. Creating 'n_counts' from X.sum(axis=1)")
        adata3.obs["n_counts"] = np.asarray(adata3.X.sum(axis=1)).ravel()
        adata3.obs["counts"] = adata3.obs["n_counts"]

# Set gene names to Ensembl IDs
A = adata3
A.var["gene_symbol"] = A.var_names.astype(str)

assert "ensembl_id" in A.var.columns, "ensembl_id column missing"
A.var_names = A.var["ensembl_id"].astype(str)
A.var_names_make_unique()
A.var_names.name = None

print("✅ Set var_names to Ensembl IDs")
print(f"\n   Sample gene mapping:")
print(A.var[["gene_symbol", "ensembl_id"]].head())


# In[28]:


# =============================================================================
# GENEFORMER COMPATIBILITY CHECK
# =============================================================================

print("\n" + "=" * 70)
print("Geneformer Compatibility Check")
print("=" * 70)

checks = []

# Check 1: n_counts column
has_n_counts = "n_counts" in A.obs.columns
checks.append(("Has n_counts column", has_n_counts))
print(f"   {'✅' if has_n_counts else '❌'} Has n_counts: {has_n_counts}")

# Check 2: Ensembl ID format
ensg = A.var_names.astype(str)
mask_ensg = ensg.str.match(r"^ENSG\d+(\.\d+)?$", na=False)
n_ensg = int(mask_ensg.sum())
checks.append(("Ensembl ID format", n_ensg == len(ensg)))
print(f"   {'✅' if n_ensg == len(ensg) else '❌'} Ensembl ID format: {n_ensg}/{len(ensg)} ({n_ensg/len(ensg)*100:.1f}%)")

# Check 3: No duplicate gene IDs
dup = int(ensg.duplicated().sum())
checks.append(("No duplicate gene IDs", dup == 0))
print(f"   {'✅' if dup == 0 else '❌'} No duplicate gene IDs: {dup == 0} (duplicates={dup})")

all_checks_passed = all(check[1] for check in checks)
print(f"\n   {'✅ All checks passed!' if all_checks_passed else '⚠️  Some checks failed'}")
adata4 = A


# In[29]:


# =============================================================================
# SAVE PROCESSED DATA
# =============================================================================

print("\n" + "=" * 70)
print("Saving Processed Data")
print("=" * 70)

# Save to h5ad
out_h5ad = f"{prefix}.geneformer_compatible.h5ad"
adata4.write(out_h5ad)
print(f"✅ Saved h5ad: {out_h5ad}")

# Save to loom
loom_file = f"{prefix}.geneformer_compatible.loom"
adata4.write_loom(loom_file, write_obsm_varm=True)
print(f"✅ Saved loom: {loom_file}")


# In[30]:


# =============================================================================
# VERIFY LOOM FILE
# =============================================================================

print("\n" + "=" * 70)
print("Verifying Loom File")
print("=" * 70)

adataL = sc.read_loom(loom_file)
print(f"✅ Successfully loaded loom file")
print(f"   Cells: {adataL.n_obs:,}, Genes: {adataL.n_vars:,}")
print(f"   .obs columns: {len(adataL.obs.columns)}")
print(f"   .var columns: {len(adataL.var.columns)}")
print(f"   Layers: {list(adataL.layers.keys())}")
print(f"   Embeddings (.obsm): {list(adataL.obsm.keys())}")


# In[ ]:





# In[31]:


# =============================================================================
# CREATE VISUALIZATIONS
# =============================================================================

print("\n" + "=" * 70)
print("Creating Visualizations")
print("=" * 70)

# Set figure directory
sc.settings.figdir = "figs"
os.makedirs("figs", exist_ok=True)

# Create multi-panel embedding plot
print("📊 Creating multi-panel embedding plot...")
fig, axs = plt.subplots(1, 3, figsize=(12, 4))
for ax, basis in zip(axs, ["X_pca", "X_umap", "X_scVI"]):
    if basis in adataL.obsm.keys():
        color_col = "class" if "class" in adataL.obs.columns else None
        sc.pl.embedding(
            adataL, basis=basis,
            color=color_col,
            ax=ax, show=False, frameon=False
        )
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        ax.text(0.5, 0.5, f"{basis} not available", 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])

plt.tight_layout()
plt.savefig(f"figs/{prefix}_embeddings_panel.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"   ✅ Saved: figs/{prefix}_embeddings_panel.png")

# Create UMAP plots
if "X_umap" in adataL.obsm.keys():
    # UMAP colored by class
    if "class" in adataL.obs.columns:
        print("📊 Creating UMAP plot (colored by class)...")
        sc.pl.embedding(
            adataL, basis="umap",
            color="class",
            frameon=False, show=False,
            save=f"_{prefix}_umap_class.png"
        )
        print(f"   ✅ Saved: figs/umap_{prefix}_umap_class.png")
    
    # UMAP colored by broad_class
    if "broad_class" in adataL.obs.columns:
        print("📊 Creating UMAP plot (colored by broad_class)...")
        sc.pl.embedding(
            adataL, basis="umap",
            color="broad_class",
            frameon=False, show=False,
            save=f"_{prefix}_umap_broad_class.png"
        )
        print(f"   ✅ Saved: figs/umap_{prefix}_umap_broad_class.png")
    
    # UMAP colored by condition
    if "condition" in adataL.obs.columns:
        print("📊 Creating UMAP plot (colored by condition)...")
        sc.pl.embedding(
            adataL, basis="umap",
            color="condition",
            frameon=False, show=False,
            save=f"_{prefix}_umap_condition.png"
        )
        print(f"   ✅ Saved: figs/umap_{prefix}_umap_condition.png")

# Create UMAP panels with multiple colors
if "X_umap" in adataL.obsm.keys():
    color_cols = []
    for col in ["leiden", "broad_class", "class", "condition"]:
        if col in adataL.obs.columns:
            color_cols.append(col)
    
    if color_cols:
        print(f"📊 Creating UMAP panels (colored by: {', '.join(color_cols)})...")
        sc.pl.umap(
            adataL, color=color_cols,
            save=f"_{prefix}_umap_panels.png"
        )
        print(f"   ✅ Saved: figs/umap_{prefix}_umap_panels.png")


# In[ ]:





# In[32]:


# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\n📁 Input file: {h5ad_file}")
print(f"   Original: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"\n📁 Output files:")
print(f"   ✅ {out_h5ad}")
print(f"   ✅ {loom_file}")
print(f"   ✅ {prefix}.cells_metadata.csv")
print(f"   ✅ {prefix}.genes_metadata.csv")
print(f"   ✅ Visualizations in figs/ directory")

print(f"\n📊 Processing summary:")
print(f"   Genes after deduplication: {adata2.n_vars:,}")
print(f"   Genes with Ensembl ID match: {adata3.n_vars:,}")
print(f"   Final dataset: {adata4.n_obs:,} cells × {adata4.n_vars:,} genes")

print(f"\n✅ Processing complete!")
print("=" * 70)


# In[ ]:




