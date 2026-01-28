#!/usr/bin/env python
# coding: utf-8

# In[1]:


print(
    
"""
Script to process h5ad file for Geneformer compatibility
Updated for: SCANVI_INTEGRATED_FILTERED.h5ad
Location: /mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026
"""

)


# In[2]:


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


# In[3]:


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


# In[4]:


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


# In[5]:


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


# In[6]:


import pandas as pd

# Show all columns
pd.set_option('display.max_columns', None)

# Also useful:
pd.set_option('display.width', None)  # Don't wrap to multiple lines
pd.set_option('display.max_colwidth', None)  # Show full column content

# Now use head()
print(adata.obs.head())


# In[9]:


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


# In[10]:


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


# In[11]:


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


# In[24]:


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





# In[17]:


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


# In[18]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Define QC variables
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

# Define categorical columns
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

# Create output directory
import os
os.makedirs('qc_violin_plots', exist_ok=True)

# Filter to existing columns
qc_vars = [q for q in qc_vars if q in adata.obs.columns]
cols = [c for c in cols if c in adata.obs.columns]

print(f"Plotting {len(qc_vars)} QC variables against {len(cols)} categorical columns...")
print(f"Total plots to generate: {len(qc_vars) * len(cols)}")

# Create a copy of data
adata_plot = adata.copy()

# Loop through each categorical column
for col in cols:
    print(f"\n{'='*80}")
    print(f"Plotting QC metrics grouped by: {col}")
    print(f"{'='*80}")
    
    # Clean the column - remove old color palette
    if f'{col}_colors' in adata_plot.uns:
        del adata_plot.uns[f'{col}_colors']
    
    # Convert to string to handle NaN
    adata_plot.obs[f'{col}_clean'] = adata_plot.obs[col].astype(str).replace('nan', 'NaN')
    
    # Get number of unique categories
    n_categories = adata_plot.obs[f'{col}_clean'].nunique()
    print(f"  Number of categories: {n_categories}")
    
    # Determine rotation based on number of categories
    rotation = 90 if n_categories > 5 else 45
    
    try:
        # Plot violin plots for this categorical column
        sc.pl.violin(
            adata_plot,
            qc_vars,
            groupby=f'{col}_clean',
            rotation=rotation,
            save=f'_{col}_all_qc.png',
            show=False
        )
        print(f"  ✅ Saved: figures/violin_{col}_all_qc.png")
        
    except Exception as e:
        print(f"  ⚠️  Error plotting {col}: {e}")
        continue

print("\n" + "="*80)
print("COMPLETE!")
print("="*80)
print(f"Plots saved in: figures/")


# In[ ]:


# to USE :

# orig.ident
# donor

# condition
# class

# predictions_scanvi
# post_filter_annots


# In[ ]:





# In[20]:


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from IPython.display import display

# Define variables
qc_vars = [
    "n_genes_by_counts",
    "total_counts",
    "total_counts_mt",
    "pct_counts_mt",
    "total_counts_ribo",
    "pct_counts_ribo",
    "total_counts_hb",
    "pct_counts_hb",
    "doublet_score",
]

# Your specified columns
cols = [
   # "orig.ident",
   # "donor",
    "condition",
    "class",
    "predictions_scanvi",
    "post_filter_annots",
]

# Filter existing
qc_vars = [q for q in qc_vars if q in adata.obs.columns]
cols = [c for c in cols if c in adata.obs.columns]

# Create output directory with a different name
output_dir = 'figures_qc_violin_by_category'
os.makedirs(output_dir, exist_ok=True)

# Prepare data
df = adata.obs.copy()

# Loop through each categorical column
for col in cols:
    print(f"\n{'='*80}")
    print(f"Plotting: {col}")
    print(f"{'='*80}")
    
    # Convert to string to handle NaN
    df[col] = df[col].astype(str).replace('nan', 'NaN')
    
    # Create grid of violin plots
    n_qc = len(qc_vars)
    n_plot_cols = 3
    n_rows = (n_qc + n_plot_cols - 1) // n_plot_cols
    
    fig, axes = plt.subplots(n_rows, n_plot_cols, figsize=(18, 5*n_rows))
    axes = axes.flatten()
    
    for idx, qc_var in enumerate(qc_vars):
        ax = axes[idx]
        
        # Create violin plot
        sns.violinplot(
            data=df,
            x=col,
            y=qc_var,
            ax=ax,
            inner='box'
        )
        
        ax.set_title(qc_var, fontsize=11, fontweight='bold')
        ax.set_xlabel(col, fontsize=10)
        ax.set_ylabel('', fontsize=9)
        ax.tick_params(axis='x', rotation=45, labelsize=8)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # Remove extra subplots
    for idx in range(len(qc_vars), len(axes)):
        fig.delaxes(axes[idx])
    
    plt.suptitle(f'QC Metrics by {col}', fontsize=16, fontweight='bold', y=1.01)
    plt.tight_layout()
    
    # Save figure
    save_path = f'{output_dir}/qc_by_{col}_grid.png'
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"✅ Saved: {save_path}")
    
    # Display in Jupyter
    plt.show()
    
    # Add spacing between plots
    print("\n")

print("="*80)
print("ALL PLOTS COMPLETE!")
print("="*80)
print(f"All figures saved in: {output_dir}/")


# In[27]:


# VIOLIN PLOTS BY SAMPLE
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             groupby='condition', rotation=45)


# In[33]:


# VIOLIN PLOTS BY SAMPLE
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             groupby='class', rotation=45)


# In[24]:


# VIOLIN PLOTS BY SAMPLE
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             groupby='predictions_scanvi', rotation=45)


# In[32]:


# VIOLIN PLOTS BY SAMPLE
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             groupby='post_filter_annots', rotation=90)


# In[ ]:





# In[34]:


# Genes / features (metadata)
print("\n🧬 Gene metadata (.var):")
print(f"   Columns: {len(adata.var.columns)}")
print(f"   Columns: {list(adata.var.columns)}")
print("\n   First few rows:")
print(adata.var.head())


# In[35]:


# Layers and embeddings
print("\n📊 Available layers:", list(adata.layers.keys()))
print("📊 Available embeddings (.obsm):", list(adata.obsm.keys()))


# In[38]:


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


# In[39]:


import scanpy as sc

# embeddings = ["X_pca", "X_scANVI", "X_scVI", "X_scvi", "X_umap"]
# bases = ["pca", "umap", "scANVI", "scVI", "scvi"]  # basis names, not X_ keys

cols = [
    "condition",
    "class",

    "orig_ident",
    "donor",

    
    "labels_scanvi",
    "predictions_scanvi",
    
    "pre_ingest_annots",
    "post_filter_annots",
    
    # "celltypes",
    # "fine_annotation",
    # "individual",
]

# Keep only bases that exist as X_<basis> in obsm
bases = [b for b in bases if f"X_{b}" in adata.obsm]

for b in bases:
    for col in cols:
        if col not in adata.obs.columns:
            continue
        sc.pl.embedding(
            adata,
            basis=b,
            color=col,
            frameon=False,
            legend_loc="on data" if adata.obs[col].nunique(dropna=True) < 15 else "right margin",
            title=f"{b}: {col}",
        )


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


# In[41]:


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
    
    if not is_cat:
        continue
    
    # Create figure with extra space for legend
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot with labels on data
    sc.pl.umap(
        adata,
        color=col,
        legend_loc="on data",
        legend_fontsize=8,
        legend_fontweight="normal",
        legend_fontoutline=2,
        frameon=False,
        title=f"UMAP: {col}",
        ax=ax,
        show=False
    )
    
    # Now add a separate legend on the right
    # Get unique categories and their colors
    if f"{col}_colors" in adata.uns:
        colors = adata.uns[f"{col}_colors"]
        categories = adata.obs[col].cat.categories if pd.api.types.is_categorical_dtype(adata.obs[col]) else adata.obs[col].unique()
        
        # Create legend handles
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=colors[i], label=cat) 
            for i, cat in enumerate(categories) if i < len(colors)
        ]
        
        # Add legend on the right
        ax.legend(
            handles=legend_elements,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            frameon=True,
            fontsize=8,
            title=col,
            title_fontsize=9
        )
    
    plt.tight_layout()
    plt.show()


# In[ ]:





# In[45]:


import matplotlib.pyplot as plt

# cols = ["post_filter_annots", "condition", "sample", "leiden"]

cols = [
 #   "condition",
    "class",
 #   "donor",
 #   "orig_ident",
    "labels_scanvi",
 #   "pre_ingest_annots",
 #   "post_filter_annots",
    "predictions_scanvi",
 #   "celltypes",
 #   "fine_annotation",
 #   "individual",
]


fig, axes = plt.subplots(
    nrows=len(cols),
    ncols=1,
    figsize=(10, 6 * len(cols)),   # 👈 key change
)

for ax, col in zip(axes, cols):
    sc.pl.umap(
        adata,
        color=col,
        ax=ax,
        show=False,
        frameon=False,
        legend_loc="on data" if adata.obs[col].nunique() < 20 else "right margin",
    )
    ax.set_title(f"UMAP: {col}", fontweight="normal")
    ax.set_aspect("equal", adjustable="box")   # keep geometry correct
    ax.axis("off")

plt.tight_layout()
plt.show()



# In[ ]:





# In[46]:


import scanpy as sc
import matplotlib.pyplot as plt

# Load data
# adata = sc.read_h5ad("SCANVI_INTEGRATED_FILTERED.h5ad")

# =============================================================================
# SET GLOBAL FONT SIZES
# =============================================================================

# Set scanpy defaults
sc.settings.set_figure_params(
    dpi=150,
    fontsize=8,           # General font size
    figsize=(6, 6),
    frameon=False
)

# Set matplotlib defaults
plt.rcParams['font.size'] = 8
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['axes.labelsize'] = 8
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7
plt.rcParams['legend.fontsize'] = 7
plt.rcParams['legend.title_fontsize'] = 8

# =============================================================================
# NOW YOUR PLOTS WITH SMALLER FONTS
# =============================================================================

# 1. UMAP with main cell type annotations
sc.pl.umap(
    adata, 
    color='post_filter_annots', 
    legend_loc='on data',
    legend_fontsize=6,          # Legend font size
    legend_fontoutline=2,       # Outline around text
    title='Cell Types'          # Title text (fontsize controlled by rcParams)
)

# 2. UMAP with multiple annotations
sc.pl.umap(
    adata, 
    color=['post_filter_annots', 'condition', 'sample', 'leiden'], 
    ncols=2,
    legend_fontsize=6
)

# 3. All embeddings side-by-side with cell types
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, emb, title in zip(axes, ['X_scVI', 'X_scANVI', 'X_umap'], ['scVI', 'scANVI', 'UMAP']):
    sc.pl.embedding(
        adata, 
        basis=emb, 
        color='post_filter_annots', 
        ax=ax, 
        show=False, 
        title=title,
        legend_fontsize=6,        # Legend font size
        legend_fontoutline=1      # Outline around legend text
    )
    # Adjust title font size manually
    ax.set_title(title, fontsize=10)

plt.tight_layout()
plt.show()


# In[47]:


import scanpy as sc
import matplotlib.pyplot as plt

# Load data
# adata = sc.read_h5ad("SCANVI_INTEGRATED_FILTERED.h5ad")

# Set small fonts globally
sc.settings.set_figure_params(dpi=150, fontsize=8, figsize=(8, 8), frameon=False)

# Create UMAP with post_filter_annots
sc.pl.umap(
    adata, 
    color='post_filter_annots',
    legend_loc='right margin',
    legend_fontsize=6,
    legend_fontoutline=1,
    title='UMAP',
    save='_post_filter_annots.png'  # Optional: saves to figures/umap_post_filter_annots.png
)


# In[ ]:





# In[49]:


import pandas as pd

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

for col in cols:
    if col not in adata.obs.columns:
        continue

    print("=" * 80)
    print(col)
    print("=" * 80)

    counts = adata.obs[col].value_counts(dropna=False)

    for val, n in counts.items():
        label = "NaN" if pd.isna(val) else val
        print(f"{label}: {n}")

    print(f"\nTotal cells: {counts.sum()}")
    print()


# In[50]:


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


# To watch : 

# orig.ident
# donor

# condition
# class

# predictions_scanvi
# post_filter_annots


# In[53]:


# Let's work with these LABELS :

cols = [
    "condition",
    "class",
    
    "post_filter_annots",
    "predictions_scanvi",
]


# In[55]:


import itertools


for col1, col2 in itertools.combinations(cols, 2):
    if col1 not in adata.obs.columns or col2 not in adata.obs.columns:
        continue

    print("=" * 90)
    print(f"{col1}  ×  {col2}")
    print("=" * 90)

    ct = pd.crosstab(
        adata.obs[col1],
        adata.obs[col2],
        dropna=False
    )

    print(ct)
    print(f"\nTotal cells: {ct.values.sum()}\n")



# In[56]:


import pandas as pd
import itertools
import matplotlib.pyplot as plt


# Only these combinations will be shown
pairs_to_plot = [
    ("condition", "class"),
    ("condition", "post_filter_annots"),
    ("condition", "predictions_scanvi"),
]

for col1, col2 in pairs_to_plot:
    # Safety checks
    if col1 not in adata.obs.columns or col2 not in adata.obs.columns:
        continue

    print("=" * 90)
    print(f"{col1} × {col2}")
    print("=" * 90)

    ct = pd.crosstab(
        adata.obs[col1],
        adata.obs[col2],
        dropna=False
    )

    # Make NaN explicit for display
    ct.index = ["NaN" if pd.isna(v) else str(v) for v in ct.index]
    ct.columns = ["NaN" if pd.isna(v) else str(v) for v in ct.columns]

    # Sort rows by total count (recommended for readability)
    ct = ct.loc[ct.sum(axis=1).sort_values().index]

    # Horizontal grouped bar plot
    ax = ct.plot(
        kind="barh",
        stacked=False,
        figsize=(12, max(4, 0.5 * len(ct))),
        width=0.85
    )

    ax.set_title(f"{col1} × {col2} (grouped counts)", fontweight="normal")
    ax.set_xlabel("Number of cells", fontweight="normal")
    ax.set_ylabel(col1, fontweight="normal")

    ax.legend(
        title=col2,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=False,
    )

    plt.tight_layout()
    plt.show()


# In[57]:


from IPython.display import display
import pandas as pd
import os

for col1, col2 in pairs_to_plot:
    if col1 not in adata.obs.columns or col2 not in adata.obs.columns:
        continue

    ct = pd.crosstab(
        adata.obs[col1],
        adata.obs[col2],
        dropna=False
    )

    ct.index = ["NaN" if pd.isna(v) else str(v) for v in ct.index]
    ct.columns = ["NaN" if pd.isna(v) else str(v) for v in ct.columns]

    ct = ct.loc[ct.sum(axis=1).sort_values().index]

    print("=" * 90)
    print(f"{col1} × {col2} (cell counts)")
    print("=" * 90)

    display(ct)

out_dir = "pairwise_cell_counts_csv"
os.makedirs(out_dir, exist_ok=True)

for col1, col2 in pairs_to_plot:
    if col1 not in adata.obs.columns or col2 not in adata.obs.columns:
        continue

    ct = pd.crosstab(
        adata.obs[col1],
        adata.obs[col2],
        dropna=False
    )

    # Make NaN explicit
    ct.index = ["NaN" if pd.isna(v) else str(v) for v in ct.index]
    ct.columns = ["NaN" if pd.isna(v) else str(v) for v in ct.columns]

    # Sort rows by total count
    ct = ct.loc[ct.sum(axis=1).sort_values().index]

    fname = f"{col1}_X_{col2}_cell_counts.csv"
    fpath = os.path.join(out_dir, fname)

    ct.to_csv(fpath)

    print(f"Saved: {fpath}")


# In[ ]:





# In[58]:


import pandas as pd
import matplotlib.pyplot as plt

# Only these combinations will be shown
pairs_to_plot = [
    ("condition", "class"),
    ("condition", "post_filter_annots"),
    ("condition", "predictions_scanvi"),
]

for col1, col2 in pairs_to_plot:
    # Safety checks
    if col1 not in adata.obs.columns or col2 not in adata.obs.columns:
        continue

    print("=" * 90)
    print(f"{col1} × {col2}")
    print("=" * 90)

    # Raw counts
    ct = pd.crosstab(
        adata.obs[col1],
        adata.obs[col2],
        dropna=False
    )

    # Make NaN explicit
    ct.index = ["NaN" if pd.isna(v) else str(v) for v in ct.index]
    ct.columns = ["NaN" if pd.isna(v) else str(v) for v in ct.columns]

    # Sort rows by total count
    ct = ct.loc[ct.sum(axis=1).sort_values().index]

    # 🔹 Convert to percentages (row-wise)
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    # Plot percentages
    ax = ct_pct.plot(
        kind="barh",
        stacked=False,
        figsize=(12, max(4, 0.5 * len(ct_pct))),
        width=0.85
    )

    ax.set_title(f"{col1} × {col2} (% of cells per {col1})", fontweight="normal")
    ax.set_xlabel("Percentage of cells (%)", fontweight="normal")
    ax.set_ylabel(col1, fontweight="normal")

    ax.legend(
        title=col2,
        bbox_to_anchor=(1.02, 0.5),
        loc="center left",
        frameon=False,
    )

    plt.tight_layout()
    plt.show()


# In[69]:


from IPython.display import display
import pandas as pd
import os

for col1, col2 in pairs_to_plot:
    if col1 not in adata.obs.columns or col2 not in adata.obs.columns:
        continue

    ct = pd.crosstab(
        adata.obs[col1],
        adata.obs[col2],
        dropna=False
    )

    ct.index = ["NaN" if pd.isna(v) else str(v) for v in ct.index]
    ct.columns = ["NaN" if pd.isna(v) else str(v) for v in ct.columns]

    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    print("=" * 90)
    print(f"{col1} × {col2} (% per {col1})")
    print("=" * 90)

    display(ct_pct.round(2))

out_dir = "pairwise_cell_counts_csv"
os.makedirs(out_dir, exist_ok=True)

for col1, col2 in pairs_to_plot:
    if col1 not in adata.obs.columns or col2 not in adata.obs.columns:
        continue

    # --- RAW COUNTS ---
    ct = pd.crosstab(
        adata.obs[col1],
        adata.obs[col2],
        dropna=False
    )

    # Make NaN explicit
    ct.index = ["NaN" if pd.isna(v) else str(v) for v in ct.index]
    ct.columns = ["NaN" if pd.isna(v) else str(v) for v in ct.columns]

    # Sort rows by total count
    ct = ct.loc[ct.sum(axis=1).sort_values().index]

    # --- PERCENTAGES (row-wise) ---
    ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

    print("=" * 90)
    print(f"{col1} × {col2} (% per {col1})")
    print("=" * 90)

    # Display percentages in Jupyter
    display(ct_pct.round(2))

    # --- SAVE COUNTS ---
    fname = f"{col1}_X_{col2}_cell_counts.csv"
    fpath = os.path.join(out_dir, fname)
    ct.to_csv(fpath)

    print(f"Saved counts to: {fpath}\n")


# In[ ]:





# In[59]:


import matplotlib.pyplot as plt

# Cell type proportions (counts)
adata.obs['post_filter_annots'].value_counts().plot(kind='barh')

plt.xlabel('Number of cells')
plt.ylabel('Post-filter annotations')
plt.tight_layout()
plt.show()


# In[60]:


import matplotlib.pyplot as plt

# Cell type proportions (counts)
adata.obs['predictions_scanvi'].value_counts().plot(kind='barh')

plt.xlabel('Number of cells')
plt.ylabel('Post-filter annotations')
plt.tight_layout()
plt.show()


# In[76]:


import matplotlib.pyplot as plt

# Cell type proportions (counts)
adata.obs['condition'].value_counts().plot(kind='barh')

plt.xlabel('Number of cells')
plt.ylabel('Condition')
plt.tight_layout()
plt.show()


# In[ ]:





# In[79]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ct_pf_pred = pd.crosstab(
    adata.obs["post_filter_annots"],
    adata.obs["predictions_scanvi"],
    dropna=False
)

# Make NaN explicit
ct_pf_pred.index = ["NaN" if pd.isna(v) else str(v) for v in ct_pf_pred.index]
ct_pf_pred.columns = ["NaN" if pd.isna(v) else str(v) for v in ct_pf_pred.columns]

# Optional: reorder columns by total abundance
ct_pf_pred = ct_pf_pred.loc[:, ct_pf_pred.sum().sort_values(ascending=False).index]

plt.figure(figsize=(1.2 * len(ct_pf_pred.columns), 0.5 * len(ct_pf_pred)))
sns.heatmap(
    ct_pf_pred,
    annot=True,
    fmt="d",
    cmap="viridis",
    linewidths=0.5,
)

plt.title("post_filter_annots × predictions_scanvi (cell counts)", fontweight="normal")
plt.xlabel("predictions_scanvi")
plt.ylabel("post_filter_annots")
plt.tight_layout()
plt.show()


# In[84]:


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

print("Row-wise z-score \n")
print("For each post_filter_annots, highlight which predictions_scanvi are enriched or depleted relative to that row’s baseline.\n")
print("“For each annotation, which predictions are enriched? Row-wise \n")

# Z-score per row
ct_z = ct_pf_pred.sub(ct_pf_pred.mean(axis=1), axis=0) \
                 .div(ct_pf_pred.std(axis=1).replace(0, np.nan), axis=0)

plt.figure(figsize=(1.2 * len(ct_z.columns), 0.5 * len(ct_z)))
sns.heatmap(
    ct_z,
    cmap="vlag",
    center=0,
    linewidths=0.5,
    annot=False,
)

plt.title("post_filter_annots × predictions_scanvi (row-wise z-score)", fontweight="normal")
plt.xlabel("predictions_scanvi")
plt.ylabel("post_filter_annots")
plt.tight_layout()
plt.show()


# In[85]:


print("Column-wise z-score")
print("For each predicted class, see which post_filter_annots contribute more or less than average.")
print("For each prediction, which annotations dominate? Column-wise")
      
# Z-score per column

ct_z_col = ct_pf_pred.sub(ct_pf_pred.mean(axis=0), axis=1) \
                     .div(ct_pf_pred.std(axis=0).replace(0, np.nan), axis=1)

plt.figure(figsize=(1.2 * len(ct_z_col.columns), 0.5 * len(ct_z_col)))
sns.heatmap(
    ct_z_col,
    cmap="vlag",
    center=0,
    linewidths=0.5,
    annot=False,
)

plt.title("post_filter_annots × predictions_scanvi (column-wise z-score)", fontweight="normal")
plt.xlabel("predictions_scanvi")
plt.ylabel("post_filter_annots")
plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


sns.heatmap(
    ct_pf_pred,
    zscore=1,        # 1 = rows, 0 = columns
    cmap="vlag",
    center=0,
    linewidths=0.5,
)


# In[80]:


ct_class_pred = pd.crosstab(
    adata.obs["class"],
    adata.obs["predictions_scanvi"],
    dropna=False
)

# Make NaN explicit
ct_class_pred.index = ["NaN" if pd.isna(v) else str(v) for v in ct_class_pred.index]
ct_class_pred.columns = ["NaN" if pd.isna(v) else str(v) for v in ct_class_pred.columns]

# Optional: reorder columns by total abundance
ct_class_pred = ct_class_pred.loc[:, ct_class_pred.sum().sort_values(ascending=False).index]

plt.figure(figsize=(1.2 * len(ct_class_pred.columns), 0.5 * len(ct_class_pred)))
sns.heatmap(
    ct_class_pred,
    annot=True,
    fmt="d",
    cmap="viridis",
    linewidths=0.5,
)

plt.title("class × predictions_scanvi (cell counts)", fontweight="normal")
plt.xlabel("predictions_scanvi")
plt.ylabel("class")
plt.tight_layout()
plt.show()


# In[ ]:





# In[ ]:





# In[ ]:




