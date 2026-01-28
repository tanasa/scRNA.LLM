#!/usr/bin/env python
# coding: utf-8

# In[1]:


print(
    
"""
Script to process h5ad file finto a Hugging Face dataset :
Updated for: SCANVI_INTEGRATED_FILTERED.geneformer_compatible.h5ad

Location: /mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026
"""

)

# reuse the old script : process_scanvi_h5ad_for_geneformer.py


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
h5ad_file = "SCANVI_INTEGRATED_FILTERED.geneformer_compatible.h5ad"
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


# In[7]:


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


# In[8]:


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


# In[ ]:





# In[9]:


# qc_vars = [
#    "n_genes_by_counts",
#    "total_counts",
#    "pct_counts_mt",
#    "total_counts_mt",
#    "log1p_total_counts_mt",
#    "pct_counts_mt",
#    "total_counts_ribo",
#    "log1p_total_counts_ribo",
#    "pct_counts_ribo",
#    "total_counts_hb",
#    "log1p_total_counts_hb",
#    "pct_counts_hb",
#    "n_genes",
#    "doublet_score",
#    "predicted_doublet",
#v]

# to USE :

# orig.ident
# donor

# condition
# class

# predictions_scanvi
# post_filter_annots


# In[10]:


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


# In[11]:


import scanpy as sc

embeddings = ["X_pca", "X_scANVI", "X_scVI", "X_scvi", "X_umap"]
bases = [e[2:] if e.startswith("X_") else e for e in embeddings]
print(bases)


# In[12]:


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





# In[13]:


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


# In[14]:


# Let's work with these LABELS :

# orig.ident
# donor
# condition
# class
# predictions_scanvi
# post_filter_annots

# Found 18 categorical columns

print("""

To include all these labels into the dataset data structure ?

['sample', 
'condition', 
'leiden', 
'orig.ident', 
'class', 
'pre_ingest_annots', 
'orig_ident', 
'labels_scanvi', 
'predictions_scanvi', 
'celltypes', 
######### 'dataset', 
'fine_annotation', 
'barcode', 
'individual', 
'batch', 
'donor', 
'post_filter_annots', 
'post_filter_annots_number'
]

""")


# In[15]:


print('''

Making the scanpy object compatible with GENEFORMER :

resuming the script : 

process_scanvi_h5ad_for_geneformer.py

''')


# In[16]:


import os
print(os.getcwd())


# In[17]:


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


# In[18]:


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


# In[19]:


cell_stats = (
    adata.obs
    .groupby(["condition", "post_filter_annots"])
    .agg(
        n_cells=("condition", "size"),
        mean_genes=("n_genes_by_counts", "mean"),
        median_genes=("n_genes_by_counts", "median"),
        std_genes=("n_genes_by_counts", "std"),
        mean_umis=("total_counts", "mean"),
        median_umis=("total_counts", "median"),
    )
    .reset_index()
    .sort_values(["condition", "post_filter_annots"])
)

cell_stats


# In[20]:


print(" The number of cells per condition :")
adata.obs["condition"].value_counts()


# In[21]:


import os
import sys

print("✅ Environment set to torch_cu")
print("Python executable:", sys.executable)
print("CONDA_PREFIX:", os.environ.get("CONDA_PREFIX"))
print("CONDA_DEFAULT_ENV:", os.environ.get("CONDA_DEFAULT_ENV"))

# Verify the setup
try:
    import torch
    print("✅ PyTorch available:", torch.__version__)
    print("✅ CUDA available:", torch.cuda.is_available())
except ImportError as e:
    print("❌ PyTorch import failed:", e)

try:
    from geneformer import TranscriptomeTokenizer
    print("✅ Geneformer available")
except ImportError as e:
    print("❌ Geneformer import failed:", e)



# In[ ]:





# In[22]:


# Set working directory
os.chdir("/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.input")
print("Files:", os.listdir('.'))


# In[23]:


# =====================================
# Paths and Directories
# =====================================
import os
from geneformer import TranscriptomeTokenizer

# Base directories
workdir   = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/"
gf_dir    = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer"  # Geneformer dictionaries


# In[24]:


# Input/output directories
h5ad_input_dir = os.path.join(workdir, "SCANVI_INTEGRATED_FILTERED.geneformer_compatible.input")
out_dir        = os.path.join(workdir, "SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset")
out_prefix     = "hypoxia.20jan"


# In[25]:


# Dictionary files
gene_median_file     = os.path.join(gf_dir, "gene_median_dictionary_gc104M.pkl")
token_dictionary_file = os.path.join(gf_dir, "token_dictionary_gc104M.pkl")
gene_mapping_file     = os.path.join(gf_dir, "ensembl_mapping_dict_gc104M.pkl")


# In[26]:


# Sanity check
print("Working dir:", workdir)
print("Input h5ad dir:", h5ad_input_dir)
print("Output dir:", out_dir)
print("Dictionary files:")
for f in [gene_median_file, token_dictionary_file, gene_mapping_file]:
    print(" -", f, "✓" if os.path.exists(f) else "❌ MISSING")


# In[28]:


# =====================================
# Load and inspect h5ad
# ===================================== : to add "cell_id" 

h5ad_file="SCANVI_INTEGRATED_FILTERED.geneformer_compatible.h5ad"
prefix = Path(h5ad_file).stem
adata = sc.read_h5ad(h5ad_file)
print(f"✅ Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

print("\n.var columns:", list(adata.var.columns))
print(adata.var.head())

print("\n.obs columns:", list(adata.obs.columns))
print(adata.obs.head())
print(f"\nCells: {adata.n_obs}")

print("\nCurrent directory:", os.getcwd())
print("Files in directory:")
for f in os.listdir():
    print(" -", f)

print(prefix)


# In[29]:


print("\n.obs columns:", list(adata.obs.columns))
print(adata.obs.head(5))
print(f"\nCells: {adata.n_obs}")


# In[30]:


# the file was modified by adding cell_id


# In[31]:


# Add cell_id column to the observations
adata.obs['cell_id'] = adata.obs.index
print(adata.obs.head())


# In[32]:


# the file was modified by adding cell_id


# In[33]:


cols = ['sample', 
'condition', 
'leiden', 
'orig.ident', 
'class', 
'pre_ingest_annots', 
'orig_ident', 
'labels_scanvi', 
'predictions_scanvi', 
'celltypes', 
######### 'dataset', 
'fine_annotation', 
######### 'barcode', 
'individual', 
'batch', 
'donor', 
'post_filter_annots', 
'post_filter_annots_number',
######### 'cell_id',
]

print(cols)

for col in cols:
    if col in adata.obs.columns:
        print(f"\n=== {col} ===")
        uniques = sorted(map(str, adata.obs[col].dropna().unique()))
        print(f"Unique ({len(uniques)}): {uniques}")
        print("\nCounts:")
        print(adata.obs[col].value_counts(dropna=False))
    else:
        print(f"\n=== {col} ===")
        print("Column not found in adata.obs")


# In[34]:


cols = ['sample', 
'condition', 
'leiden', 
'orig.ident', 
'class', 
'pre_ingest_annots', 
'orig_ident', 
'labels_scanvi', 
'predictions_scanvi', 
'celltypes', 
######### 'dataset', 
'fine_annotation', 
'barcode', 
'individual', 
'batch', 
'donor', 
'post_filter_annots', 
'post_filter_annots_number',
'cell_id',
]


# In[35]:


cols = [c for c in cols if c in adata.obs.columns]

for c in cols:
    # make missing explicit and force string
    adata.obs[c] = adata.obs[c].astype("string")      # keeps <NA> properly
    adata.obs[c] = adata.obs[c].fillna("NA")          # replace missing
    adata.obs[c] = adata.obs[c].astype(str)           # ensure python str, not pandas scalar

# write fixed file


# In[36]:


# Save with a modified filename
output_file = h5ad_file.replace('.h5ad', '_with_cell_id.h5ad')
adata.write_h5ad(output_file)
print(f"✅ Saved modified h5ad with cell_id: {output_file}")


# In[37]:


# =====================================
# Tokenizer (Geneformer V2 ready)
# =====================================
tokenizer = TranscriptomeTokenizer(
    custom_attr_name_dict={
        # "cell_id": "obs_names",  # 👈 From AnnData.obs.index (cell barcodes)
        "cell_id": "cell_id",  # 👈 From AnnData.obs.index (cell barcodes)
        "sample": "sample",
        "condition": "condition",
        "leiden": "leiden",
        "orig.ident": "orig.ident",  # Note: maps to column with dot
        "class": "class",
        "pre_ingest_annots": "pre_ingest_annots",
        "labels_scanvi": "labels_scanvi",
        "predictions_scanvi": "predictions_scanvi",
        "celltypes": "celltypes",
        "fine_annotation": "fine_annotation",
        "individual": "individual",
        "batch": "batch",
        "donor": "donor",
        "post_filter_annots": "post_filter_annots",
        "post_filter_annots_number": "post_filter_annots_number",
    },
    gene_median_file = gene_median_file,
    token_dictionary_file = token_dictionary_file,
    gene_mapping_file = gene_mapping_file,
    collapse_gene_ids = True,
    model_input_size = 4096,
    special_token = True,
    model_version = "V2",
)


# In[38]:


# Run tokenization
tokenizer.tokenize_data(
    data_directory = h5ad_input_dir,
    output_directory = out_dir,
    output_prefix = out_prefix,
    file_format = "h5ad"
)

print("✅ Tokenization complete:", out_prefix)


# In[40]:


import os, glob

print("Input dir exists:", os.path.isdir(h5ad_input_dir))
print("Input dir contents:", os.listdir(h5ad_input_dir) if os.path.isdir(h5ad_input_dir) else "N/A")

print("Out dir exists:", os.path.isdir(out_dir))
print("Out dir contents:", os.listdir(out_dir) if os.path.isdir(out_dir) else "N/A")

# Look for any *.dataset subfolders right under out_dir
candidates = [p for p in glob.glob(os.path.join(out_dir, "*.dataset")) if os.path.isdir(p)]
print("Direct .dataset candidates:", candidates)


# In[41]:


# =====================================
# Load tokenized dataset (robust: auto-detect .dataset folder)
# =====================================
from datasets import load_from_disk
import os
import glob  # Add this import

# Preferred path if naming matches:
preferred_dataset_dir = os.path.join(out_dir, f"{out_prefix}.dataset")

def find_dataset_dir(base_dir, preferred=None):
    # 1) Use preferred if it exists
    if preferred and os.path.isdir(preferred):
        return preferred
    # 2) Otherwise, search for any "*.dataset" folder with HF signature files
    candidates = [p for p in glob.glob(os.path.join(base_dir, "*.dataset")) if os.path.isdir(p)]
    for cand in candidates:
        files = set(os.listdir(cand))
        if {"dataset_info.json", "state.json"} <= files:
            return cand
    return None

dataset_dir = find_dataset_dir(out_dir, preferred_dataset_dir)
print("Looking for dataset at:", preferred_dataset_dir)
if dataset_dir is None:
    # Give helpful diagnostics and stop early
    print("Contents of out_dir:", os.listdir(out_dir))
    raise FileNotFoundError(
        f"No HuggingFace dataset folder found under: {out_dir}\n"
        f"Expected something like: {out_prefix}.dataset with dataset_info.json and state.json"
    )

print("Found dataset at:", dataset_dir)
token_data = load_from_disk(dataset_dir)


# In[42]:


print("\n• Dataset structure:")
print(f"Type: {type(token_data)}")
print(f"Features: {token_data.features}")
print(f"• Number of tokenized samples: {len(token_data)}")

# Show a few samples (guard in case dataset is small)
def safe_show(ds, i, label=None):
    if i < len(ds):
        print(f"\n• {label or f'Sample #{i}'}:", ds[i])

safe_show(token_data, 0, "First sample")
safe_show(token_data, 1, "Second sample")
safe_show(token_data, 2, "Third sample")
safe_show(token_data, 2697, "Sample #2698")
safe_show(token_data, 2698, "Sample #2699")


# In[ ]:





# In[43]:


# =====================================
# Convert tokenized dataset → Pandas
# =====================================
import pandas as pd

# Inspect columns
print("Available columns:", token_data.column_names)

# Example: inspect one item
i = 0
item = token_data[i]


# In[44]:


custom_attr_name_dict = {
      # "cell_id": "obs_names",  # 👈 From AnnData.obs.index (cell barcodes)
        "cell_id": "cell_id",      
        "sample": "sample",
        "condition": "condition",
        "leiden": "leiden",
        "orig.ident": "orig.ident",  # Note: maps to column with dot
        "class": "class",
        "pre_ingest_annots": "pre_ingest_annots",
        "labels_scanvi": "labels_scanvi",
        "predictions_scanvi": "predictions_scanvi",
        "celltypes": "celltypes",
        "fine_annotation": "fine_annotation",
        "individual": "individual",
        "batch": "batch",
        "donor": "donor",
        "post_filter_annots": "post_filter_annots",
        "post_filter_annots_number": "post_filter_annots_number",
}

for k in custom_attr_name_dict.keys():
    if k in item:
        print(f"{k}:", item[k])


# In[45]:


# Select metadata columns into a DataFrame
meta_cols = custom_attr_name_dict.keys()
meta_cols = [c for c in meta_cols if c in token_data.column_names]

df_meta = token_data.select_columns(meta_cols).to_pandas()
print("\nMetadata head:\n", df_meta.head())
print("Rows:", len(df_meta))

# Convert full dataset (⚠ may be large!)
df_all = token_data.to_pandas()
print("\nFull DataFrame head:\n", df_all.head())
print("Rows:", len(df_all))


# In[46]:


# =====================================
# Save to external files
# =====================================
prefix = out_prefix if "out_prefix" in locals() else "token_data"

csv_path = f"{prefix}.all.csv"
txt_path = f"{prefix}.preview.txt"

df_all.to_csv(csv_path, index=False)

with open(txt_path, "w", encoding="utf-8") as f:
    f.write(df_all.head().to_string())
    f.write(f"\n\nrows: {len(df_all)}\n")

print("✅ Saved CSV:", csv_path)
print("✅ Saved preview TXT:", txt_path)


# Check if 'label' exists in tokenized dataset
if 'label' in token_data.column_names:
    print("✓ Found 'label' column in tokenized data")
    print("Sample values:", token_data[0]['label'])
else:
    print("✗ No 'label' column found in tokenized data")

# Check all available columns
print("Available columns:", token_data.column_names)


# In[47]:


# Claude AI : about label

# You do NOT need a "label" column. 
# The target column is specified later during training with the cell_state_dict parameter. 
# Your current custom_attr_name_dict is correct!


# In[ ]:




