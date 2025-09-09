#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# In[2]:


import sys

# ==============================
# Environment setup (if needed)
# ==============================
env_path = "/home/btanasa/miniconda3/envs/torch-cu"
os.environ["CONDA_PREFIX"] = env_path
os.environ["CONDA_DEFAULT_ENV"] = "torch-cu"
os.environ["RETICULATE_PYTHON"] = f"{env_path}/bin/python"
print("✅ Environment set to torch-cu")
print("Python executable:", sys.executable)


# In[3]:


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


# In[4]:


# Set working directory
os.chdir("/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.input")
print("Files:", os.listdir('.'))


# In[5]:


# =====================================
# Paths and Directories
# =====================================
import os
from geneformer import TranscriptomeTokenizer

# Base directories
workdir   = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/"
gf_dir    = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer"  # Geneformer dictionaries

# Input/output directories
h5ad_input_dir = os.path.join(workdir, "test_geneformer_3_full_h5ad.reprocess.sep10.input")
out_dir        = os.path.join(workdir, "test_geneformer_3_full_h5ad.reprocess.sep10.dataset")
out_prefix     = "hypoxia_all_samples_tokens"

# Dictionary files
gene_median_file     = os.path.join(gf_dir, "gene_median_dictionary_gc104M.pkl")
token_dictionary_file = os.path.join(gf_dir, "token_dictionary_gc104M.pkl")
gene_mapping_file     = os.path.join(gf_dir, "ensembl_mapping_dict_gc104M.pkl")

# Sanity check
print("Working dir:", workdir)
print("Input h5ad dir:", h5ad_input_dir)
print("Output dir:", out_dir)
print("Dictionary files:")
for f in [gene_median_file, token_dictionary_file, gene_mapping_file]:
    print(" -", f, "✓" if os.path.exists(f) else "❌ MISSING")


# In[7]:


# Previously, we did add the cell_id to the metadata : Sep9th vs Sep10th

# h5ad_file ="hypoxia_anndata_annotated.adata4_geneformer.h5ad"
# adata = sc.read_h5ad("h5ad_file")

# Add cell_id column equal to the rownames (barcodes) of .obs
# adata.obs["cell_id"] = adata.obs.index.astype(str)

# Quick check
# print(adata.obs.head()[["cell_id"]])

# write the updated AnnData into the tokenizer INPUT directory

# updated_h5ad = os.path.join(h5ad_input_dir, f"{prefix}.input_with_cell_id.h5ad")
# adata.write(updated_h5ad)

# print("Wrote updated h5ad for tokenization:", updated_h5ad)


# In[9]:


# =====================================
# Load and inspect h5ad
# =====================================

# h5ad_file="hypoxia_anndata_annotated.adata4_geneformer.h5ad" 
# the file was modified by adding cell_id

h5ad_file="hypoxia_anndata_annotated.adata5_geneformer.h5ad"
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


# In[10]:


print(prefix)


# In[11]:


print("\n.obs columns:", list(adata.obs.columns))
print(adata.obs.head(5))
print(f"\nCells: {adata.n_obs}")


# In[12]:


cols = ["class", "broad_class", "samples", "condition"]

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


# In[13]:


# =====================================
# Tokenizer (Geneformer V2 ready)
# =====================================
tokenizer = TranscriptomeTokenizer(
    custom_attr_name_dict={
        "cell_id": "obs_names",  # 👈 This preserves cell barcodes from AnnData.obs.index
        "leiden": "leiden",
        "class": "class",
        "condition": "condition",
        "broad_class": "broad_class",
    },
    gene_median_file = gene_median_file,
    token_dictionary_file = token_dictionary_file,
    gene_mapping_file = gene_mapping_file,
    collapse_gene_ids = True,
    model_input_size = 4096,
    special_token = True,
    model_version = "V2",
)

# Run tokenization
tokenizer.tokenize_data(
    data_directory = h5ad_input_dir,
    output_directory = out_dir,
    output_prefix = out_prefix,
    file_format = "h5ad"
)

print("✅ Tokenization complete:", out_prefix)


# In[14]:


import os, glob

print("Input dir exists:", os.path.isdir(h5ad_input_dir))
print("Input dir contents:", os.listdir(h5ad_input_dir) if os.path.isdir(h5ad_input_dir) else "N/A")

print("Out dir exists:", os.path.isdir(out_dir))
print("Out dir contents:", os.listdir(out_dir) if os.path.isdir(out_dir) else "N/A")

# Look for any *.dataset subfolders right under out_dir
candidates = [p for p in glob.glob(os.path.join(out_dir, "*.dataset")) if os.path.isdir(p)]
print("Direct .dataset candidates:", candidates)


# In[ ]:





# In[15]:


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


# In[16]:


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


# In[17]:


# =====================================
# Convert tokenized dataset → Pandas
# =====================================
import pandas as pd

# Inspect columns
print("Available columns:", token_data.column_names)

# Example: inspect one item
i = 0
item = token_data[i]
for k in ["cell_id", "leiden", "class", "broad_class", "condition"]:
    if k in item:
        print(f"{k}:", item[k])

# Select metadata columns into a DataFrame
meta_cols = ["cell_id", "leiden", "class", "broad_class", "condition"]
meta_cols = [c for c in meta_cols if c in token_data.column_names]

df_meta = token_data.select_columns(meta_cols).to_pandas()
print("\nMetadata head:\n", df_meta.head())
print("Rows:", len(df_meta))

# Convert full dataset (⚠ may be large!)
df_all = token_data.to_pandas()
print("\nFull DataFrame head:\n", df_all.head())
print("Rows:", len(df_all))

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


# In[18]:


# =====================================
# Save cell and gene metadata
# =====================================
import scanpy as sc
import pandas as pd

# Reload the AnnData object
adata = sc.read_h5ad(h5ad_file)

print("\n==============================")
print("📊 Final AnnData Metadata Check")
print("==============================")
print(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}\n")

# Cell metadata
obs_path = f"{prefix}.cells_metadata.csv"
adata.obs.to_csv(obs_path)
print(f"✅ Saved cell metadata (.obs) → {obs_path}")

# Gene/feature metadata
var_path = f"{prefix}.genes_metadata.csv"
adata.var.to_csv(var_path)
print(f"✅ Saved gene metadata (.var) → {var_path}")

# Also preview in console
print("\n• .obs (cell metadata) columns:")
print(list(adata.obs.columns))
print(adata.obs.head())

print("\n• .var (gene metadata) columns:")
print(list(adata.var.columns))
print(adata.var.head())


# In[19]:


# 1. hypoxia_all_samples_tokens.all.csv (~1.46M rows)
# This is the full token dataset exported from HuggingFace datasets.
# Each cell is represented not as one row, but as a sequence of tokens (gene IDs mapped to token IDs).
# If each of your ~104,439 cells gets ~14 tokens (on average), that already gives ~1.46 million rows in the .all.csv.
# This is expected, because the dataset is at the token level, not the cell level.

# 2. hypoxia_all_samples_tokens.cells_metadata.csv (104,439 rows)
# This matches exactly your number of cells.
# Each row = one cell.

# Columns = .obs metadata (sample, class, broad_class, condition, n_counts, etc.).
# ✅ This is the file you’d usually join with embeddings or labels.

# 3. hypoxia_all_samples_tokens.genes_metadata.csv (21,863 rows)
# This matches your number of genes (adata.n_vars).
# Each row = one gene.

# Columns = .var metadata (gene IDs, gene symbols, flags like mt/ribo/hb, etc.).
# ✅ This is what you’d use to check features (genes).

# ⚖️ Why the mismatch?
# Cells metadata file = one row per cell (104,439).
# Genes metadata file = one row per gene (21,863).
# Tokens file = one row per token, so it scales with both cells × expressed genes per cell. 
# That’s why it explodes in size.


# In[20]:


#!/usr/bin/env bash

# Usage:
#   ./extract_seqs.sh input.txt > sequences.txt
#   # or: cat input.txt | ./extract_seqs.sh

# Split fields at either "]," or "," and print the token right after "],"
# awk -F'],|,' '
#  NF >= 2 {
#    print $2
#  }
# ' "$@"

# ./script_extract_seqs.sh hypoxia_all_samples_tokens.all.csv | sort -u | wc -l


# In[ ]:





# In[ ]:




