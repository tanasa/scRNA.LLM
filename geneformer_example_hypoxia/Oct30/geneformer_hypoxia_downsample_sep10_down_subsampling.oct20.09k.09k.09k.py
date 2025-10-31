#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
from pathlib import Path

# Set working directory

WORKDIR = Path("")
if not WORKDIR.exists():
    raise FileNotFoundError(f"Workdir not found: {WORKDIR}")
os.chdir(WORKDIR)
print("📁 Working dir:", WORKDIR)
print("📄 Files:", ", ".join(sorted(p.name for p in WORKDIR.iterdir())))
print("Files:", os.listdir('.'))


# In[2]:


import os, sys, gc, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import torch

from datasets import load_from_disk, Dataset, DatasetDict, concatenate_datasets

# --- Warnings & GC ---
warnings.filterwarnings("ignore", category=FutureWarning)
torch.cuda.empty_cache()
gc.collect();

# --- Runtime / GPU info ---
print(f"🐍 Python: {sys.executable}")
print(f"🧠 PyTorch: {torch.__version__}")
cuda_ok = torch.cuda.is_available()
print(f"⚙️  CUDA available: {cuda_ok}")

if cuda_ok:
    dev = torch.cuda.current_device()
    name = torch.cuda.get_device_name(dev)
    alloc_gb = torch.cuda.memory_allocated() / (1024**3)
    reserv_gb = torch.cuda.memory_reserved() / (1024**3)
    print(f"🖥️  GPU[{dev}]: {name}")
    print(f"   • mem allocated: {alloc_gb:.2f} GB | reserved: {reserv_gb:.2f} GB")
    torch.cuda.empty_cache()
    gc.collect()
    alloc_gb = torch.cuda.memory_allocated() / (1024**3)
    reserv_gb = torch.cuda.memory_reserved() / (1024**3)
    print(f"   • after empty_cache → allocated: {alloc_gb:.2f} GB | reserved: {reserv_gb:.2f} GB")

# --- Optional: echo intended conda env (note: this does NOT switch the running kernel) ---
ENV_PREFIX = "/home/btanasa/miniconda3/envs/torch-cu"
os.environ["CONDA_PREFIX"] = ENV_PREFIX
os.environ["CONDA_DEFAULT_ENV"] = "torch-cu"
os.environ["RETICULATE_PYTHON"] = str(Path(ENV_PREFIX) / "bin" / "python")
print("ENV → CONDA_PREFIX:", os.environ.get("CONDA_PREFIX"))
print("ENV → CONDA_DEFAULT_ENV:", os.environ.get("CONDA_DEFAULT_ENV"))

# --- Geneformer availability (optional) ---
try:
    from geneformer import Classifier  # or TranscriptomeTokenizer
    print("✅ Geneformer import: OK")
except Exception as e:
    print("❌ Geneformer import failed:", e)


# In[ ]:





# In[3]:


from pathlib import Path
from collections import Counter
from datasets import load_from_disk, Dataset, DatasetDict, concatenate_datasets
import pandas as pd
import matplotlib.pyplot as plt

# --- Workdir + dataset path ---
WORKDIR = Path("/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.dataset")
DATASET_NAME = "hypoxia_all_samples_tokens.dataset"   

hypoxia_dataset = (WORKDIR / DATASET_NAME).resolve()
assert hypoxia_dataset.exists(), (
    f"Not found: {hypoxia_dataset}\n"
    f"Tip: available *.dataset dirs here: {[p.name for p in WORKDIR.glob('*.dataset')]}"
)

print("📁 WORKDIR:", WORKDIR)
print("📦 Dataset path:", hypoxia_dataset)

# --- Load ---
token_dataset = load_from_disk(str(hypoxia_dataset))

# Define ds_view (merge if DatasetDict)
if isinstance(token_dataset, DatasetDict):
    ds_view = concatenate_datasets(list(token_dataset.values()))
else:
    ds_view = token_dataset

# --- Basic info ---
print("Type of loaded dataset:", type(token_dataset))
print("Dataset summary:")
print(token_dataset)

# --- Features / schema ---
print("\nDataset features (schema):")
print(ds_view.features)

# --- Quick peek (first 2 records) ---
n_peek = min(2, len(ds_view))
print(f"\nFirst {n_peek} records:")
print(ds_view[:n_peek])

# --- Small pandas preview (safe) ---
HEAD_N = min(5, len(ds_view))
df_head = ds_view.select(range(HEAD_N)).to_pandas()
print("\nHead of the dataset as pandas DataFrame:")
print(df_head)

# Keep for later steps
original_dataset = token_dataset

# --- Optional: quick distributions if present ---
cols_to_check = ["class", "condition"]
for col in cols_to_check:
    if col in ds_view.column_names:
        c = Counter(ds_view[col])
        total = sum(c.values())
        print(f"\n=== {col} distribution ===")
        for k, v in c.most_common(15):
            print(f"{k}: {v} ({v/total*100:.2f}%)")

        # Small bar plot (top 15) — shown inline
        labels, values = zip(*c.most_common(15))
        plt.figure(figsize=(4, 4))
        plt.bar(range(len(values)), values)
        plt.xticks(range(len(values)), labels, rotation=90, fontsize=8)
        plt.title(f"{col} distribution")
        plt.tight_layout()
        plt.show()
    else:
        print(f"\nColumn '{col}' not found. Available: {ds_view.column_names}")

# --- Optional: input_ids shape sanity ---
if "input_ids" in ds_view.column_names and len(ds_view) > 0:
    ex0 = ds_view[0]["input_ids"]
    print(f"\nExample 'input_ids' length for row 0: {len(ex0)} (first 10): {ex0[:10]}")


# In[4]:


# === Per-condition distributions for the ORIGINAL dataset `original_dataset` ===
from datasets import Dataset, DatasetDict, concatenate_datasets
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- config ---
OUTDIR = "original_per_condition"
TOP_N  = 25        # limit to top-N classes overall for readability; set None for all
DPI    = 150
os.makedirs(OUTDIR, exist_ok=True)

# Normalize to a single Dataset
orig_ds = original_dataset
if isinstance(orig_ds, DatasetDict):
    orig_ds = concatenate_datasets(list(orig_ds.values()))
elif not isinstance(orig_ds, Dataset):
    raise TypeError(f"Unexpected type: {type(orig_ds)}")

# Build tidy DataFrame
df_orig = pd.DataFrame({"class": orig_ds["class"], "condition": orig_ds["condition"]})

# Counts: rows = class, cols = condition
counts = pd.crosstab(df_orig["class"], df_orig["condition"]).sort_index()
conditions = counts.columns.tolist()

# Choose a consistent class order (top-N by overall total)
if TOP_N is not None:
    class_order = counts.sum(axis=1).sort_values(ascending=False).head(TOP_N).index.tolist()
    counts = counts.loc[class_order]
else:
    class_order = counts.index.tolist()

# Column % (within each condition)
colpct = (counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1) * 100.0).fillna(0.0).round(2)

# --- Print tables ---
from IPython.display import display
print("Cells × Condition — COUNTS (original, limited to displayed classes)")
display(counts)
print("\nCells × Condition — COLUMN % (within condition)")
display(colpct)

# --- Plot: per-condition barh (COUNTS) ---
n_cond  = len(conditions)
n_class = len(class_order)
fig_h   = max(6, 0.28 * n_class)

fig, axes = plt.subplots(1, n_cond, figsize=(4.0 * n_cond, fig_h), sharey=True, dpi=DPI)
if n_cond == 1:
    axes = [axes]

ypos = np.arange(n_class)
for ax, cond in zip(axes, conditions):
    vals = counts[cond].reindex(class_order).values
    ax.barh(ypos, vals)
    ax.set_title(f"{cond} (n={int(vals.sum())})", fontsize=10)
    ax.set_yticks(ypos)
    ax.set_yticklabels(class_order, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Cells")
    ax.tick_params(axis="x", labelsize=8)

fig.suptitle("Original — cell type distribution per condition (counts)", fontsize=12, y=1.02)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(OUTDIR, "per_condition_counts_barh.png"), dpi=DPI, bbox_inches="tight")
plt.close(fig)

# --- Plot: per-condition barh (COLUMN %) ---
fig, axes = plt.subplots(1, n_cond, figsize=(4.0 * n_cond, fig_h), sharey=True, dpi=DPI)
if n_cond == 1:
    axes = [axes]

for ax, cond in zip(axes, conditions):
    vals = colpct[cond].reindex(class_order).values
    ax.barh(ypos, vals)
    ax.set_title(f"{cond} (% within condition)", fontsize=10)
    ax.set_yticks(ypos)
    ax.set_yticklabels(class_order, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("% of cells in condition")
    ax.tick_params(axis="x", labelsize=8)

fig.suptitle("Original — cell type distribution per condition (percent)", fontsize=12, y=1.02)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(OUTDIR, "per_condition_colpct_barh.png"), dpi=DPI, bbox_inches="tight")
plt.close(fig)

print(f"Saved figures to: {OUTDIR}/")


# In[5]:


print("Which cell types are in percentage higher than 5% ?")


# In[6]:


# === Per-condition distributions for the ORIGINAL dataset `original_dataset` ===
from datasets import Dataset, DatasetDict, concatenate_datasets
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- config ---
OUTDIR = "original_per_condition"
TOP_N  = 25        # limit to top-N classes overall for readability; set None for all
DPI    = 150
MIN_PERCENTAGE = 5.0  # 5% threshold per condition
os.makedirs(OUTDIR, exist_ok=True)

# Normalize to a single Dataset
orig_ds = original_dataset
if isinstance(orig_ds, DatasetDict):
    orig_ds = concatenate_datasets(list(orig_ds.values()))
elif not isinstance(orig_ds, Dataset):
    raise TypeError(f"Unexpected type: {type(orig_ds)}")

# Build tidy DataFrame
df_orig = pd.DataFrame({"class": orig_ds["class"], "condition": orig_ds["condition"]})

# Counts: rows = class, cols = condition
counts = pd.crosstab(df_orig["class"], df_orig["condition"]).sort_index()
conditions = counts.columns.tolist()

# Choose a consistent class order (top-N by overall total)
if TOP_N is not None:
    class_order = counts.sum(axis=1).sort_values(ascending=False).head(TOP_N).index.tolist()
    counts = counts.loc[class_order]
else:
    class_order = counts.index.tolist()

# Column % (within each condition)
colpct = (counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1) * 100.0).fillna(0.0).round(2)

# =============================================================================
# 5% THRESHOLD ANALYSIS
# =============================================================================
print("="*80)
print("5% THRESHOLD ANALYSIS PER CONDITION")
print("="*80)

# Check which cell types meet 5% threshold in each condition
meets_threshold = colpct >= MIN_PERCENTAGE

print(f"Cell types meeting {MIN_PERCENTAGE}% threshold in each condition:")
print("-" * 60)

for condition in conditions:
    condition_cell_types = meets_threshold[condition]
    qualifying_types = condition_cell_types[condition_cell_types].index.tolist()
    print(f"\n{condition}:")
    print(f"  Total cell types: {len(condition_cell_types)}")
    print(f"  Meeting {MIN_PERCENTAGE}% threshold: {len(qualifying_types)}")
    print(f"  Qualifying types: {qualifying_types}")
    
    # Show percentages for qualifying types
    for cell_type in qualifying_types:
        pct = colpct.loc[cell_type, condition]
        print(f"    - {cell_type}: {pct:.2f}%")

# Find cell types that meet 5% threshold in ALL conditions
print(f"\n" + "="*60)
print("CELL TYPES MEETING 5% THRESHOLD IN ALL CONDITIONS")
print("="*60)

all_conditions_meet = meets_threshold.all(axis=1)
qualifying_all_conditions = all_conditions_meet[all_conditions_meet].index.tolist()

print(f"Cell types meeting {MIN_PERCENTAGE}% in ALL conditions: {len(qualifying_all_conditions)}")
print(f"Qualifying types: {qualifying_all_conditions}")

# Show detailed breakdown for each qualifying type
print(f"\nDetailed breakdown for qualifying types:")
for cell_type in qualifying_all_conditions:
    print(f"\n{cell_type}:")
    for condition in conditions:
        pct = colpct.loc[cell_type, condition]
        status = "✓" if pct >= MIN_PERCENTAGE else "✗"
        print(f"  - {condition}: {pct:.2f}% {status}")

# Find cell types that meet 5% threshold in ANY condition
print(f"\n" + "="*60)
print("CELL TYPES MEETING 5% THRESHOLD IN ANY CONDITION")
print("="*60)

any_condition_meets = meets_threshold.any(axis=1)
qualifying_any_condition = any_condition_meets[any_condition_meets].index.tolist()

print(f"Cell types meeting {MIN_PERCENTAGE}% in ANY condition: {len(qualifying_any_condition)}")
print(f"Qualifying types: {qualifying_any_condition}")

# Show which conditions each type qualifies for
print(f"\nDetailed breakdown for any-condition qualifying types:")
for cell_type in qualifying_any_condition:
    print(f"\n{cell_type}:")
    qualifying_conditions = []
    for condition in conditions:
        pct = colpct.loc[cell_type, condition]
        if pct >= MIN_PERCENTAGE:
            qualifying_conditions.append(condition)
            print(f"  - {condition}: {pct:.2f}% ✓")
        else:
            print(f"  - {condition}: {pct:.2f}% ✗")
    print(f"  → Qualifies in: {qualifying_conditions}")


# In[7]:


print("Which cell types are in percentage higher than 1% ?")


# In[8]:


# === Per-condition distributions for the ORIGINAL dataset `original_dataset` ===
from datasets import Dataset, DatasetDict, concatenate_datasets
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# --- config ---
OUTDIR = "original_per_condition"
TOP_N  = 25        # limit to top-N classes overall for readability; set None for all
DPI    = 150
MIN_PERCENTAGE = 1.0  # 1% threshold per condition
os.makedirs(OUTDIR, exist_ok=True)

# Normalize to a single Dataset
orig_ds = original_dataset
if isinstance(orig_ds, DatasetDict):
    orig_ds = concatenate_datasets(list(orig_ds.values()))
elif not isinstance(orig_ds, Dataset):
    raise TypeError(f"Unexpected type: {type(orig_ds)}")

# Build tidy DataFrame
df_orig = pd.DataFrame({"class": orig_ds["class"], "condition": orig_ds["condition"]})

# Counts: rows = class, cols = condition
counts = pd.crosstab(df_orig["class"], df_orig["condition"]).sort_index()
conditions = counts.columns.tolist()

# Choose a consistent class order (top-N by overall total)
if TOP_N is not None:
    class_order = counts.sum(axis=1).sort_values(ascending=False).head(TOP_N).index.tolist()
    counts = counts.loc[class_order]
else:
    class_order = counts.index.tolist()

# Column % (within each condition)
colpct = (counts.div(counts.sum(axis=0).replace(0, np.nan), axis=1) * 100.0).fillna(0.0).round(2)

# =============================================================================
# 1% THRESHOLD ANALYSIS
# =============================================================================
print("="*80)
print("1% THRESHOLD ANALYSIS PER CONDITION")
print("="*80)

# Check which cell types meet 1% threshold in each condition
meets_threshold = colpct >= MIN_PERCENTAGE

print(f"Cell types meeting {MIN_PERCENTAGE}% threshold in each condition:")
print("-" * 60)

for condition in conditions:
    condition_cell_types = meets_threshold[condition]
    qualifying_types = condition_cell_types[condition_cell_types].index.tolist()
    print(f"\n{condition}:")
    print(f"  Total cell types: {len(condition_cell_types)}")
    print(f"  Meeting {MIN_PERCENTAGE}% threshold: {len(qualifying_types)}")
    print(f"  Qualifying types: {qualifying_types}")
    
    # Show percentages for qualifying types
    for cell_type in qualifying_types:
        pct = colpct.loc[cell_type, condition]
        print(f"    - {cell_type}: {pct:.2f}%")

# Find cell types that meet 1% threshold in ALL conditions
print(f"\n" + "="*60)
print("CELL TYPES MEETING 1% THRESHOLD IN ALL CONDITIONS")
print("="*60)

all_conditions_meet = meets_threshold.all(axis=1)
qualifying_all_conditions = all_conditions_meet[all_conditions_meet].index.tolist()

print(f"Cell types meeting {MIN_PERCENTAGE}% in ALL conditions: {len(qualifying_all_conditions)}")
print(f"Qualifying types: {qualifying_all_conditions}")

# Show detailed breakdown for each qualifying type
print(f"\nDetailed breakdown for qualifying types:")
for cell_type in qualifying_all_conditions:
    print(f"\n{cell_type}:")
    for condition in conditions:
        pct = colpct.loc[cell_type, condition]
        status = "✓" if pct >= MIN_PERCENTAGE else "✗"
        print(f"  - {condition}: {pct:.2f}% {status}")

# Find cell types that meet 1% threshold in ANY condition
print(f"\n" + "="*60)
print("CELL TYPES MEETING 1% THRESHOLD IN ANY CONDITION")
print("="*60)

any_condition_meets = meets_threshold.any(axis=1)
qualifying_any_condition = any_condition_meets[any_condition_meets].index.tolist()

print(f"Cell types meeting {MIN_PERCENTAGE}% in ANY condition: {len(qualifying_any_condition)}")
print(f"Qualifying types: {qualifying_any_condition}")

# Show which conditions each type qualifies for
print(f"\nDetailed breakdown for any-condition qualifying types:")
for cell_type in qualifying_any_condition:
    print(f"\n{cell_type}:")
    qualifying_conditions = []
    for condition in conditions:
        pct = colpct.loc[cell_type, condition]
        if pct >= MIN_PERCENTAGE:
            qualifying_conditions.append(condition)
            print(f"  - {condition}: {pct:.2f}% ✓")
        else:
            print(f"  - {condition}: {pct:.2f}% ✗")
    print(f"  → Qualifies in: {qualifying_conditions}")

# Summary statistics
print(f"\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)

total_cell_types = len(counts)
all_conditions_count = len(qualifying_all_conditions)
any_condition_count = len(qualifying_any_condition)

print(f"Total cell types in dataset: {total_cell_types}")
print(f"Cell types meeting {MIN_PERCENTAGE}% in ALL conditions: {all_conditions_count} ({all_conditions_count/total_cell_types*100:.1f}%)")
print(f"Cell types meeting {MIN_PERCENTAGE}% in ANY condition: {any_condition_count} ({any_condition_count/total_cell_types*100:.1f}%)")
print(f"Cell types NOT meeting {MIN_PERCENTAGE}% in ANY condition: {total_cell_types - any_condition_count} ({(total_cell_types - any_condition_count)/total_cell_types*100:.1f}%)")

# Comparison with 5% threshold
print(f"\n" + "="*60)
print("COMPARISON: 1% vs 5% THRESHOLD")
print("="*60)

# Calculate 5% threshold for comparison
meets_5pct = colpct >= 5.0
all_conditions_5pct = meets_5pct.all(axis=1)
qualifying_all_5pct = all_conditions_5pct[all_conditions_5pct].index.tolist()
any_condition_5pct = meets_5pct.any(axis=1)
qualifying_any_5pct = any_condition_5pct[any_condition_5pct].index.tolist()

print(f"5% threshold - ALL conditions: {len(qualifying_all_5pct)} cell types")
print(f"1% threshold - ALL conditions: {len(qualifying_all_conditions)} cell types")
print(f"Difference (1% vs 5%): +{len(qualifying_all_conditions) - len(qualifying_all_5pct)} cell types")

print(f"\n5% threshold - ANY condition: {len(qualifying_any_5pct)} cell types")
print(f"1% threshold - ANY condition: {len(qualifying_any_condition)} cell types")
print(f"Difference (1% vs 5%): +{len(qualifying_any_condition) - len(qualifying_any_5pct)} cell types")

# Show which cell types are gained by lowering threshold to 1%
gained_all_conditions = set(qualifying_all_conditions) - set(qualifying_all_5pct)
gained_any_condition = set(qualifying_any_condition) - set(qualifying_any_5pct)

if gained_all_conditions:
    print(f"\nCell types gained by lowering threshold to 1% (ALL conditions):")
    for cell_type in sorted(gained_all_conditions):
        print(f"  - {cell_type}")

if gained_any_condition:
    print(f"\nCell types gained by lowering threshold to 1% (ANY condition):")
    for cell_type in sorted(gained_any_condition):
        print(f"  - {cell_type}")



# In[9]:


print("""
Types of Stratification:

1. Single Stratification:
   One variable (e.g., cell type only)
   Simple proportional sampling

2. Joint Stratification:
   Multiple variables (e.g., cell type + condition)
   Preserves interactions between variables

3. Balanced Stratification:
   EQUAL09k representation across groups
   Ignores original proportions

""")


# In[10]:


import os
print(os.getcwd())


# In[11]:


# =============================================================================
# PROPORTIONAL SUBSAMPLING: 9000 CELLS PER CONDITION
# =============================================================================

from datasets import Dataset, DatasetDict, concatenate_datasets
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
import math, random, os, shutil

# ========= CONFIG =========
CELLS_PER_CONDITION = 9000     # 9000 cells per condition
TOTAL_TARGET = 27000           # 9000 * 3 conditions
SEED = 42
SAVE_PATH = "./hypoxia_EQUAL09k_9000_per_condition.dataset"
NUM_PROC = 8
# ==========================

print("="*80)
print("STRATIFIED EQUAL SUBSAMPLING: 9000 CELLS PER CONDITION")
print("="*80)
print(f"Target: {CELLS_PER_CONDITION} cells per condition")
print(f"Total target: {TOTAL_TARGET} cells")

# 0) Normalize to a single Dataset
ds = original_dataset
if isinstance(ds, DatasetDict):
    ds = concatenate_datasets(list(ds.values()))
elif not isinstance(ds, Dataset):
    raise TypeError(f"Unexpected type: {type(ds)}")

# 1) Sanity checks
for col in ["class", "condition"]:
    if col not in ds.column_names:
        raise KeyError(f"Missing column '{col}'. Available: {ds.column_names}")

print(f"Start: {len(ds):,} cells")

# 2) Get condition data and calculate percentages
# FIX: Use pandas to get unique values
df_temp = pd.DataFrame({"condition": ds["condition"]})
conditions = sorted(df_temp["condition"].unique())
print(f"Conditions: {conditions}")

# Calculate percentages for each condition
condition_percentages = {}
for condition in conditions:
    condition_data = ds.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    
    # Create cross-tabulation
    df_condition = pd.DataFrame({"class": condition_data["class"], "condition": condition_data["condition"]})
    crosstab = pd.crosstab(df_condition["class"], df_condition["condition"])
    
    # Calculate percentages within this condition
    percentages = (crosstab / crosstab.sum(axis=0)) * 100
    condition_percentages[condition] = percentages[condition].to_dict()
    
    print(f"\n{condition} - Cell type percentages:")
    for cell_type, pct in sorted(condition_percentages[condition].items(), key=lambda x: x[1], reverse=True):
        print(f"  {cell_type}: {pct:.2f}%")

# 3) Subsampling per condition
print(f"\n" + "="*60)
print("STRATIFIED EQUAL SUBSAMPLING PER CONDITION")
print("="*60)

sampled_parts = []

for condition in conditions:
    print(f"\nProcessing condition: {condition}")
    
    # Filter data for this condition
    condition_data = ds.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    available_cells = len(condition_data)
    
    print(f"  Available cells: {available_cells:,}")
    
    if available_cells == 0:
        print(f"  ⚠️ No cells found for condition {condition}")
        continue
    
    # Get percentages for this condition
    percentages = condition_percentages[condition]
    
    # Calculate target cells per cell type (proportional)
    cell_type_targets = {}
    total_target_cells = 0
    
    print(f"  Target cells per cell type (proportional):")
    for cell_type, percentage in percentages.items():
        target_cells_for_type = int(CELLS_PER_CONDITION * (percentage / 100))
        cell_type_targets[cell_type] = target_cells_for_type
        total_target_cells += target_cells_for_type
        
        if target_cells_for_type > 0:
            print(f"    {cell_type}: {target_cells_for_type} cells ({percentage:.2f}%)")
    
    # Adjust for rounding errors to get exactly 9000 cells
    if total_target_cells != CELLS_PER_CONDITION:
        difference = CELLS_PER_CONDITION - total_target_cells
        print(f"  Adjusting for rounding: {difference} cells")
        
        # Add/subtract from the largest cell type
        largest_type = max(cell_type_targets.items(), key=lambda x: x[1])
        cell_type_targets[largest_type[0]] += difference
        print(f"    Adjusted {largest_type[0]}: {cell_type_targets[largest_type[0]]} cells")
    
    # Sample from each cell type within this condition
    condition_parts = []
    actual_sampled = 0
    
    for cell_type, target_count in cell_type_targets.items():
        if target_count <= 0:
            continue
            
        # Filter for this cell type within this condition
        cell_type_data = condition_data.filter(
            lambda ex, cls=cell_type: ex["class"] == cls, 
            num_proc=NUM_PROC
        )
        
        if len(cell_type_data) > 0:
            # Sample the target number
            actual_sample = min(target_count, len(cell_type_data))
            sampled_class = cell_type_data.shuffle(seed=SEED).select(range(actual_sample))
            condition_parts.append(sampled_class)
            actual_sampled += actual_sample
            
            print(f"    Sampled {actual_sample:,} cells from {cell_type}")
        else:
            print(f"    ⚠️ No cells found for {cell_type} in {condition}")
    
    # Combine all cell types for this condition
    if condition_parts:
        condition_sampled = concatenate_datasets(condition_parts).shuffle(seed=SEED)
        sampled_parts.append(condition_sampled)
        print(f"  ✅ Total sampled from {condition}: {len(condition_sampled):,} cells")
    else:
        print(f"  ⚠️ No cells sampled from {condition}")

# 4) Combine all sampled data
if sampled_parts:
    ds_sampled = concatenate_datasets(sampled_parts).shuffle(seed=SEED)
    print(f"\n✅ Final sampled dataset: {len(ds_sampled):,} cells")
else:
    raise ValueError("No cells were sampled!")

# 5) Verification of sampling
print(f"\n" + "="*60)
print("VERIFICATION OF STRATIFIED SUBSAMPLING")
print("="*60)

# Check cells per condition
sampled_condition_counts = Counter(ds_sampled["condition"])
print(f"Cells per condition after sampling:")
for condition in conditions:
    count = sampled_condition_counts.get(condition, 0)
    print(f"  {condition}: {count:,} cells")

# Check cell type distribution per condition
print(f"\nCell type distribution per condition (percentages):")
for condition in conditions:
    print(f"\n{condition}:")
    condition_sampled = ds_sampled.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    
    if len(condition_sampled) > 0:
        # Calculate percentages
        df_sampled = pd.DataFrame({"class": condition_sampled["class"], "condition": condition_sampled["condition"]})
        crosstab_sampled = pd.crosstab(df_sampled["class"], df_sampled["condition"])
        percentages_sampled = (crosstab_sampled / crosstab_sampled.sum(axis=0)) * 100
        
        # Show top cell types
        for cell_type, pct in sorted(percentages_sampled[condition].items(), key=lambda x: x[1], reverse=True):
            if pct > 0.1:  # Only show types with >0.1%
                print(f"  {cell_type}: {pct:.2f}%")

# 6) Create cross-tabulation
df_sampled = pd.DataFrame({
    "class": ds_sampled["class"], 
    "condition": ds_sampled["condition"]
})

ct_counts = pd.crosstab(df_sampled["class"], df_sampled["condition"]).sort_index()
ct_rowpct = (ct_counts.div(ct_counts.sum(axis=1).replace(0, np.nan), axis=0) * 100).fillna(0).round(2)

print(f"\nCross-tabulation (counts):")
print(ct_counts)

print(f"\nCross-tabulation (row percentages):")
print(ct_rowpct)

# 7) Save the sampled dataset
if SAVE_PATH:
    if os.path.exists(SAVE_PATH):
        shutil.rmtree(SAVE_PATH)
    ds_sampled.save_to_disk(SAVE_PATH)
    print(f"\n✅ Saved sampled dataset to: {SAVE_PATH}")

# 8) Summary statistics
print(f"\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)
print(f"Original dataset: {len(ds):,} cells")
print(f"Sampled dataset: {len(ds_sampled):,} cells")
print(f"Reduction: {((len(ds) - len(ds_sampled)) / len(ds) * 100):.1f}%")
print(f"Target achieved: {len(ds_sampled):,} / {TOTAL_TARGET} cells")

# Show final distribution
print(f"\nFinal distribution per condition:")
for condition in conditions:
    count = sampled_condition_counts.get(condition, 0)
    pct = (count / len(ds_sampled)) * 100
    print(f"  {condition}: {count:,} cells ({pct:.1f}%)")

print(f"\n" + "="*60)
print("STRATIFIED SUBSAMPLING COMPLETE")
print("="*60)


# In[12]:


# Row % within class
# What it answers: “Within this cell class, how are its cells split across conditions?”

# Column % within condition
# What it answers: “Within this condition, how are its cells split across classes?”

# Condition-only stratification: you down-sample by columns (conditions) only. Column proportions are preserved; row/class proportions may shift.
# Class-only stratification: you down-sample by rows (classes) only. Row proportions are preserved; column/condition proportions may shift.
# Joint (class × condition) stratification (what we implemented): you down-sample by cells of the contingency table (each class–condition pair). 
# This effectively preserves both row and column distributions (marginals) as a consequence of preserving the joint distribution—allowing 
# for tiny rounding differences and any buckets with too few cells (which can end up at 0 after downsampling)


# In[13]:


print(OUTDIR)


# In[19]:


# Add this section after the sampling is complete:

# =============================================================================
# COMPARISON: BEFORE vs AFTER EQUAL SUBSAMPLING
# =============================================================================

print(f"\n" + "="*80)
print("COMPARISON: BEFORE vs AFTER EQUAL SUBSAMPLING")
print("="*80)

# Create comparison tables
comparison_data = []

for condition in conditions:
    print(f"\n{'='*60}")
    print(f"CONDITION: {condition}")
    print(f"{'='*60}")
    
    # Original data for this condition
    orig_condition = ds.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    orig_df = pd.DataFrame({"class": orig_condition["class"], "condition": orig_condition["condition"]})
    orig_crosstab = pd.crosstab(orig_df["class"], orig_df["condition"])
    orig_percentages = (orig_crosstab / orig_crosstab.sum(axis=0)) * 100
    orig_pct = orig_percentages[condition].to_dict()
    
    # Sampled data for this condition
    sampled_condition = ds_sampled.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    sampled_df = pd.DataFrame({"class": sampled_condition["class"], "condition": sampled_condition["condition"]})
    sampled_crosstab = pd.crosstab(sampled_df["class"], sampled_df["condition"])
    sampled_percentages = (sampled_crosstab / sampled_crosstab.sum(axis=0)) * 100
    sampled_pct = sampled_percentages[condition].to_dict()
    
    # Get all unique cell types
    all_cell_types = set(orig_pct.keys()) | set(sampled_pct.keys())
    
    print(f"\nOriginal dataset: {len(orig_condition):,} cells")
    print(f"Sampled dataset: {len(sampled_condition):,} cells")
    print(f"Reduction: {((len(orig_condition) - len(sampled_condition)) / len(orig_condition) * 100):.1f}%")
    
    print(f"\n{'Cell Type':<25} {'Original %':<12} {'Sampled %':<12} {'Difference':<12} {'Status':<10}")
    print("-" * 80)
    
    for cell_type in sorted(all_cell_types):
        orig_val = orig_pct.get(cell_type, 0)
        sampled_val = sampled_pct.get(cell_type, 0)
        difference = sampled_val - orig_val
        
        # Determine status
        if abs(difference) < 0.5:
            status = "✓ Good"
        elif abs(difference) < 1.0:
            status = "⚠ Fair"
        else:
            status = "✗ Poor"
        
        print(f"{cell_type:<25} {orig_val:>10.2f}% {sampled_val:>10.2f}% {difference:>+10.2f}% {status:<10}")
        
        # Store for summary
        comparison_data.append({
            'condition': condition,
            'cell_type': cell_type,
            'original_pct': orig_val,
            'sampled_pct': sampled_val,
            'difference': difference,
            'status': status
        })

# Create summary DataFrame
comparison_df = pd.DataFrame(comparison_data)

# Save comparison to CSV
comparison_df.to_csv("EQUAL09k_subsampling_comparison.csv", index=False)
print(f"\n✅ Comparison data saved to: EQUAL09k_subsampling_comparison.csv")


# In[15]:


# =============================================================================
# BARPLOTS: BEFORE vs AFTER EQUAL09k SUBSAMPLING
# =============================================================================

import matplotlib.pyplot as plt
import seaborn as sns

# Set style
plt.style.use('default')
sns.set_palette("husl")

# Create comparison plots
fig, axes = plt.subplots(2, len(conditions), figsize=(6*len(conditions), 10))
if len(conditions) == 1:
    axes = axes.reshape(2, 1)

fig.suptitle('Cell Type Percentages: Before vs After EQUAL09k Downsampling', fontsize=16, fontweight='bold')

for i, condition in enumerate(conditions):
    # Get data for this condition
    orig_condition = ds.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    sampled_condition = ds_sampled.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    
    # Calculate percentages
    orig_df = pd.DataFrame({"class": orig_condition["class"], "condition": orig_condition["condition"]})
    orig_crosstab = pd.crosstab(orig_df["class"], orig_df["condition"])
    orig_percentages = (orig_crosstab / orig_crosstab.sum(axis=0)) * 100
    orig_pct = orig_percentages[condition].to_dict()
    
    sampled_df = pd.DataFrame({"class": sampled_condition["class"], "condition": sampled_condition["condition"]})
    sampled_crosstab = pd.crosstab(sampled_df["class"], sampled_df["condition"])
    sampled_percentages = (sampled_crosstab / sampled_crosstab.sum(axis=0)) * 100
    sampled_pct = sampled_percentages[condition].to_dict()
    
    # Get all cell types and sort by original percentage
    all_cell_types = sorted(set(orig_pct.keys()) | set(sampled_pct.keys()), 
                           key=lambda x: orig_pct.get(x, 0), reverse=True)
    
    # Prepare data for plotting
    orig_values = [orig_pct.get(ct, 0) for ct in all_cell_types]
    sampled_values = [sampled_pct.get(ct, 0) for ct in all_cell_types]
    
    # Plot 1: Original percentages
    axes[0, i].bar(range(len(all_cell_types)), orig_values, alpha=0.7, color='skyblue', edgecolor='navy')
    axes[0, i].set_title(f'{condition} - Original\n({len(orig_condition):,} cells)', fontweight='bold')
    axes[0, i].set_ylabel('Percentage (%)')
    axes[0, i].set_xticks(range(len(all_cell_types)))
    axes[0, i].set_xticklabels(all_cell_types, rotation=45, ha='right')
    axes[0, i].grid(True, alpha=0.3)
    
    # Add value labels on bars
    for j, v in enumerate(orig_values):
        if v > 1:  # Only label bars > 1%
            axes[0, i].text(j, v + 0.5, f'{v:.1f}%', ha='center', va='bottom', fontsize=8)
    
    # Plot 2: Sampled percentages
    axes[1, i].bar(range(len(all_cell_types)), sampled_values, alpha=0.7, color='lightcoral', edgecolor='darkred')
    axes[1, i].set_title(f'{condition} - After Subsampling\n({len(sampled_condition):,} cells)', fontweight='bold')
    axes[1, i].set_ylabel('Percentage (%)')
    axes[1, i].set_xlabel('Cell Type')
    axes[1, i].set_xticks(range(len(all_cell_types)))
    axes[1, i].set_xticklabels(all_cell_types, rotation=45, ha='right')
    axes[1, i].grid(True, alpha=0.3)
    
    # Add value labels on bars
    for j, v in enumerate(sampled_values):
        if v > 1:  # Only label bars > 1%
            axes[1, i].text(j, v + 0.5, f'{v:.1f}%', ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig('EQUAL09k_subsampling_barplots.png', dpi=300, bbox_inches='tight')
plt.show()

# =============================================================================
# SIDE-BY-SIDE COMPARISON BARPLOTS
# =============================================================================

# Create side-by-side comparison plots
fig, axes = plt.subplots(1, len(conditions), figsize=(8*len(conditions), 6))
if len(conditions) == 1:
    axes = [axes]

fig.suptitle('Cell Type Percentages: Side-by-Side Comparison', fontsize=16, fontweight='bold')

for i, condition in enumerate(conditions):
    # Get data for this condition
    orig_condition = ds.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    sampled_condition = ds_sampled.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    
    # Calculate percentages
    orig_df = pd.DataFrame({"class": orig_condition["class"], "condition": orig_condition["condition"]})
    orig_crosstab = pd.crosstab(orig_df["class"], orig_df["condition"])
    orig_percentages = (orig_crosstab / orig_crosstab.sum(axis=0)) * 100
    orig_pct = orig_percentages[condition].to_dict()
    
    sampled_df = pd.DataFrame({"class": sampled_condition["class"], "condition": sampled_condition["condition"]})
    sampled_crosstab = pd.crosstab(sampled_df["class"], sampled_df["condition"])
    sampled_percentages = (sampled_crosstab / sampled_crosstab.sum(axis=0)) * 100
    sampled_pct = sampled_percentages[condition].to_dict()
    
    # Get all cell types and sort by original percentage
    all_cell_types = sorted(set(orig_pct.keys()) | set(sampled_pct.keys()), 
                           key=lambda x: orig_pct.get(x, 0), reverse=True)
    
    # Prepare data for plotting
    orig_values = [orig_pct.get(ct, 0) for ct in all_cell_types]
    sampled_values = [sampled_pct.get(ct, 0) for ct in all_cell_types]
    
    # Create side-by-side bars
    x = np.arange(len(all_cell_types))
    width = 0.35
    
    bars1 = axes[i].bar(x - width/2, orig_values, width, label='Original', alpha=0.7, color='skyblue', edgecolor='navy')
    bars2 = axes[i].bar(x + width/2, sampled_values, width, label='After Subsampling', alpha=0.7, color='lightcoral', edgecolor='darkred')
    
    axes[i].set_title(f'{condition}\nOriginal: {len(orig_condition):,} cells → Sampled: {len(sampled_condition):,} cells', fontweight='bold')
    axes[i].set_ylabel('Percentage (%)')
    axes[i].set_xlabel('Cell Type')
    axes[i].set_xticks(x)
    axes[i].set_xticklabels(all_cell_types, rotation=45, ha='right')
    axes[i].legend()
    axes[i].grid(True, alpha=0.3)
    
    # Add value labels on bars
    for j, (v1, v2) in enumerate(zip(orig_values, sampled_values)):
        if v1 > 1:  # Only label bars > 1%
            axes[i].text(j - width/2, v1 + 0.5, f'{v1:.1f}%', ha='center', va='bottom', fontsize=8)
        if v2 > 1:
            axes[i].text(j + width/2, v2 + 0.5, f'{v2:.1f}%', ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig('EQUAL09k_subsampling_side_by_side.png', dpi=300, bbox_inches='tight')
plt.show()

# =============================================================================
# DIFFERENCE BARPLOT
# =============================================================================

# Create difference plot
fig, axes = plt.subplots(1, len(conditions), figsize=(8*len(conditions), 6))
if len(conditions) == 1:
    axes = [axes]

fig.suptitle('Cell Type Percentage Differences: After - Before Downsampling', fontsize=16, fontweight='bold')

for i, condition in enumerate(conditions):
    # Get data for this condition
    orig_condition = ds.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    sampled_condition = ds_sampled.filter(lambda ex, cond=condition: ex["condition"] == cond, num_proc=NUM_PROC)
    
    # Calculate percentages
    orig_df = pd.DataFrame({"class": orig_condition["class"], "condition": orig_condition["condition"]})
    orig_crosstab = pd.crosstab(orig_df["class"], orig_df["condition"])
    orig_percentages = (orig_crosstab / orig_crosstab.sum(axis=0)) * 100
    orig_pct = orig_percentages[condition].to_dict()
    
    sampled_df = pd.DataFrame({"class": sampled_condition["class"], "condition": sampled_condition["condition"]})
    sampled_crosstab = pd.crosstab(sampled_df["class"], sampled_df["condition"])
    sampled_percentages = (sampled_crosstab / sampled_crosstab.sum(axis=0)) * 100
    sampled_pct = sampled_percentages[condition].to_dict()
    
    # Get all cell types and sort by original percentage
    all_cell_types = sorted(set(orig_pct.keys()) | set(sampled_pct.keys()), 
                           key=lambda x: orig_pct.get(x, 0), reverse=True)
    
    # Calculate differences
    differences = [sampled_pct.get(ct, 0) - orig_pct.get(ct, 0) for ct in all_cell_types]
    
    # Color bars based on difference
    colors = ['red' if d < -0.5 else 'orange' if d < 0.5 else 'green' for d in differences]
    
    bars = axes[i].bar(range(len(all_cell_types)), differences, color=colors, alpha=0.7, edgecolor='black')
    axes[i].set_title(f'{condition} - Percentage Differences', fontweight='bold')
    axes[i].set_ylabel('Difference (%)')
    axes[i].set_xlabel('Cell Type')
    axes[i].set_xticks(range(len(all_cell_types)))
    axes[i].set_xticklabels(all_cell_types, rotation=45, ha='right')
    axes[i].axhline(y=0, color='black', linestyle='-', alpha=0.3)
    axes[i].grid(True, alpha=0.3)
    
    # Add value labels on bars
    for j, v in enumerate(differences):
        if abs(v) > 0.1:  # Only label significant differences
            axes[i].text(j, v + (0.1 if v >= 0 else -0.1), f'{v:+.1f}%', ha='center', va='bottom' if v >= 0 else 'top', fontsize=8)
    
    # Add legend
    axes[i].text(0.02, 0.98, 'Red: < -0.5%\nOrange: -0.5% to +0.5%\nGreen: > +0.5%', 
                transform=axes[i].transAxes, verticalalignment='top', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig('EQUAL09k_subsampling_differences.png', dpi=300, bbox_inches='tight')
plt.show()

print("✅ Barplots saved:")
print("  - EQUAL09k_subsampling_barplots.png")
print("  - EQUAL09k_subsampling_side_by_side.png") 
print("  - EQUAL09k_subsampling_differences.png")


# In[16]:


# In the currant implementation, the EQUAL09k sampling does NOT filter cells based on a 1% threshold.


# In[17]:


# The resulting folder is : hypoxia_EQUAL09k_9000_per_condition.dataset
# /mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.dataset


# In[ ]:




