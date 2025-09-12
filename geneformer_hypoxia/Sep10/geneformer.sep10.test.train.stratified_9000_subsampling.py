#!/usr/bin/env python
# coding: utf-8

# In[1]:


# 🧱 Standard library
import os
import sys
import subprocess
import datetime
import warnings
import gc
import glob
import pickle
from collections import Counter, defaultdict
from pathlib import Path

# 📦 Third-party libraries
import numpy as np
import pandas as pd
import seaborn as sns
import torch

# 📚 Machine learning & NLP
from sklearn.metrics import accuracy_score, f1_score
from sklearn.model_selection import train_test_split
from transformers import BertForSequenceClassification, Trainer
from transformers.training_args import TrainingArguments

# 🧬 Geneformer & HuggingFace datasets
from datasets import Dataset, DatasetDict, load_from_disk
from geneformer import (
    Classifier,
    DataCollatorForCellClassification,
    EmbExtractor
)


# In[2]:


import os
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Set working directory
# os.chdir("/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset")
os.chdir("/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.stratified.sampling.9000.dataset/")

print("Files:", os.listdir('.'))


# In[3]:


# ==============================
# Environment setup (if needed)
# ==============================
env_path = "/home/btanasa/miniconda3/envs/torch-cu"
os.environ["CONDA_PREFIX"] = env_path
os.environ["CONDA_DEFAULT_ENV"] = "torch-cu"
os.environ["RETICULATE_PYTHON"] = f"{env_path}/bin/python"
print("✅ Environment set to torch-cu")
print("Python executable:", sys.executable)


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

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Clear GPU cache
torch.cuda.empty_cache()
gc.collect()

# Check GPU memory
print(f"GPU memory allocated: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"GPU memory reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")

import torch
torch.cuda.empty_cache()
print(f"CUDA memory after empty_cache: {torch.cuda.memory_allocated() / (1024**3):.2f} GB")


import os
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

print("✅ Environment variable set!")
print(f"PYTORCH_CUDA_ALLOC_CONF: {os.environ.get('PYTORCH_CUDA_ALLOC_CONF')}")


# In[4]:


import torch
import gc

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# Clear GPU cache
torch.cuda.empty_cache()
gc.collect()

# Check GPU memory
print(f"GPU memory allocated: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"GPU memory reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")

import torch
torch.cuda.empty_cache()
print(f"CUDA memory after empty_cache: {torch.cuda.memory_allocated() / (1024**3):.2f} GB")


# In[5]:


# mkdir /mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.stratified.sampling.9000.dataset.out

token_dictionary = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl"
model_v1 = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/Geneformer-V1-10M/"
model_v2 = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/Geneformer-V2-104M"

# 🚀 OPTIMIZED CONFIGURATION FOR HIGHER ACCURACY
CONFIG = {
    
    # "dataset_path": "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset",]
    "dataset_path": "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.stratified.sampling.9000.dataset/",
    "model_path": "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/Geneformer-V2-104M",
    "output_base_dir": "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.stratified.sampling.9000.dataset.out",

    # Data splitting
    "train_ratio": 0.7,
    "eval_ratio": 0.15,
    "test_ratio": 0.15,
    "random_state": 42,

    # 🎯 OPTIMIZED TRAINING PARAMETERS FOR HIGHER ACCURACY
    "num_train_epochs": 8.0,            # ✅ Increased from 0.9 to 8.0 for better learning
    "learning_rate": 0.0001,            # ✅ Reduced from 0.0008 for stability
    "lr_scheduler_type": "cosine",      # ✅ Changed from polynomial to cosine (better for longer training)
    "warmup_ratio": 0.1,                # ✅ Use ratio instead of absolute steps
    "weight_decay": 0.01,               # ✅ Reduced from 0.03 to 0.01 (less aggressive regularization)
    "per_device_train_batch_size": 2,   # ✅ Reduced from 4 to 2 for stability
    "per_device_eval_batch_size": 2,    # ✅ Match training batch size
    "forward_batch_size": 2,            # ✅ Match batch size
    "freeze_layers": 1,                 # ✅ Reduced from 2 to 1 (less restrictive)
    "nproc": 16,
    "ngpu": 1,

    # Cell type filtering (None = use all cell types)
    "filter_cell_types": None, 

    # Plotting
    "plot_dpi": 300,
    "plot_figsize": (10, 8),
    "color_palette": "Set3"
}


# In[6]:


import gc 
import datetime 

print("🚀 Geneformer Hypoxia Condition Classification - OPTIMIZED VERSION")
print("="*80)
print("🎯 OPTIMIZED FOR HIGHER ACCURACY:")
print(f"  • Training epochs: {CONFIG['num_train_epochs']} (was 0.9)")
print(f"  • Learning rate: {CONFIG['learning_rate']} (was 0.0008)")
print(f"  • Scheduler: {CONFIG['lr_scheduler_type']} (was polynomial)")
print(f"  • Batch size: {CONFIG['per_device_train_batch_size']} (was 4)")
print(f"  • Freeze layers: {CONFIG['freeze_layers']} (was 2)")
print("="*80)

# GPU memory management
print(f"\n🔧 GPU Memory Management:")
torch.cuda.empty_cache()
gc.collect()
print(f"GPU memory allocated: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"GPU memory reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")

# Setup output directory
current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
output_dir = f"{CONFIG['output_base_dir']}/{datestamp}"
output_prefix = "hypoxia_condition_classifier_optimized"  # Updated for optimized condition classification


# In[7]:


os.makedirs(output_dir, exist_ok=True)
print(f"📁 Output: {output_dir}")


# In[8]:


from pathlib import Path
from collections import Counter
from datasets import load_from_disk, Dataset, DatasetDict, concatenate_datasets

print("\n📊 STEP 1: Loading and exploring dataset...")

# --- validate path ---
ds_path = Path(CONFIG["dataset_path"])
if not ds_path.exists():
    raise FileNotFoundError(f"❌ Dataset not found: {ds_path}\nPlease check CONFIG['dataset_path'].")

# --- load ---
try:
    raw = load_from_disk(str(ds_path))
except Exception as e:
    raise RuntimeError(f"❌ Failed to load dataset: {e}")

# --- normalize to a single Dataset ---
if isinstance(raw, DatasetDict):
    dataset = concatenate_datasets(list(raw.values()))
elif isinstance(raw, Dataset):
    dataset = raw
else:
    raise TypeError(f"Unexpected dataset type: {type(raw)}")

n = len(dataset)
print(f"✅ Loaded {n:,} cells")

# --- schema/columns ---
print(f"\n📋 Dataset structure:")
print("  Columns:", dataset.column_names)
try:
    print("  Features:", dataset.features)
except Exception:
    pass  # features may not be present in some edge cases

# --- required columns ---
required = ["condition", "input_ids"]  # 
missing = [c for c in required if c not in dataset.column_names]
if missing:
    raise KeyError(f"❌ Missing required columns: {missing}")

# --- condition distribution ---
condition_counts = Counter(dataset["condition"])  
print(f"\n🏷️  Conditions ({len(condition_counts)} total):")  
for cond_name, count in condition_counts.most_common():
    pct = 100 * count / n
    print(f"  {cond_name}: {count} ({pct:.1f}%)")

# --- sample row ---
if n > 0:
    sample = dataset[0]
    print(
        f"\n🔎 Sample cell → condition: {sample.get('condition')} | " 
        f"tokens: {len(sample.get('input_ids', []))}"
    )


# In[9]:


# =============================================================================
# STEP 2: CREATE DATA SPLITS
# =============================================================================

print(f"\n📊 STEP 2: Creating data splits...")

def create_simple_splits(dataset, train_ratio=0.7, eval_ratio=0.15, test_ratio=0.15, random_state=42):
    """
    Create simple train/eval/test splits using the datasets library
    This is the user's preferred simple split method
    """
    print(f"Splits: {train_ratio*100:.0f}% train, {eval_ratio*100:.0f}% eval, {test_ratio*100:.0f}% test")

    # Shuffle and split
    dataset = dataset.shuffle(seed=random_state)
    n_total = len(dataset)
    n_train = int(n_total * train_ratio)
    n_eval = int(n_total * eval_ratio)
    n_test = n_total - n_train - n_eval

    # Create split labels
    split_labels = ["train"] * n_train + ["eval"] * n_eval + ["test"] * n_test
    dataset = dataset.map(lambda example, idx: {"split": split_labels[idx]}, with_indices=True)

    return dataset, {
        'train_size': n_train,
        'eval_size': n_eval,
        'test_size': n_test,
        'total_size': n_total
    }

# Create splits
dataset, split_info = create_simple_splits(
    dataset,
    CONFIG["train_ratio"],
    CONFIG["eval_ratio"],
    CONFIG["test_ratio"],
    CONFIG["random_state"]
)

# Show split distribution
split_counts = Counter([example["split"] for example in dataset])
print(f"\n📊 Split sizes:")
for split, count in split_counts.items():
    print(f"  {split}: {count} cells")

# Save dataset with splits
complete_dataset_path = f"{output_dir}/hypoxia_condition_with_splits_ready.dataset"
dataset.save_to_disk(complete_dataset_path)
print(f"✅ Saved dataset with splits")


# In[10]:


# =============================================================================
# STEP 3: SETUP CLASSIFIER WITH OPTIMIZED TRAINING ARGUMENTS
# =============================================================================

print(f"\n🔧 STEP 3: Setting up classifier with OPTIMIZED training arguments...")

# 🚀 ENHANCED TRAINING ARGUMENTS FOR HIGHER ACCURACY
training_args = {
    "num_train_epochs": CONFIG["num_train_epochs"],
    "learning_rate": CONFIG["learning_rate"],
    "lr_scheduler_type": CONFIG["lr_scheduler_type"],
    "warmup_ratio": CONFIG["warmup_ratio"],  # Use ratio instead of absolute steps
    "weight_decay": CONFIG["weight_decay"],
    "per_device_train_batch_size": CONFIG["per_device_train_batch_size"],
    "per_device_eval_batch_size": CONFIG["per_device_eval_batch_size"],
    "seed": CONFIG["random_state"],
    
    # ✅ NEW: Enhanced training for better accuracy
    "gradient_accumulation_steps": 4,    # Effective batch size = 2 * 4 = 8
    "max_grad_norm": 1.0,               # Gradient clipping for stability
    "save_steps": 100,                   # Save checkpoints during training
    "eval_steps": 100,                   # Evaluate during training
    "logging_steps": 50,                 # Log progress frequently
    "save_total_limit": 3,               # Keep best 3 checkpoints
    "load_best_model_at_end": True,     # Load best model at end
    "metric_for_best_model": "accuracy", # Use accuracy to select best model
    
    # ✅ NEW: Better regularization and stability
    "label_smoothing_factor": 0.1,      # Prevent overconfidence
    "dataloader_pin_memory": False,     # Save memory
    "remove_unused_columns": True,      # Clean up data
    
    # Memory optimization
    "gradient_checkpointing": True,      # Save memory during training
    "fp16": False,                       # Keep as False for stability
}

print("🎯 OPTIMIZED TRAINING ARGUMENTS:")
for key, value in training_args.items():
    print(f"  {key}: {value}")

# Initialize classifier with optimized settings
cc = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "condition", "states": "all"},  
    filter_data=CONFIG["filter_cell_types"],
    training_args=training_args,
    max_ncells=None,
    freeze_layers=CONFIG["freeze_layers"],
    num_crossval_splits=1,
    forward_batch_size=CONFIG["forward_batch_size"],
    nproc=CONFIG["nproc"],
    ngpu=CONFIG["ngpu"]
)

print("✅ OPTIMIZED Classifier ready")


# In[11]:


# =============================================================================
# STEP 4: PREPARE DATA FOR TRAINING
# =============================================================================

print(f"\n📋 STEP 4: Preparing data for training...")

# Add debugging code to check dataset structure
print(f"\n🔍 DEBUGGING DATASET STRUCTURE:")
print(f"Dataset type: {type(dataset)}")
print(f"Dataset length: {len(dataset)}")
print(f"Columns: {dataset.column_names}")
print(f"Features: {dataset.features}")

# Check first few rows
print(f"\n�� SAMPLE DATA:")
for i in range(min(3, len(dataset))):
    print(f"Row {i}: {dataset[i]}")

# Check split distribution
print(f"\n�� SPLIT DISTRIBUTION:")
split_counts = Counter([example["split"] for example in dataset])
for split, count in split_counts.items():
    print(f"  {split}: {count} cells")

# Split dictionaries for Geneformer
train_test_split_dict = {"attr_key": "split", "train": ["train", "eval"], "test": ["test"]}
train_eval_split_dict = {"attr_key": "split", "train": ["train"], "eval": ["eval"]}

# =============================================================================
# CALCULATE SPLIT PERCENTAGES
# =============================================================================

print(f"\n📊 SPLIT PERCENTAGE ANALYSIS:")
print("="*60)

total_cells = len(dataset)
print(f"Total cells: {total_cells:,}")

# Calculate percentages for each split
print(f"\n📈 INDIVIDUAL SPLIT PERCENTAGES:")
for split, count in split_counts.items():
    percentage = (count / total_cells) * 100
    print(f"  {split.upper()}: {count:,} cells ({percentage:.2f}%)")

# Calculate combined percentages (train+eval vs test)
train_eval_combined = split_counts.get('train', 0) + split_counts.get('eval', 0)
test_count = split_counts.get('test', 0)

print(f"\n�� COMBINED SPLIT PERCENTAGES:")
print(f"  TRAIN + EVAL: {train_eval_combined:,} cells ({(train_eval_combined/total_cells)*100:.2f}%)")
print(f"  TEST: {test_count:,} cells ({(test_count/total_cells)*100:.2f}%)")

# Calculate training vs evaluation percentages
train_count = split_counts.get('train', 0)
eval_count = split_counts.get('eval', 0)

print(f"\n📈 TRAINING vs EVALUATION PERCENTAGES:")
print(f"  TRAIN: {train_count:,} cells ({(train_count/total_cells)*100:.2f}%)")
print(f"  EVAL: {eval_count:,} cells ({(eval_count/total_cells)*100:.2f}%)")

# Calculate effective training percentage (train + eval for actual training)
effective_training = train_count + eval_count
print(f"\n�� EFFECTIVE TRAINING PERCENTAGE:")
print(f"  TRAIN + EVAL (for training): {effective_training:,} cells ({(effective_training/total_cells)*100:.2f}%)")
print(f"  TEST (held out): {test_count:,} cells ({(test_count/total_cells)*100:.2f}%)")

# Verify the percentages add up to 100%
total_percentage = sum([(count/total_cells)*100 for count in split_counts.values()])
print(f"\n✅ VERIFICATION: Total percentage = {total_percentage:.2f}%")

# Save split analysis to file
split_analysis_file = f"{output_dir}/{output_prefix}_split_analysis.txt"
with open(split_analysis_file, 'w') as f:
    f.write("DATASET SPLIT ANALYSIS\n")
    f.write("="*50 + "\n\n")
    f.write(f"Total cells: {total_cells:,}\n\n")
    
    f.write("INDIVIDUAL SPLIT PERCENTAGES:\n")
    f.write("-" * 30 + "\n")
    for split, count in split_counts.items():
        percentage = (count / total_cells) * 100
        f.write(f"{split.upper()}: {count:,} cells ({percentage:.2f}%)\n")
    
    f.write(f"\nCOMBINED SPLIT PERCENTAGES:\n")
    f.write("-" * 30 + "\n")
    f.write(f"TRAIN + EVAL: {train_eval_combined:,} cells ({(train_eval_combined/total_cells)*100:.2f}%)\n")
    f.write(f"TEST: {test_count:,} cells ({(test_count/total_cells)*100:.2f}%)\n")
    
    f.write(f"\nTRAINING vs EVALUATION PERCENTAGES:\n")
    f.write("-" * 30 + "\n")
    f.write(f"TRAIN: {train_count:,} cells ({(train_count/total_cells)*100:.2f}%)\n")
    f.write(f"EVAL: {eval_count:,} cells ({(eval_count/total_cells)*100:.2f}%)\n")
    
    f.write(f"\nEFFECTIVE TRAINING PERCENTAGE:\n")
    f.write("-" * 30 + "\n")
    f.write(f"TRAIN + EVAL (for training): {effective_training:,} cells ({(effective_training/total_cells)*100:.2f}%)\n")
    f.write(f"TEST (held out): {test_count:,} cells ({(test_count/total_cells)*100:.2f}%)\n")
    
    f.write(f"\nVERIFICATION: Total percentage = {total_percentage:.2f}%\n")

print(f"📝 Split analysis saved to: {split_analysis_file}")


# In[12]:


# =============================================================================
# CELL TYPE DISTRIBUTION PER SPLIT
# =============================================================================

print(f"\n🔬 CELL TYPE DISTRIBUTION PER SPLIT:")
print("="*80)

# Get all unique cell types from the dataset
all_cell_types = set()
for example in dataset:
    if 'class' in example:
        all_cell_types.add(example['class'])
    elif 'cell_type' in example:
        all_cell_types.add(example['cell_type'])

all_cell_types = sorted(list(all_cell_types))
print(f"Found {len(all_cell_types)} cell types: {all_cell_types}")

# Analyze each split
for split_name in ['train', 'eval', 'test']:
    print(f"\n📊 {split_name.upper()} SPLIT ANALYSIS:")
    print("-" * 50)
    
    # Filter data for this split
    split_data = [example for example in dataset if example['split'] == split_name]
    split_total = len(split_data)
    
    if split_total == 0:
        print(f"  No data found for {split_name} split")
        continue
    
    print(f"  Total cells: {split_total:,}")
    
    # Count cell types in this split
    split_cell_counts = Counter()
    for example in split_data:
        if 'class' in example:
            split_cell_counts[example['class']] += 1
        elif 'cell_type' in example:
            split_cell_counts[example['cell_type']] += 1
    
    # Calculate and display percentages
    print(f"  Cell type distribution:")
    for cell_type in all_cell_types:
        count = split_cell_counts.get(cell_type, 0)
        percentage = (count / split_total) * 100
        print(f"    {cell_type}: {count:,} cells ({percentage:.2f}%)")
    
    # Show top 5 cell types
    top_cell_types = split_cell_counts.most_common(5)
    print(f"  Top 5 cell types:")
    for i, (cell_type, count) in enumerate(top_cell_types, 1):
        percentage = (count / split_total) * 100
        print(f"    {i}. {cell_type}: {count:,} cells ({percentage:.2f}%)")

# =============================================================================
# COMPARATIVE ANALYSIS
# =============================================================================

print(f"\n📈 COMPARATIVE ANALYSIS:")
print("="*80)

# Create a comprehensive comparison table
comparison_data = []
for split_name in ['train', 'eval', 'test']:
    split_data = [example for example in dataset if example['split'] == split_name]
    split_total = len(split_data)
    
    if split_total == 0:
        continue
    
    split_cell_counts = Counter()
    for example in split_data:
        if 'class' in example:
            split_cell_counts[example['class']] += 1
        elif 'cell_type' in example:
            split_cell_counts[example['cell_type']] += 1
    
    for cell_type in all_cell_types:
        count = split_cell_counts.get(cell_type, 0)
        percentage = (count / split_total) * 100
        comparison_data.append({
            'split': split_name,
            'cell_type': cell_type,
            'count': count,
            'percentage': percentage
        })

# Create DataFrame for better visualization
import pandas as pd
df_comparison = pd.DataFrame(comparison_data)

# Pivot table for better comparison
pivot_table = df_comparison.pivot_table(
    index='cell_type', 
    columns='split', 
    values='percentage', 
    fill_value=0
)

print(f"\n📊 CELL TYPE PERCENTAGES BY SPLIT:")
print(pivot_table.round(2))

# =============================================================================
# SAVE DETAILED ANALYSIS
# =============================================================================

# Save detailed cell type analysis to file
cell_type_analysis_file = f"{output_dir}/{output_prefix}_cell_type_analysis.txt"
with open(cell_type_analysis_file, 'w') as f:
    f.write("CELL TYPE DISTRIBUTION PER SPLIT\n")
    f.write("="*60 + "\n\n")
    
    f.write(f"Total cell types found: {len(all_cell_types)}\n")
    f.write(f"Cell types: {', '.join(all_cell_types)}\n\n")
    
    for split_name in ['train', 'eval', 'test']:
        split_data = [example for example in dataset if example['split'] == split_name]
        split_total = len(split_data)
        
        if split_total == 0:
            f.write(f"\n{split_name.upper()} SPLIT: No data found\n")
            continue
        
        f.write(f"\n{split_name.upper()} SPLIT:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Total cells: {split_total:,}\n\n")
        
        split_cell_counts = Counter()
        for example in split_data:
            if 'class' in example:
                split_cell_counts[example['class']] += 1
            elif 'cell_type' in example:
                split_cell_counts[example['cell_type']] += 1
        
        f.write("Cell type distribution:\n")
        for cell_type in all_cell_types:
            count = split_cell_counts.get(cell_type, 0)
            percentage = (count / split_total) * 100
            f.write(f"  {cell_type}: {count:,} cells ({percentage:.2f}%)\n")
        
        # Top 5 cell types
        top_cell_types = split_cell_counts.most_common(5)
        f.write(f"\nTop 5 cell types:\n")
        for i, (cell_type, count) in enumerate(top_cell_types, 1):
            percentage = (count / split_total) * 100
            f.write(f"  {i}. {cell_type}: {count:,} cells ({percentage:.2f}%)\n")
    
    # Add comparison table
    f.write(f"\n\nCOMPARATIVE TABLE:\n")
    f.write("="*60 + "\n")
    f.write(pivot_table.round(2).to_string())

print(f"�� Cell type analysis saved to: {cell_type_analysis_file}")

# =============================================================================
# BALANCE CHECK
# =============================================================================

print(f"\n⚖️  BALANCE CHECK:")
print("="*60)

# Check if cell type distribution is similar across splits
for cell_type in all_cell_types:
    percentages = []
    for split_name in ['train', 'eval', 'test']:
        split_data = [example for example in dataset if example['split'] == split_name]
        split_total = len(split_data)
        
        if split_total == 0:
            continue
        
        split_cell_counts = Counter()
        for example in split_data:
            if 'class' in example:
                split_cell_counts[example['class']] += 1
            elif 'cell_type' in example:
                split_cell_counts[example['cell_type']] += 1
        
        count = split_cell_counts.get(cell_type, 0)
        percentage = (count / split_total) * 100
        percentages.append(percentage)
    
    if len(percentages) >= 2:
        max_pct = max(percentages)
        min_pct = min(percentages)
        difference = max_pct - min_pct
        
        if difference < 1.0:
            status = "✅ Balanced"
        elif difference < 5.0:
            status = "⚠️  Fairly balanced"
        else:
            status = "❌ Imbalanced"
        
        print(f"  {cell_type}: {min_pct:.2f}% - {max_pct:.2f}% (diff: {difference:.2f}%) {status}")

print(f"\n✅ Cell type distribution analysis complete!")


# In[13]:


# =============================================================================
# CELL TYPE DISTRIBUTION PER CONDITION
# =============================================================================

print(f"\n🔬 CELL TYPE DISTRIBUTION PER CONDITION:")
print("="*80)

# Get all unique conditions from the dataset
all_conditions = set()
for example in dataset:
    if 'condition' in example:
        all_conditions.add(example['condition'])

all_conditions = sorted(list(all_conditions))
print(f"Found {len(all_conditions)} conditions: {all_conditions}")

# Analyze each condition
for condition_name in all_conditions:
    print(f"\n📊 {condition_name.upper()} CONDITION ANALYSIS:")
    print("-" * 50)
    
    # Filter data for this condition
    condition_data = [example for example in dataset if example['condition'] == condition_name]
    condition_total = len(condition_data)
    
    if condition_total == 0:
        print(f"  No data found for {condition_name} condition")
        continue
    
    print(f"  Total cells: {condition_total:,}")
    
    # Count cell types in this condition
    condition_cell_counts = Counter()
    for example in condition_data:
        if 'class' in example:
            condition_cell_counts[example['class']] += 1
        elif 'cell_type' in example:
            condition_cell_counts[example['cell_type']] += 1
    
    # Calculate and display percentages
    print(f"  Cell type distribution:")
    for cell_type in all_cell_types:
        count = condition_cell_counts.get(cell_type, 0)
        percentage = (count / condition_total) * 100
        print(f"    {cell_type}: {count:,} cells ({percentage:.2f}%)")
    
    # Show top 5 cell types
    top_cell_types = condition_cell_counts.most_common(5)
    print(f"  Top 5 cell types:")
    for i, (cell_type, count) in enumerate(top_cell_types, 1):
        percentage = (count / condition_total) * 100
        print(f"    {i}. {cell_type}: {count:,} cells ({percentage:.2f}%)")

# =============================================================================
# COMBINED SPLIT + CONDITION ANALYSIS
# =============================================================================

print(f"\n🔬 COMBINED SPLIT + CONDITION ANALYSIS:")
print("="*80)

# Create comprehensive analysis for each split-condition combination
for split_name in ['train', 'eval', 'test']:
    print(f"\n📊 {split_name.upper()} SPLIT - BY CONDITION:")
    print("-" * 50)
    
    split_data = [example for example in dataset if example['split'] == split_name]
    split_total = len(split_data)
    
    if split_total == 0:
        print(f"  No data found for {split_name} split")
        continue
    
    print(f"  Total cells in {split_name}: {split_total:,}")
    
    # Analyze each condition within this split
    for condition_name in all_conditions:
        condition_in_split = [example for example in split_data if example['condition'] == condition_name]
        condition_count = len(condition_in_split)
        
        if condition_count == 0:
            continue
        
        percentage_in_split = (condition_count / split_total) * 100
        print(f"    {condition_name}: {condition_count:,} cells ({percentage_in_split:.2f}% of {split_name})")
        
        # Cell type distribution within this condition in this split
        condition_cell_counts = Counter()
        for example in condition_in_split:
            if 'class' in example:
                condition_cell_counts[example['class']] += 1
            elif 'cell_type' in example:
                condition_cell_counts[example['cell_type']] += 1
        
        # Show top 3 cell types for this condition in this split
        top_cell_types = condition_cell_counts.most_common(3)
        if top_cell_types:
            print(f"      Top cell types: {', '.join([f'{ct}({count})' for ct, count in top_cell_types])}")

# =============================================================================
# CONDITION BALANCE CHECK
# =============================================================================

print(f"\n⚖️  CONDITION BALANCE CHECK:")
print("="*60)

# Check if condition distribution is similar across splits
for condition_name in all_conditions:
    percentages = []
    for split_name in ['train', 'eval', 'test']:
        split_data = [example for example in dataset if example['split'] == split_name]
        split_total = len(split_data)
        
        if split_total == 0:
            continue
        
        condition_in_split = [example for example in split_data if example['condition'] == condition_name]
        condition_count = len(condition_in_split)
        percentage = (condition_count / split_total) * 100
        percentages.append(percentage)
    
    if len(percentages) >= 2:
        max_pct = max(percentages)
        min_pct = min(percentages)
        difference = max_pct - min_pct
        
        if difference < 1.0:
            status = "✅ Balanced"
        elif difference < 5.0:
            status = "⚠️  Fairly balanced"
        else:
            status = "❌ Imbalanced"
        
        print(f"  {condition_name}: {min_pct:.2f}% - {max_pct:.2f}% (diff: {difference:.2f}%) {status}")

# =============================================================================
# DETAILED CROSS-TABULATION
# =============================================================================

print(f"\n�� DETAILED CROSS-TABULATION:")
print("="*80)

# Create cross-tabulation of split vs condition
cross_tab_data = []
for example in dataset:
    cross_tab_data.append({
        'split': example['split'],
        'condition': example['condition'],
        'cell_type': example.get('class', example.get('cell_type', 'unknown'))
    })

df_cross = pd.DataFrame(cross_tab_data)

# Split vs Condition counts
print(f"\n�� SPLIT vs CONDITION COUNTS:")
split_condition_counts = pd.crosstab(df_cross['split'], df_cross['condition'])
print(split_condition_counts)

# Split vs Condition percentages
print(f"\n�� SPLIT vs CONDITION PERCENTAGES (row-wise):")
split_condition_pct = pd.crosstab(df_cross['split'], df_cross['condition'], normalize='index') * 100
print(split_condition_pct.round(2))

# Cell type vs Condition counts
print(f"\n📈 CELL TYPE vs CONDITION COUNTS:")
celltype_condition_counts = pd.crosstab(df_cross['cell_type'], df_cross['condition'])
print(celltype_condition_counts)

# Cell type vs Condition percentages
print(f"\n📈 CELL TYPE vs CONDITION PERCENTAGES (row-wise):")
celltype_condition_pct = pd.crosstab(df_cross['cell_type'], df_cross['condition'], normalize='index') * 100
print(celltype_condition_pct.round(2))

# =============================================================================
# SAVE COMPREHENSIVE ANALYSIS
# =============================================================================

# Save comprehensive analysis to file
comprehensive_analysis_file = f"{output_dir}/{output_prefix}_comprehensive_analysis.txt"
with open(comprehensive_analysis_file, 'w') as f:
    f.write("COMPREHENSIVE CELL TYPE AND CONDITION ANALYSIS\n")
    f.write("="*70 + "\n\n")
    
    f.write(f"Total cell types found: {len(all_cell_types)}\n")
    f.write(f"Cell types: {', '.join(all_cell_types)}\n\n")
    
    f.write(f"Total conditions found: {len(all_conditions)}\n")
    f.write(f"Conditions: {', '.join(all_conditions)}\n\n")
    
    # Condition analysis
    f.write("CONDITION ANALYSIS:\n")
    f.write("="*30 + "\n")
    for condition_name in all_conditions:
        condition_data = [example for example in dataset if example['condition'] == condition_name]
        condition_total = len(condition_data)
        
        if condition_total == 0:
            f.write(f"\n{condition_name.upper()}: No data found\n")
            continue
        
        f.write(f"\n{condition_name.upper()}:\n")
        f.write(f"Total cells: {condition_total:,}\n")
        
        condition_cell_counts = Counter()
        for example in condition_data:
            if 'class' in example:
                condition_cell_counts[example['class']] += 1
            elif 'cell_type' in example:
                condition_cell_counts[example['cell_type']] += 1
        
        f.write("Cell type distribution:\n")
        for cell_type in all_cell_types:
            count = condition_cell_counts.get(cell_type, 0)
            percentage = (count / condition_total) * 100
            f.write(f"  {cell_type}: {count:,} cells ({percentage:.2f}%)\n")
    
    # Split analysis
    f.write(f"\n\nSPLIT ANALYSIS:\n")
    f.write("="*30 + "\n")
    for split_name in ['train', 'eval', 'test']:
        split_data = [example for example in dataset if example['split'] == split_name]
        split_total = len(split_data)
        
        if split_total == 0:
            f.write(f"\n{split_name.upper()}: No data found\n")
            continue
        
        f.write(f"\n{split_name.upper()}:\n")
        f.write(f"Total cells: {split_total:,}\n")
        
        # Condition distribution in this split
        f.write("Condition distribution:\n")
        for condition_name in all_conditions:
            condition_in_split = [example for example in split_data if example['condition'] == condition_name]
            condition_count = len(condition_in_split)
            percentage = (condition_count / split_total) * 100
            f.write(f"  {condition_name}: {condition_count:,} cells ({percentage:.2f}%)\n")
    
    # Cross-tabulations
    f.write(f"\n\nCROSS-TABULATIONS:\n")
    f.write("="*30 + "\n")
    f.write("Split vs Condition counts:\n")
    f.write(split_condition_counts.to_string())
    f.write("\n\nSplit vs Condition percentages:\n")
    f.write(split_condition_pct.round(2).to_string())
    f.write("\n\nCell type vs Condition counts:\n")
    f.write(celltype_condition_counts.to_string())
    f.write("\n\nCell type vs Condition percentages:\n")
    f.write(celltype_condition_pct.round(2).to_string())

print(f"📝 Comprehensive analysis saved to: {comprehensive_analysis_file}")

print(f"\n✅ Comprehensive cell type and condition analysis complete!")


# In[ ]:





# In[14]:


try:
    print(f"\n🔍 PREPARING DATA:")
    print(f"Input file: {complete_dataset_path}")
    print(f"Output dir: {output_dir}")
    print(f"Output prefix: {output_prefix}")
    print(f"Split dict: {train_test_split_dict}")
    
    cc.prepare_data(
        input_data_file=complete_dataset_path,
        output_directory=output_dir,
        output_prefix=output_prefix,
        split_id_dict=train_test_split_dict
    )
    print("✅ Data prepared")
except Exception as e:
    print(f"❌ Data preparation failed: {e}")
    print(f"Exception type: {type(e)}")
    import traceback
    traceback.print_exc()
    # exit()


# In[15]:


# Add this before training:
import torch
import gc

# Clear everything
torch.cuda.empty_cache()
gc.collect()

# Set memory fraction
torch.cuda.set_per_process_memory_fraction(0.8)  # Use only 80% of GPU memory

# Check memory before training
print(f"GPU memory before training: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")


# In[14]:


# Add this right before your training code:
import torch
import gc

print("🧹 Clearing GPU memory...")
torch.cuda.empty_cache()
gc.collect()

# Force garbage collection multiple times
for _ in range(3):
    gc.collect()
    torch.cuda.empty_cache()

# Set memory fraction to be conservative
torch.cuda.set_per_process_memory_fraction(0.7)  # Use only 70% of GPU

print(f"GPU memory after cleanup: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"GPU memory reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")



# Add this debugging code:
print("🔍 Memory usage breakdown:")
print(f"PyTorch allocated: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"PyTorch reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")
print(f"Total GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.2f} GB")


# In[ ]:


# =============================================================================
# STEP 5: TRAIN MODEL WITH OPTIMIZED SETTINGS
# =============================================================================

print(f"\n🚀 STEP 5: Training model with OPTIMIZED settings for higher accuracy...")

# Enhanced memory management before training
print("🧹 Enhanced GPU memory cleanup...")
torch.cuda.empty_cache()
gc.collect()

# Force garbage collection multiple times
for _ in range(3):
    gc.collect()
    torch.cuda.empty_cache()

# Set memory fraction to be conservative
torch.cuda.set_per_process_memory_fraction(0.7)  # Use only 70% of GPU

print(f"GPU memory after cleanup: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"GPU memory reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")

# Check required files exist
prepared_dataset_file = f"{output_dir}/{output_prefix}_labeled_train.dataset"
id_class_dict_file = f"{output_dir}/{output_prefix}_id_class_dict.pkl"

if not os.path.exists(prepared_dataset_file):
    print(f"❌ Prepared dataset not found: {prepared_dataset_file}")
    exit()

if not os.path.exists(id_class_dict_file):
    print(f"❌ ID class dictionary not found: {id_class_dict_file}")
    exit()

# Train model with enhanced monitoring and timing
try:
    print(f"�� Starting OPTIMIZED training...")
    print(f"  Model directory: {CONFIG['model_path']}")
    print(f"  Prepared data: {prepared_dataset_file}")
    print(f"  Output: {output_dir}")
    print(f"  Epochs: {CONFIG['num_train_epochs']}")
    print(f"  Learning rate: {CONFIG['learning_rate']}")
    print(f"  Batch size: {CONFIG['per_device_train_batch_size']}")
    print(f"  Effective batch size: {CONFIG['per_device_train_batch_size'] * training_args['gradient_accumulation_steps']}")
    
    # ⏱️ START TIMING
    import time
    training_start_time = time.time()
    training_start_datetime = datetime.datetime.now()
    print(f"⏰ Training started at: {training_start_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
    
    training_metrics = cc.validate(
        model_directory=CONFIG["model_path"],
        prepared_input_data_file=prepared_dataset_file,
        id_class_dict_file=id_class_dict_file,
        output_directory=output_dir,
        output_prefix=output_prefix,
        split_id_dict=train_eval_split_dict
    )
    
    # ⏱️ END TIMING
    training_end_time = time.time()
    training_end_datetime = datetime.datetime.now()
    training_duration = training_end_time - training_start_time
    
    print("✅ OPTIMIZED Training completed!")
    print(f"⏰ Training ended at: {training_end_datetime.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"⏱️  Training duration: {training_duration:.2f} seconds ({training_duration/60:.2f} minutes)")
    
    if training_metrics:
        acc = training_metrics.get('accuracy')
        f1 = training_metrics.get('macro_f1')
    
        if acc is not None and isinstance(acc, (int, float)):
              print(f"�� Training accuracy: {acc:.4f}")
        else:
              print(f"�� Training accuracy: {acc}")
        
        if f1 is not None and isinstance(f1, (int, float)):
              print(f"�� Training F1: {f1:.4f}")
        else:
              print(f"🎯 Training F1: {f1}")
    else:
           print("�� Training completed but no metrics returned")

    # 📝 SAVE TIMING INFORMATION TO FILE
    timing_file = f"{output_dir}/{output_prefix}_training_timing.txt"
    with open(timing_file, 'w') as f:
        f.write("GENEFORMER TRAINING TIMING REPORT\n")
        f.write("="*50 + "\n\n")
        f.write(f"Training Start: {training_start_datetime.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Training End: {training_end_datetime.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Duration: {training_duration:.2f} seconds\n")
        f.write(f"Duration: {training_duration/60:.2f} minutes\n")
        f.write(f"Duration: {training_duration/3600:.2f} hours\n\n")
        
        f.write("TRAINING CONFIGURATION:\n")
        f.write("-" * 30 + "\n")
        f.write(f"Epochs: {CONFIG['num_train_epochs']}\n")
        f.write(f"Learning Rate: {CONFIG['learning_rate']}\n")
        f.write(f"Batch Size: {CONFIG['per_device_train_batch_size']}\n")
        f.write(f"Effective Batch Size: {CONFIG['per_device_train_batch_size'] * training_args['gradient_accumulation_steps']}\n")
        f.write(f"Freeze Layers: {CONFIG['freeze_layers']}\n")
        f.write(f"Scheduler: {CONFIG['lr_scheduler_type']}\n\n")
        
        f.write("PERFORMANCE METRICS:\n")
        f.write("-" * 30 + "\n")
        if training_metrics:
            if 'accuracy' in training_metrics:
                f.write(f"Accuracy: {training_metrics['accuracy']:.4f}\n")
            if 'macro_f1' in training_metrics:
                f.write(f"Macro F1-Score: {training_metrics['macro_f1']:.4f}\n")
        
        f.write(f"\nTraining completed successfully!\n")
        f.write(f"Output directory: {output_dir}\n")
    
    print(f"�� Timing information saved to: {timing_file}")

except Exception as e:
    # ⏱️ TIMING FOR FAILED TRAINING
    training_end_time = time.time()
    training_duration = training_end_time - training_start_time
    
    print(f"❌ Training failed: {e}")
    print(f"⏰ Training failed at: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"⏱️  Training duration before failure: {training_duration:.2f} seconds ({training_duration/60:.2f} minutes)")
    
    # �� SAVE FAILED TIMING INFORMATION
    timing_file = f"{output_dir}/{output_prefix}_training_timing_FAILED.txt"
    with open(timing_file, 'w') as f:
        f.write("GENEFORMER TRAINING TIMING REPORT (FAILED)\n")
        f.write("="*50 + "\n\n")
        f.write(f"Training Start: {training_start_datetime.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Training Failed: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Duration before failure: {training_duration:.2f} seconds\n")
        f.write(f"Duration before failure: {training_duration/60:.2f} minutes\n\n")
        f.write(f"Error: {e}\n")
        f.write(f"Error type: {type(e)}\n")
    
    print(f"📝 Failed timing information saved to: {timing_file}")
    print(f"Exception type: {type(e)}")
    import traceback
    traceback.print_exc()
    # exit()


# In[ ]:


# =============================================================================
# STEP 6: EVALUATE MODEL
# =============================================================================

print(f"\n🧪 STEP 6: Evaluating model...")

# Find trained model
model_dirs = [d for d in os.listdir(output_dir) if "geneformer_cellClassifier" in d and output_prefix in d]
if model_dirs:
    trained_model_dir = os.path.join(output_dir, model_dirs[0], "ksplit1")
    print(f"✅ Found model: {trained_model_dir}")
else:
    print("❌ Model not found")
    exit()

# Create evaluation classifier
eval_cc = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "condition", "states": "all"},  
    forward_batch_size=CONFIG["forward_batch_size"],
    nproc=CONFIG["nproc"]
)

# Test files
test_data_file = f"{output_dir}/{output_prefix}_labeled_test.dataset"

if not os.path.exists(test_data_file):
    print(f"❌ Test dataset not found: {test_data_file}")
    exit()


# In[14]:


# Evaluate model
try:
    all_metrics_test = eval_cc.evaluate_saved_model(
        model_directory=trained_model_dir,
        id_class_dict_file=id_class_dict_file,
        test_data_file=test_data_file,
        output_directory=output_dir,
        output_prefix=f"{output_prefix}_test"
    )
    print("✅ Evaluation completed!")

    # Show results
    print(f"\n📊 TEST RESULTS:")
    for key, value in all_metrics_test.items():
        if isinstance(value, (int, float)):
            print(f"  {key}: {value:.4f}")

except Exception as e:
    print(f"❌ Evaluation failed: {e}")
    all_metrics_test = None



# In[ ]:


# =============================================================================
# STEP 7: CREATE OPTIMIZED VISUALIZATIONS WITH ROTATED LABELS
# =============================================================================

if all_metrics_test:
    print(f"\n📊 STEP 7: Creating OPTIMIZED visualizations with rotated labels...")

    import matplotlib.pyplot as plt
    import seaborn as sns

    # Get unique conditions (not cell types)  # ✅ FIXED: Changed from cell types to conditions
    unique_conditions = list(condition_counts.keys())  # ✅ FIXED: Changed from class_counts to condition_counts

    # Custom confusion matrix plot with rotated labels
    try:
        conf_matrix = all_metrics_test["conf_matrix"]
        
        plt.figure(figsize=(12, 10))
        sns.heatmap(
            conf_matrix, 
            annot=True, 
            fmt='d', 
            cmap='Blues',
            xticklabels=unique_conditions,
            yticklabels=unique_conditions
        )
        
        # Rotate x-axis labels for readability
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        
        plt.title(f"Confusion Matrix - Hypoxia Conditions (OPTIMIZED)", fontsize=16, pad=20)
        plt.xlabel("Predicted", fontsize=14)
        plt.ylabel("Actual", fontsize=14)
        
        # Adjust layout to prevent label cutoff
        plt.tight_layout()
        
        # Save the plot
        conf_mat_file = f"{output_dir}/{output_prefix}_custom_conf_mat_optimized.pdf"
        plt.savefig(conf_mat_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print("✅ Custom confusion matrix saved with rotated labels")
        
    except Exception as e:
        print(f"❌ Custom confusion matrix failed: {e}")

    # Custom predictions plot with rotated labels
    try:
        # Load predictions data
        with open(f"{output_dir}/{output_prefix}_test_pred_dict.pkl", 'rb') as f:
            pred_dict = pickle.load(f)
        
        # Create custom predictions visualization
        plt.figure(figsize=(16, 10))
        
        # Example: Create a heatmap of prediction probabilities
        if 'pred_proba' in pred_dict:
            # Get prediction probabilities for first few samples
            sample_probs = pred_dict['pred_proba'][:15]  # First 15 samples
            
            # Create a heatmap of prediction probabilities
            plt.imshow(sample_probs, cmap='viridis', aspect='auto')
            plt.colorbar(label='Prediction Probability')
            
            # Set labels with rotation
            plt.xticks(range(len(unique_conditions)), unique_conditions, rotation=45, ha='right')
            plt.yticks(range(min(15, len(sample_probs))), [f"Sample {i+1}" for i in range(min(15, len(sample_probs)))])
            
            plt.title("Prediction Probabilities - Hypoxia Conditions (OPTIMIZED)", fontsize=16, pad=20)
            plt.xlabel("Conditions", fontsize=14)
            plt.ylabel("Samples", fontsize=14)
            
            # Adjust layout
            plt.tight_layout()
            
            # Save the plot
            pred_plot_file = f"{output_dir}/{output_prefix}_custom_predictions_optimized.pdf"
            plt.savefig(pred_plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            print("✅ Custom predictions plot saved with rotated labels")
        else:
            print("⚠️  No prediction probabilities found for custom plot")
            
    except Exception as e:
        print(f"❌ Custom predictions plot failed: {e}")

    # Simple performance summary
    print(f"\n📊 PERFORMANCE SUMMARY:")
    print("="*50)

    if 'acc' in all_metrics_test:
        print(f"Accuracy: {all_metrics_test['acc']:.4f}")
    if 'macro_f1' in all_metrics_test:
        print(f"Macro F1-Score: {all_metrics_test['macro_f1']:.4f}")

    # Save enhanced summary to file
    summary_file = f"{output_dir}/{output_prefix}_enhanced_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("HYPOXIA CONDITION CLASSIFICATION - ENHANCED SUMMARY (OPTIMIZED)\n")
        f.write("="*70 + "\n\n")
        f.write(f"Dataset: {len(dataset)} cells\n")
        f.write(f"Conditions: {len(unique_conditions)}\n")  
        f.write(f"Training epochs: {CONFIG['num_train_epochs']} (OPTIMIZED)\n")
        f.write(f"Learning rate: {CONFIG['learning_rate']} (OPTIMIZED)\n")
        f.write(f"Batch size: {CONFIG['per_device_train_batch_size']} (OPTIMIZED)\n")
        f.write(f"Effective batch size: {CONFIG['per_device_train_batch_size'] * training_args['gradient_accumulation_steps']}\n")
        f.write(f"Freeze layers: {CONFIG['freeze_layers']} (OPTIMIZED)\n")
        f.write(f"Scheduler: {CONFIG['lr_scheduler_type']} (OPTIMIZED)\n\n")

        f.write("OVERALL METRICS:\n")
        f.write("-" * 30 + "\n")
        if 'acc' in all_metrics_test:
            f.write(f"Accuracy: {all_metrics_test['acc']:.4f}\n")
        if 'macro_f1' in all_metrics_test:
            f.write(f"Macro F1-Score: {all_metrics_test['macro_f1']:.4f}\n")

    print(f"✅ Enhanced summary saved to: {summary_file}")


# In[ ]:


# =============================================================================
# FINAL SUMMARY
# =============================================================================

print(f"\n🎉 OPTIMIZED ANALYSIS COMPLETE!")
print("="*80)
print(f"📁 Output: {output_dir}")
print(f"🎯 OPTIMIZATION SUMMARY:")
print(f"  • Increased epochs: {CONFIG['num_train_epochs']} (was 0.9)")
print(f"  • Reduced learning rate: {CONFIG['learning_rate']} (was 0.0008)")
print(f"  • Better scheduler: {CONFIG['lr_scheduler_type']} (was polynomial)")
print(f"  • Gradient accumulation: {training_args['gradient_accumulation_steps']} steps")
print(f"  • Enhanced regularization: label smoothing, gradient clipping")
print(f"  • Memory optimization: gradient checkpointing, conservative memory usage")

# List generated files
print(f"\n📊 Generated files:")
plot_files = []
data_files = []
model_files = []
summary_files = []

for file in os.listdir(output_dir):
    if output_prefix in file:
        if any(x in file for x in ['.png', '.pdf']):
            plot_files.append(file)
        elif any(x in file for x in ['.dataset', '.pkl']):
            data_files.append(file)
        elif any(x in file for x in ['summary']):
            summary_files.append(file)
        else:
            model_files.append(file)

if plot_files:
    print(f"\n🎨 Plot files ({len(plot_files)}):")
    for file in plot_files:
        print(f"  📊 {file}")

if data_files:
    print(f"\n📋 Data files ({len(data_files)}):")
    for file in data_files:
        print(f"  📄 {file}")

if summary_files:
    print(f"\n📈 Summary files ({len(summary_files)}):")
    for file in summary_files:
        print(f"  📊 {file}")

if model_files:
    print(f"\n🤖 Model files ({len(model_files)}):")
    for file in model_files:
        print(f"  🧠 {file}")

# Show best metrics
if 'all_metrics_test' in locals() and all_metrics_test:
    print(f"\n📈 BEST PERFORMANCE (OPTIMIZED):")
    accuracy_keys = [k for k in all_metrics_test.keys() if 'acc' in k.lower()]
    f1_keys = [k for k in all_metrics_test.keys() if 'f1' in k.lower()]

    if accuracy_keys:
        acc_key = accuracy_keys[0]
        print(f"  Accuracy: {all_metrics_test[acc_key]:.4f}")
    if f1_keys:
        f1_key = f1_keys[0]
        print(f"  F1-Score: {all_metrics_test[f1_key]:.4f}")

print(f"\n✅ OPTIMIZED training complete! Check {output_dir} for results.")
print("🎯 Expected improvements: Higher accuracy, better generalization, more stable training")
print("⏱️  Training time: ~1.5-2 hours (8 epochs) vs. ~10-15 minutes (0.9 epochs)")
print("="*80)


# In[ ]:


# Add this analysis
hypoxia_analysis = pd.crosstab(
    [x["condition"] for x in dataset],   
    [x.get("condition", "unknown") for x in dataset]  # if you have this field
)
print(hypoxia_analysis)



# In[ ]:




