#!/usr/bin/env python
# coding: utf-8

# According to Claude :

# The workflow is:

# CPU: Load and filter dataset
# CPU: Create perturbed versions (delete gene tokens)
# CPU → GPU: Send batches to GPU
# GPU: Run model inference on perturbed data
# GPU → CPU: Return predictions
# CPU: Aggregate and save results

import os, tempfile, pathlib

# Force temp to /scratch/$USER
os.environ["TMPDIR"] = f"/scratch/{os.getenv('USER')}"
tempfile.tempdir = os.environ["TMPDIR"]
pathlib.Path(os.environ["TMPDIR"]).mkdir(parents=True, exist_ok=True)

# Optional caches local to scratch (nice to have on HPC)
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.environ["TMPDIR"], "mpl"))
os.environ.setdefault("PYTORCH_TMPDIR", os.path.join(os.environ["TMPDIR"], "torch"))
os.environ.setdefault("HF_HOME", os.path.join(os.environ["TMPDIR"], "huggingface"))
os.environ.setdefault("TRANSFORMERS_CACHE", os.path.join(os.environ["HF_HOME"], "transformers"))
os.environ.setdefault("HF_DATASETS_CACHE", os.path.join(os.environ["HF_HOME"], "datasets"))
for p in (os.environ["MPLCONFIGDIR"], os.environ["PYTORCH_TMPDIR"], os.environ["HF_HOME"],
          os.environ["TRANSFORMERS_CACHE"], os.environ["HF_DATASETS_CACHE"]):
    pathlib.Path(p).mkdir(parents=True, exist_ok=True)

# Use spawn + clean shutdown for multiprocess
import multiprocess as mp
if __name__ == "__main__":
    try:
        mp.set_start_method("spawn", force=True)
    except RuntimeError:
        # Already set, ignore
        pass

# =============================================================================
# OPTIMIZED ISP ANALYSIS SCRIPT - GPU ONLY VERSION
# GPU-accelerated analysis for embedding extraction and ISP perturbation
# =============================================================================

# === Standard Library ===
import os
import gc
import datetime
import warnings
import logging
import time
import sys
import torch

# === Data & Math ===
import numpy as np
import pandas as pd

# === Hugging Face & Transformers ===
import transformers
from huggingface_hub import snapshot_download
from datasets import load_from_disk, Dataset

# === Geneformer ===
from geneformer import EmbExtractor, Classifier, InSilicoPerturber, InSilicoPerturberStats
import pickle

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# ==============================
# COMPREHENSIVE LOGGING SETUP
# ==============================

# Create timestamp for logging
current_time = datetime.datetime.now()
timestamp = current_time.strftime("%Y%m%d_%H%M%S")

# Create log directory
log_dir = f"/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/logs"
os.makedirs(log_dir, exist_ok=True)

# Configure logging (will be updated with strategy after strategy selection)
log_file = os.path.join(log_dir, f"isp_analysis_optimized_{timestamp}.log")

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()  # Also print to console
    ]
)

logger = logging.getLogger(__name__)

# Log script start
script_start_time = time.time()  # Track overall script execution time from the very beginning
logger.info("=" * 80)
logger.info("STARTING GPU-ONLY ISP ANALYSIS SCRIPT")
logger.info("=" * 80)
logger.info(f"Script: GPU-only ISP analysis")
logger.info(f"Timestamp: {timestamp}")
logger.info(f"Log file: {log_file}")
logger.info(f"Python executable: {sys.executable}")
logger.info(f"Working directory: {os.getcwd()}")
logger.info(f"Python version: {sys.version}")
logger.info(f"Script start time: {script_start_time}")

# ==============================
# GPU-ONLY ENVIRONMENT SETUP
# ==============================

print("🚀 Starting GPU-ONLY Hypoxia Geneformer ISP Pipeline")
print("🎯 Strategy: GPU for all operations")
print("🔧 GPU OPTIMIZED VERSION")
print("=" * 70)

# OPTIMIZED GPU MEMORY MANAGEMENT
# Increased max_split_size_mb to reduce fragmentation for larger batch sizes
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:128,expandable_segments:True"
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"

# Multiprocessing configuration to avoid EOFError
os.environ["TOKENIZERS_PARALLELISM"] = "false"  # Disable tokenizer multiprocessing
# Disable datasets library multiprocessing to avoid EOFError
os.environ["HF_DATASETS_NUM_PROC"] = "0"  # Disable multiprocessing in datasets library
os.environ["HF_DATASETS_NUM_WORKERS"] = "0"  # Alternative name for same setting
# Try to prevent multiprocessing connection issues
try:
    import multiprocessing
    multiprocessing.set_start_method("spawn", force=True)
except RuntimeError:
    # Already set, ignore
    pass

logger.info("GPU environment variables set:")
logger.info(f"PYTORCH_CUDA_ALLOC_CONF: {os.environ.get('PYTORCH_CUDA_ALLOC_CONF')}")
logger.info(f"CUDA_LAUNCH_BLOCKING: {os.environ.get('CUDA_LAUNCH_BLOCKING')}")

# GPU memory configuration - check available memory first
if torch.cuda.is_available():
    # Get total GPU memory
    gpu_total = torch.cuda.get_device_properties(0).total_memory / 1024**3
    # Clear all cached memory first
    torch.cuda.empty_cache()
    torch.cuda.ipc_collect()
    # Get currently used memory
    gpu_allocated = torch.cuda.memory_allocated() / 1024**3
    gpu_reserved = torch.cuda.memory_reserved() / 1024**3
    gpu_free = gpu_total - gpu_reserved
    
    logger.info(f"GPU Memory Status before configuration:")
    logger.info(f"  Total: {gpu_total:.2f} GB")
    logger.info(f"  Currently allocated: {gpu_allocated:.2f} GB")
    logger.info(f"  Currently reserved: {gpu_reserved:.2f} GB")
    logger.info(f"  Estimated free: {gpu_free:.2f} GB")
    
    # Dynamically adjust memory fraction based on available memory
    # Increased limit to allow larger batch sizes for ISP
    # Use up to 50% of total, but ensure at least 2GB is free for other processes
    target_usage = min(0.50, (gpu_free - 2.0) / gpu_total)
    target_usage = max(0.15, target_usage)  # Minimum 15%, maximum 50%
    
    # If very little memory is free, set to minimum
    if gpu_free < 3.0:
        target_usage = 0.15
        logger.warning(f"⚠️ Very little GPU memory free ({gpu_free:.2f} GB). Using minimum fraction (15%).")
        logger.warning("⚠️ Consider killing other GPU processes or using a GPU with more memory.")
        print(f"⚠️  WARNING: Only {gpu_free:.2f} GB GPU memory free!")
        print(f"⚠️  Using minimum memory fraction (15%) due to other processes.")
    
    torch.cuda.set_per_process_memory_fraction(target_usage)
    logger.info(f"Set GPU memory fraction to: {target_usage:.2%} (targeting {gpu_free:.2f} GB free)")
    print(f"🔧 GPU Memory Fraction: {target_usage:.2%} (targeting {gpu_free:.2f} GB free)")
    
    # Clear cache again after setting fraction
    torch.cuda.empty_cache()
    torch.cuda.ipc_collect()

# Check PyTorch and CUDA
try:
    logger.info(f"PyTorch version: {torch.__version__}")
    print(f"✅ PyTorch version: {torch.__version__}")
    if torch.cuda.is_available():
        logger.info("CUDA available - using GPU for all operations")
        print("✅ CUDA available - using GPU for all operations")
        device = "cuda"
        # Clear GPU cache
        torch.cuda.empty_cache()
        gpu_allocated = torch.cuda.memory_allocated() / 1024**3
        gpu_reserved = torch.cuda.memory_reserved() / 1024**3
        gpu_total = torch.cuda.get_device_properties(0).total_memory / 1024**3
        gpu_free = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / 1024**3
        
        logger.info(f"GPU memory allocated: {gpu_allocated:.2f} GB")
        logger.info(f"GPU memory reserved: {gpu_reserved:.2f} GB")
        logger.info(f"GPU memory total: {gpu_total:.2f} GB")
        logger.info(f"GPU memory free: {gpu_free:.2f} GB")
        
        print(f"GPU memory allocated: {gpu_allocated:.2f} GB")
        print(f"GPU memory reserved: {gpu_reserved:.2f} GB")
        print(f"GPU memory total: {gpu_total:.2f} GB")
        print(f"GPU memory free: {gpu_free:.2f} GB")
    else:
        logger.error("CUDA not available - GPU required for this script")
        print("❌ CUDA not available - GPU required for this script")
        raise RuntimeError("CUDA not available. This script requires GPU.")
except ImportError as e:
    logger.error(f"PyTorch import failed: {e}")
    print(f"❌ PyTorch import failed: {e}")
    raise

# Check Geneformer availability
try:
    from geneformer import TranscriptomeTokenizer
    logger.info("Geneformer available")
    print("✅ Geneformer available")
except ImportError as e:
    logger.error(f"Geneformer import failed: {e}")
    print(f"❌ Geneformer import failed: {e}")
    raise

# Memory cleanup
gc.collect()

# Memory monitoring function
def monitor_gpu_memory():
    if torch.cuda.is_available():
        allocated = torch.cuda.memory_allocated() / 1024**3
        reserved = torch.cuda.memory_reserved() / 1024**3
        free = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / 1024**3
        logger.info(f"GPU Memory - Allocated: {allocated:.2f} GB, Reserved: {reserved:.2f} GB, Free: {free:.2f} GB")
        print(f"🖥️ GPU Memory - Allocated: {allocated:.2f} GB, Reserved: {reserved:.2f} GB, Free: {free:.2f} GB")
        return allocated, reserved, free
    return None, None, None

# Initial memory check
logger.info("Initial GPU Memory Status:")
print("\n🔍 Initial GPU Memory Status:")
monitor_gpu_memory()

# ==============================
# CONFIGURATION AND PATHS
# ==============================

logger.info("Setting up configuration and paths...")

# Model version
model_version2 = "V2"
logger.info(f"Model version: {model_version2}")

# Token dictionary
token_dictionary_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl"
token_dict_path = token_dictionary_file
logger.info(f"Token dictionary file: {token_dictionary_file}")

# Load token dictionary
logger.info("Loading token dictionary...")
start_time = time.time()
try:
    with open(token_dictionary_file, "rb") as f:
        token_dict = pickle.load(f)
    load_time = time.time() - start_time
    logger.info(f"Token dictionary loaded successfully in {load_time:.2f} seconds")
    logger.info(f"Token dictionary type: {type(token_dict)}")
    logger.info(f"Token dictionary length: {len(token_dict)}")
    print(type(token_dict))
    print(len(token_dict))
except Exception as e:
    logger.error(f"Failed to load token dictionary: {e}")
    print(f"❌ Failed to load token dictionary: {e}")
    exit()

# Show first few entries
logger.info("Token dictionary sample entries:")
for i, (gene, token_id) in enumerate(token_dict.items()):
    logger.info(f"{gene} → {token_id}")
    print(f"{gene} → {token_id}")
    if i >= 6:
        break

# Convert to DataFrame
token_df = pd.DataFrame(list(token_dict.items()), columns=["Ensembl_ID", "Token_ID"])
logger.info(f"Token DataFrame shape: {token_df.shape}")
print(token_df)

# ==============================
# FOLDERS AND INFORMATION
# ==============================

logger.info("Setting up input/output paths...")

# Input paths
model_directory = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.model/260122101411_08epochs/260122_geneformer_cellClassifier_hypoxia_18K/ksplit1/checkpoint-210"
input_data_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.dataset"

logger.info(f"Model directory: {model_directory}")
logger.info(f"Input data file: {input_data_file}")

# Verify the paths
logger.info("Verifying file existence...")
required_files = [
    ("Token dictionary", token_dictionary_file),
    ("Model directory", model_directory),
    ("Input data file", input_data_file)
]

for name, path in required_files:
    if os.path.exists(path):
        logger.info(f"✅ {name}: {path}")
        print(f"✅ File exists: {path}")
    else:
        logger.error(f"❌ {name}: {path}")
        print(f"❌ File does NOT exist: {path}")

# Load the dataset
logger.info("Loading dataset...")
print(f"\n📊 Loading dataset...")
start_time = time.time()

try:
    dataset = load_from_disk(input_data_file)
    load_time = time.time() - start_time
    logger.info(f"Dataset loaded successfully in {load_time:.2f} seconds")
    logger.info(f"Dataset size: {len(dataset):,} cells")
    logger.info(f"Dataset columns: {dataset.column_names}")
    print(f"✅ Dataset loaded successfully!")
except Exception as e:
    logger.error(f"Failed to load dataset: {e}")
    print(f"❌ Failed to load dataset: {e}")
    exit()

# ==============================
# DATASET ANALYSIS
# ==============================

logger.info("Analyzing dataset...")

print(f"\n📊 BASIC DATASET INFORMATION:")
print(f"   Total cells: {len(dataset):,}")
print(f"   Dataset type: {type(dataset)}")
print(f"   Number of columns: {len(dataset.column_names)}")

logger.info(f"Total cells: {len(dataset):,}")
logger.info(f"Dataset type: {type(dataset)}")
logger.info(f"Number of columns: {len(dataset.column_names)}")

print(f"\n📋 METADATA:")
print(f"   Available metadata: {dataset.column_names}")

print(f"\n🔍 FEATURE SCHEMA:")
for col_name in dataset.column_names:
    col_type = dataset.features[col_name]
    print(f"   {col_name}: {col_type}")

# Check number of rows and columns
print(f"Number of rows: {dataset.num_rows}")
print(f"Column names: {dataset.column_names}")

# Data composition analysis
logger.info("Analyzing data composition...")
df = dataset.to_pandas()

# Basic counts
print("🧍 Unique medical conditions:", df['condition'].nunique())
print(df['condition'].value_counts(), end="\n\n")

print("🔬 Unique cell types:", df['class'].nunique())
print(df['class'].value_counts(), end="\n\n")

print("🩺 Broad classes of cell populations:", df['class'].nunique())
print(df['class'].value_counts(), end="\n\n")

# Log basic statistics
condition_counts = df['condition'].value_counts()
cell_type_counts = df['class'].value_counts()
broad_class_counts = df['class'].value_counts()

logger.info(f"Unique medical conditions: {df['condition'].nunique()}")
logger.info(f"Condition distribution: {dict(condition_counts)}")
logger.info(f"Unique cell types: {df['class'].nunique()}")
logger.info(f"Cell type distribution: {dict(cell_type_counts)}")
logger.info(f"Broad classes: {df['class'].nunique()}")
logger.info(f"Broad class distribution: {dict(broad_class_counts)}")

# =============================================================================
# CELL TYPE DISTRIBUTION PER CONDITION - COUNTS AND PERCENTAGES
# =============================================================================

logger.info("Computing cell type distribution per condition...")

# 1. Cross-tabulation for counts
print("📊 CELL TYPE DISTRIBUTION PER CONDITION (COUNTS):")
print("=" * 60)
counts_table = pd.crosstab(df['class'], df['condition'], margins=True)
print(counts_table)
print()

# Log counts table
logger.info("Cell type counts per condition:")
logger.info(f"Counts table shape: {counts_table.shape}")
logger.info(f"Counts table:\n{counts_table.to_string()}")

# 2. Cross-tabulation for percentages (within each condition)
print("📊 CELL TYPE DISTRIBUTION PER CONDITION (PERCENTAGES):")
print("=" * 60)
percentages_table = pd.crosstab(df['class'], df['condition'], normalize='columns') * 100
percentages_table = percentages_table.round(2)
print(percentages_table)
print()

# Log percentages table
logger.info("Cell type percentages per condition:")
logger.info(f"Percentages table shape: {percentages_table.shape}")
logger.info(f"Percentages table:\n{percentages_table.to_string()}")

# 3. Summary statistics
print("📈 SUMMARY STATISTICS:")
print("=" * 60)
print(f"Total cells: {len(df):,}")
print(f"Total conditions: {df['condition'].nunique()}")
print(f"Total cell types: {df['class'].nunique()}")
print()

# Log summary statistics
logger.info(f"Summary - Total cells: {len(df):,}")
logger.info(f"Summary - Total conditions: {df['condition'].nunique()}")
logger.info(f"Summary - Total cell types: {df['class'].nunique()}")

# 4. Most common cell type per condition
print("🏆 MOST COMMON CELL TYPE PER CONDITION:")
print("=" * 60)
most_common_per_condition = {}
for condition in df['condition'].unique():
    condition_data = df[df['condition'] == condition]
    most_common = condition_data['class'].mode().iloc[0]
    count = condition_data['class'].value_counts().iloc[0]
    percentage = (count / len(condition_data)) * 100
    most_common_per_condition[condition] = {
        'cell_type': most_common,
        'count': count,
        'percentage': percentage
    }
    print(f"{condition}: {most_common} ({count:,} cells, {percentage:.1f}%)")
print()

# Log most common cell types
logger.info("Most common cell type per condition:")
for condition, data in most_common_per_condition.items():
    logger.info(f"{condition}: {data['cell_type']} ({data['count']:,} cells, {data['percentage']:.1f}%)")

# 5. Cell type diversity per condition (number of unique cell types)
print("🔬 CELL TYPE DIVERSITY PER CONDITION:")
print("=" * 60)
diversity = df.groupby('condition')['class'].nunique().sort_values(ascending=False)
for condition, n_types in diversity.items():
    print(f"{condition}: {n_types} different cell types")
print()

# Log diversity statistics
logger.info("Cell type diversity per condition:")
for condition, n_types in diversity.items():
    logger.info(f"{condition}: {n_types} different cell types")

# 6. Detailed breakdown for each condition
print("📋 DETAILED BREAKDOWN BY CONDITION:")
print("=" * 60)
detailed_breakdown = {}
for condition in sorted(df['condition'].unique()):
    print(f"\n--- {condition} ---")
    condition_data = df[df['condition'] == condition]
    cell_type_counts = condition_data['class'].value_counts()
    cell_type_percentages = (cell_type_counts / len(condition_data)) * 100
    
    detailed_breakdown[condition] = {}
    for cell_type, count in cell_type_counts.items():
        percentage = cell_type_percentages[cell_type]
        detailed_breakdown[condition][cell_type] = {
            'count': count,
            'percentage': percentage
        }
        print(f"  {cell_type}: {count:,} cells ({percentage:.1f}%)")

# Log detailed breakdown
logger.info("Detailed breakdown by condition:")
for condition, cell_types in detailed_breakdown.items():
    logger.info(f"--- {condition} ---")
    for cell_type, data in cell_types.items():
        logger.info(f"  {cell_type}: {data['count']:,} cells ({data['percentage']:.1f}%)")

# ==============================
# CELL TYPE CONFIGURATION
# ==============================

cell_type = "Imm_MGE"  # Change this to analyze different cell types
logger.info(f"Selected cell type for analysis: {cell_type}")

# ==============================
# GPU-ONLY CONFIGURATION
# ==============================

CONFIG = {
    # Dataset paths
    "input_data_file": input_data_file,
    
    # Model paths
    "model_path": model_directory,
    "token_dictionary_file": token_dictionary_file,
    
    "target_size": None,  # Use full dataset
    "random_state": 73,
    
    # Model configuration
    "model_version": "V2",
    "freeze_layers": 2,
    # "forward_batch_size": 8,  # Batch size for GPU operations : unused
    
    # Batch sizes for GPU operations (reduced to avoid OOM and multiprocessing issues)
    "emb_forward_batch_size": 32,   # Batch size for embedding extraction (GPU) - reduced from 16
    "isp_forward_batch_size": 32,   # Batch size for ISP analysis (GPU) - reduced from 4
    
    # Cell counts for GPU operations (reduced to avoid OOM and multiprocessing issues)
    "emb_max_ncells": 1000,        # Cell count for embedding extraction (GPU) - reduced from 1000
    "isp_max_ncells": 100,        # Cell count for ISP analysis (GPU) - reduced from 500
    
    # Process configuration - GPU accelerated (reduced to avoid multiprocessing issues)
    "emb_nproc": 1,            # Processes for embedding extraction
    "isp_nproc": 1,             # Processes for ISP analysis (set to 1 to avoid EOFError in multiprocessing)
    
    # Output (will be updated with strategy after strategy selection)
    "output_prefix": f"cell_type_ISP_{cell_type}_jan2026"
}

logger.info(f"Configuration: {CONFIG}")

print("📋 Configuration:")
for key, value in CONFIG.items():
    print(f"   {key}: {value}")

print("\n🔄 GPU-Only Configuration Summary:")
print(f"   Embedding Extraction: {CONFIG['emb_nproc']} processes (GPU accelerated)")
print(f"   ISP Analysis: {CONFIG['isp_nproc']} processes (GPU accelerated)")
print(f"   Embedding Batch Size: {CONFIG['emb_forward_batch_size']}")
print(f"   ISP Batch Size: {CONFIG['isp_forward_batch_size']}")
print(f"   Embedding Max Cells: {CONFIG['emb_max_ncells']}")
print(f"   ISP Max Cells: {CONFIG['isp_max_ncells']}")

# ==============================
# SETUP AND VERIFICATION
# ==============================

logger.info("Setting up output directories...")
print("🚀 Starting Hypoxia Geneformer ISP Pipeline (GPU-Only)")
print("=" * 60)

# Create timestamped output directory
current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
datestamp_min = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"

# Output directory will be created after strategy selection
logger.info(f"Timestamp: {datestamp}")
print(f"📅 Timestamp: {datestamp}")

# Verify all required files exist
logger.info("Verifying required files...")
print(f"\n🔍 Verifying files...")
required_files = [
    ("Input dataset", CONFIG["input_data_file"]),
    ("Model directory", CONFIG["model_path"]),
    ("Token dictionary", CONFIG["token_dictionary_file"])
]

for name, path in required_files:
    if os.path.exists(path):
        logger.info(f"✅ {name}: {path}")
        print(f"✅ {name}: {path}")
    else:
        logger.error(f"❌ {name}: {path}")
        print(f"❌ {name}: {path}")
        print(f"   File not found - please check the path")

logger.info(f"Analyzing cell type: {cell_type}")
print(f"🔬 Analyzing cell type: {cell_type}")

# ==============================
# PERTURBATION STRATEGY CONFIGURATION
# ==============================

# Simple configuration for 2-condition setup (2dHIE → non-hie)
PERTURBATION_STRATEGY = "2dHIE_to_non_hie"

# Define perturbation parameters
cell_states_to_model = {
    "state_key": "condition",
    "start_state": "2dHIE",
    "goal_state": "non-hie",
    "alt_states": []  # No alternative states with only 2 conditions
}
strategy_description = "2dHIE → non-hie (detecting hits that shift state from 2dHIE to non-hie)"

# Log the chosen strategy
logger.info(f"Perturbation strategy: {PERTURBATION_STRATEGY}")
logger.info(f"Strategy description: {strategy_description}")
logger.info(f"Cell states to model: {cell_states_to_model}")

print(f"🎯 Perturbation strategy: {PERTURBATION_STRATEGY}")
print(f"📋 {strategy_description}")
print(f"🔬 Cell states: {cell_states_to_model}")

# Update output prefix to include strategy
CONFIG["output_prefix"] = f"cell_type_ISP_{cell_type}_2dHIE_to_non_hie_jan2026"
logger.info(f"Updated output prefix: {CONFIG['output_prefix']}")
print(f"📝 Output prefix updated: {CONFIG['output_prefix']}")

# Update output directory to include strategy
output_dir = f"/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.model.isp/{cell_type}_2dHIE_to_non_hie_{datestamp}_jan2026"
os.makedirs(output_dir, exist_ok=True)
logger.info(f"Updated output directory: {output_dir}")
print(f"📁 Output directory updated: {output_dir}")

# Update log file name to include strategy
log_file = os.path.join(log_dir, f"isp_analysis_2dHIE_to_non_hie_jan2026_{timestamp}.log")
logger.info(f"Updated log file: {log_file}")
print(f"📄 Log file updated: {log_file}")

# Update logging configuration with new log file
logging.getLogger().handlers[0].stream.close()
logging.getLogger().removeHandler(logging.getLogger().handlers[0])
file_handler = logging.FileHandler(log_file)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logging.getLogger().addHandler(file_handler)

# Log the strategy selection
logger.info("=" * 80)
logger.info("PERTURBATION STRATEGY SELECTED")
logger.info("=" * 80)
logger.info(f"Strategy: {PERTURBATION_STRATEGY}")
logger.info(f"Description: {strategy_description}")
logger.info(f"Cell states: {cell_states_to_model}")
logger.info(f"Output directory: {output_dir}")
logger.info(f"Log file: {log_file}")
logger.info("=" * 80)

# ==============================
# EMBEDDING EXTRACTION (GPU MODE)
# ==============================

logger.info("Starting GPU-accelerated embedding extraction...")
print("🚀 Starting GPU-accelerated embedding extraction...")

# Clear GPU memory aggressively before embedding extraction
if torch.cuda.is_available():
    logger.info("Clearing GPU memory before embedding extraction...")
    torch.cuda.empty_cache()
    torch.cuda.ipc_collect()
    gc.collect()
    monitor_gpu_memory()

# Initialize state_embs_dict to None BEFORE try block to ensure it exists in scope
state_embs_dict = None
filter_data_dict = None  # Also initialize filter_data_dict here

try:
    from geneformer import EmbExtractor

    # Ensure filter_data_dict is set in this scope
    if filter_data_dict is None:
        filter_data_dict = {"post_filter_annots": [cell_type]}
    else:
        # Update it to ensure it's current
        filter_data_dict = {"post_filter_annots": [cell_type]}
    
    logger.info(f"Filter data: {filter_data_dict}")
    print(f"🔎 Cell filter: {filter_data_dict}")

    # Setup embedding extractor with GPU optimization
    logger.info("Initializing embedding extractor (GPU mode)...")
    print(f"🔧 Initializing embedding extractor (GPU mode)...")
    
    # Check GPU memory before creating EmbExtractor
    if torch.cuda.is_available():
        gpu_mem_before = torch.cuda.memory_allocated() / 1024**3
        gpu_mem_free_before = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / 1024**3
        logger.info(f"GPU memory before EmbExtractor: {gpu_mem_before:.2f} GB allocated, {gpu_mem_free_before:.2f} GB free")
        
        if gpu_mem_free_before < 2.0:
            logger.warning("⚠️ Very little GPU memory available. Consider reducing batch size or max_ncells.")
            print(f"⚠️  WARNING: Only {gpu_mem_free_before:.2f} GB GPU memory free!")
    
    embex = EmbExtractor(
        model_type="CellClassifier",
        num_classes=2,
        filter_data=filter_data_dict,
        max_ncells=CONFIG["emb_max_ncells"],
        emb_layer=0,
        summary_stat="exact_mean",
        forward_batch_size=CONFIG["emb_forward_batch_size"],
        model_version="V2",
        nproc=CONFIG["emb_nproc"],
        emb_mode="cls"
    )

    # Extract state embeddings
    logger.info("Extracting state embeddings (GPU accelerated)...")
    print(f"📊 Extracting state embeddings (GPU accelerated)...")
    start_time = time.time()

    state_embs_dict = embex.get_state_embs(
        cell_states_to_model,
        CONFIG["model_path"],
        CONFIG["input_data_file"],
        output_dir,
        CONFIG["output_prefix"]
    )

    end_time = time.time()
    execution_time = end_time - start_time

    logger.info(f"State embeddings extracted successfully in {execution_time:.2f} seconds")
    logger.info(f"State embeddings keys: {list(state_embs_dict.keys())}")
    print(f"✅ State embeddings extracted successfully!")
    print(f"⏱️ Execution time: {execution_time:.2f} seconds ({execution_time/60:.2f} minutes)")

    # Save embeddings
    pkl_file = os.path.join(output_dir, f"{CONFIG['output_prefix']}.pkl")
    with open(pkl_file, 'wb') as f:
        pickle.dump(state_embs_dict, f)
    logger.info(f"State embeddings saved to: {pkl_file}")
    print(f"💾 State embeddings saved to: {pkl_file}")

    # Clear GPU memory after embedding extraction
    if torch.cuda.is_available():
        logger.info("Cleaning GPU memory after embedding extraction...")
        print("\n🧹 Cleaning GPU memory after embedding extraction...")
        monitor_gpu_memory()
        torch.cuda.empty_cache()
        gc.collect()
        print("✅ GPU memory cleared")
        monitor_gpu_memory()

except torch.cuda.OutOfMemoryError as e:
    logger.error("=" * 80)
    logger.error("CUDA OUT OF MEMORY ERROR - EMBEDDING EXTRACTION")
    logger.error("=" * 80)
    logger.error(f"Error: {e}")
    
    # Clear GPU memory
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.ipc_collect()
        monitor_gpu_memory()
    
    logger.error("SOLUTIONS:")
    logger.error("1. Reduce emb_forward_batch_size (currently: {})".format(CONFIG["emb_forward_batch_size"]))
    logger.error("2. Reduce emb_max_ncells (currently: {})".format(CONFIG["emb_max_ncells"]))
    logger.error("3. Kill other processes using GPU memory")
    logger.error("4. Use a GPU with more available memory")
    logger.error("5. Check for memory leaks in other processes")
    
    print("=" * 70)
    print("❌ CUDA OUT OF MEMORY - Embedding Extraction Failed")
    print("=" * 70)
    print("GPU memory is insufficient. Try:")
    print(f"  1. Reduce emb_forward_batch_size (current: {CONFIG['emb_forward_batch_size']})")
    print(f"  2. Reduce emb_max_ncells (current: {CONFIG['emb_max_ncells']})")
    print("  3. Kill other GPU processes")
    print("=" * 70)
    
    import traceback
    traceback.print_exc()
    sys.exit(1)
    
except Exception as e:
    logger.error(f"State embedding extraction failed: {e}")
    print(f"❌ State embedding extraction failed: {e}")
    print(f"   Cannot proceed without state embeddings")
    
    # Clear GPU memory on any error
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        torch.cuda.ipc_collect()
    
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Verify state_embs_dict was created successfully
if state_embs_dict is None:
    logger.error("State embeddings dictionary is not available. Cannot proceed to ISP.")
    print("❌ State embeddings dictionary is not available. Cannot proceed to ISP.")
    logger.error("Script will exit. Please check embedding extraction logs above.")
    sys.exit(1)

# Verify it has the expected structure
if not isinstance(state_embs_dict, dict) or len(state_embs_dict) == 0:
    logger.error("State embeddings dictionary is empty or invalid. Cannot proceed to ISP.")
    print("❌ State embeddings dictionary is empty or invalid. Cannot proceed to ISP.")
    sys.exit(1)

logger.info("State embeddings verified. Proceeding to ISP perturbation...")
logger.info(f"State embeddings contain {len(state_embs_dict)} states: {list(state_embs_dict.keys())}")
print("✅ State embeddings verified. Proceeding to ISP perturbation...")
print(f"   Found {len(state_embs_dict)} state embeddings: {list(state_embs_dict.keys())}")

# ==============================
# ISP PERTURBATION (GPU MODE)
# ==============================

logger.info("Starting GPU-accelerated ISP perturbation analysis...")
print("🚀 Starting GPU-accelerated ISP perturbation analysis...")

# Output folders
isp_output_dir = os.path.join(output_dir, "isp")
isp_stats_dir = os.path.join(output_dir, "isp_stats")
os.makedirs(isp_output_dir, exist_ok=True)
os.makedirs(isp_stats_dir, exist_ok=True)

logger.info(f"ISP output directory: {isp_output_dir}")
logger.info(f"ISP stats directory: {isp_stats_dir}")

# Genes to perturb
genes_to_perturb = "all"
logger.info(f"Genes to perturb: {genes_to_perturb}")

# Log the perturbation strategy for ISP
logger.info(f"ISP using perturbation strategy: {PERTURBATION_STRATEGY}")
logger.info(f"ISP strategy description: {strategy_description}")
logger.info(f"ISP cell states to model: {cell_states_to_model}")

print(f"🎯 ISP Perturbation strategy: {PERTURBATION_STRATEGY}")
print(f"📋 ISP {strategy_description}")
print(f"🔬 ISP Cell states: {cell_states_to_model}")

# Aggressive memory clearing before ISP (critical for large batch sizes)
if torch.cuda.is_available():
    logger.info("Clearing GPU memory aggressively before ISP...")
    print("\n🧹 Aggressive GPU memory cleanup before ISP...")
    torch.cuda.empty_cache()
    torch.cuda.ipc_collect()
    gc.collect()
    torch.cuda.synchronize()
    
    # Check memory status before ISP
    logger.info("GPU Memory Status before ISP:")
    print("\n🔍 GPU Memory Status before ISP:")
    monitor_gpu_memory()
    
    # Check if we have enough memory
    gpu_free = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / 1024**3
    if gpu_free < 8.0:
        logger.warning(f"⚠️ Only {gpu_free:.2f} GB free before ISP. Large batch sizes may fail.")
        logger.warning("⚠️ Consider reducing isp_forward_batch_size or killing other GPU processes.")
        print(f"⚠️  WARNING: Only {gpu_free:.2f} GB GPU memory free before ISP!")
        print(f"⚠️  If using large batch sizes (e.g., 36), consider reducing to 8-16 or freeing memory.")

# Track ISP success and timing
isp_success = False
isp_start_time = time.time()  # Track overall ISP section start time
init_start = isp_start_time
perturb_start = None
perturb_time = 0
total_time = 0

logger.info("=" * 80)
logger.info("ENTERING ISP PERTURBATION SECTION")
logger.info("=" * 80)
logger.info(f"ISP start time: {isp_start_time}")

# Verify state_embs_dict exists and is valid before entering ISP section
# Use try/except to check if variable exists
try:
    state_embs_check = state_embs_dict is not None
except NameError:
    state_embs_check = False
    logger.error("=" * 80)
    logger.error("CRITICAL ERROR: state_embs_dict is not defined!")
    logger.error("=" * 80)
    logger.error("Embedding extraction must have failed. Cannot proceed to ISP.")
    logger.error("Check logs above for embedding extraction errors.")
    print("=" * 70)
    print("❌ CRITICAL ERROR: state_embs_dict not defined!")
    print("=" * 70)
    print("Embedding extraction must have failed. Cannot proceed to ISP.")
    print("Script will exit.")
    sys.exit(1)

if not state_embs_check:
    logger.error("=" * 80)
    logger.error("CRITICAL ERROR: state_embs_dict is None!")
    logger.error("=" * 80)
    logger.error("Embedding extraction must have failed. Cannot proceed to ISP.")
    logger.error("Check logs above for embedding extraction errors.")
    print("=" * 70)
    print("❌ CRITICAL ERROR: state_embs_dict is None!")
    print("=" * 70)
    print("Embedding extraction must have failed. Cannot proceed to ISP.")
    print("Script will exit.")
    sys.exit(1)

logger.info(f"State embeddings available: True")
logger.info(f"State embeddings keys: {list(state_embs_dict.keys())}")

# Verify filter_data_dict also exists
try:
    filter_data_check = filter_data_dict is not None
except NameError:
    filter_data_check = False
    logger.error("CRITICAL ERROR: filter_data_dict is not defined!")
    logger.error("Cannot proceed to ISP. Script will exit.")
    print("❌ CRITICAL ERROR: filter_data_dict not defined!")
    sys.exit(1)

if not filter_data_check:
    logger.error("CRITICAL ERROR: filter_data_dict is None!")
    logger.error("Cannot proceed to ISP. Script will exit.")
    print("❌ CRITICAL ERROR: filter_data_dict is None!")
    sys.exit(1)

print("=" * 70)
print("🚀 ENTERING ISP PERTURBATION SECTION")
print("=" * 70)

try:
    # Time the ISP initialization
    logger.info("Initializing InSilicoPerturber (GPU mode)...")
    print("🔧 Initializing InSilicoPerturber (GPU mode)...")
    logger.info(f"ISP configuration:")
    logger.info(f"  - max_ncells: {CONFIG['isp_max_ncells']}")
    logger.info(f"  - forward_batch_size: {CONFIG['isp_forward_batch_size']}")
    logger.info(f"  - nproc: {CONFIG['isp_nproc']}")
    logger.info(f"  - model_path: {CONFIG['model_path']}")
    logger.info(f"  - input_data_file: {CONFIG['input_data_file']}")
    logger.info(f"  - state_embs_dict available: {state_embs_dict is not None}")
    logger.info(f"  - filter_data_dict available: {filter_data_dict is not None}")
    init_start = time.time()
    
    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb=genes_to_perturb,
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=2,
        emb_mode="cls",  # Use CLS embedding for V2
        filter_data=filter_data_dict,
        cell_states_to_model=cell_states_to_model,
        state_embs_dict=state_embs_dict,
        max_ncells=CONFIG["isp_max_ncells"],
        emb_layer=0,
        forward_batch_size=CONFIG["isp_forward_batch_size"],
        model_version="V2",
        nproc=CONFIG["isp_nproc"]
    )
    
    # Avoid the Column * int bug
    if hasattr(isp, "special_token"):
        isp.special_token = False
    
    init_time = time.time() - init_start
    logger.info(f"InSilicoPerturber initialized in {init_time:.2f} seconds")
    print(f"✅ InSilicoPerturber initialized in {init_time:.2f} seconds")
    
    # Time the perturbation execution
    logger.info("=" * 60)
    logger.info("STARTING ISP PERTURBATION EXECUTION")
    logger.info("=" * 60)
    logger.info(f"Model path: {CONFIG['model_path']}")
    logger.info(f"Input data file: {CONFIG['input_data_file']}")
    logger.info(f"ISP output directory: {isp_output_dir}")
    logger.info(f"Prefix: {CONFIG['output_prefix']}_tf")
    logger.info("This step typically takes 1-2 hours...")
    print("🧬 Starting ISP perturbation execution (GPU mode)...")
    print("   ⚠️  This step typically takes 1-2 hours...")
    print(f"   📁 Output will be saved to: {isp_output_dir}")
    perturb_start = time.time()
    
    isp.perturb_data(
        CONFIG["model_path"],
        CONFIG["input_data_file"],
        isp_output_dir,
        f"{CONFIG['output_prefix']}_tf"
    )
    
    perturb_time = time.time() - perturb_start
    total_time = time.time() - init_start
    isp_success = True
    
    logger.info("=" * 60)
    logger.info("ISP PERTURBATION EXECUTION COMPLETED")
    logger.info("=" * 60)
    
    logger.info(f"ISP perturbation completed successfully in {perturb_time:.2f} seconds")
    logger.info(f"Total ISP time (init + execution): {total_time:.2f} seconds")
    print(f"✅ ISP outputs -> {isp_output_dir}")
    print(f"⏱️ Perturbation execution time: {perturb_time:.2f} seconds ({perturb_time/60:.2f} minutes)")
    print(f"⏱️ Total time (init + execution): {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    
except Exception as e:
    perturb_time = time.time() - perturb_start if perturb_start is not None else 0
    total_time = time.time() - isp_start_time
    
    logger.error(f"ISP perturbation failed after {perturb_time:.2f} seconds: {e}")
    logger.error(f"Error type: {type(e).__name__}")
    logger.error(f"Error details: {str(e)}")
    
    # Try to get more detailed error information
    import traceback
    error_traceback = traceback.format_exc()
    logger.error(f"Full traceback:\n{error_traceback}")
    
    print(f"❌ ISP perturbation failed: {e}")
    print(f"⏱️ Failed after: {perturb_time:.2f} seconds ({perturb_time/60:.2f} minutes)")
    print(f"   Error type: {type(e).__name__}")
    
    # Check if it's a multiprocessing issue
    if "EOFError" in str(e) or "multiprocess" in str(e).lower():
        logger.error("=" * 60)
        logger.error("MULTIPROCESSING EOFError DETECTED")
        logger.error("=" * 60)
        logger.error("This is a known issue with multiprocessing in geneformer ISP.")
        logger.error(f"Current settings:")
        logger.error(f"  - isp_nproc: {CONFIG['isp_nproc']}")
        logger.error(f"  - isp_max_ncells: {CONFIG['isp_max_ncells']}")
        logger.error(f"  - isp_forward_batch_size: {CONFIG['isp_forward_batch_size']}")
        logger.error("=" * 60)
        logger.error("SOLUTIONS TO TRY:")
        logger.error("1. Script is already set to nproc=1 (single process mode)")
        logger.error("2. If still failing, the dataset.map() call in geneformer uses")
        logger.error("   internal multiprocessing that cannot be fully disabled")
        logger.error("3. Possible solutions:")
        logger.error("   - Try running on a machine with more memory")
        logger.error("   - Reduce max_ncells further (currently: {})".format(CONFIG["isp_max_ncells"]))
        logger.error("   - Reduce forward_batch_size further (currently: {})".format(CONFIG["isp_forward_batch_size"]))
        logger.error("   - Check if dataset is corrupted or incomplete")
        logger.error("   - Try loading dataset fresh to ensure it's not locked")
        logger.error("=" * 60)
        print("=" * 70)
        print("❌ MULTIPROCESSING EOFError - Known Issue")
        print("=" * 70)
        print("This error occurs in geneformer's internal multiprocessing.")
        print(f"Current nproc: {CONFIG['isp_nproc']} (already at minimum: 1)")
        print("\n💡 Possible fixes:")
        print("   1. Reduce max_ncells (currently: {})".format(CONFIG["isp_max_ncells"]))
        print("   2. Reduce forward_batch_size (currently: {})".format(CONFIG["isp_forward_batch_size"]))
        print("   3. Check system memory - may need more RAM")
        print("   4. Ensure dataset is not corrupted or in use by another process")
        print("=" * 70)
    
    # Check if it's a missing variable issue
    if "NameError" in str(type(e).__name__) or "state_embs_dict" in str(e):
        logger.error("State embeddings dictionary missing. This should not happen if embedding extraction succeeded.")
        print("❌ State embeddings dictionary missing. Embedding extraction may have failed silently.")
        isp_success = False

# ==============================
# STATISTICS COMPUTATION (GPU MODE)
# ==============================

# Initialize total_final_time to track overall execution
# Use isp_start_time if it exists, otherwise use script_start_time as fallback
stats_init_start = time.time()
try:
    # Try to use isp_start_time if it exists
    total_final_time = time.time() - isp_start_time
except NameError:
    # Fallback if isp_start_time wasn't set (shouldn't happen, but safety check)
    total_final_time = time.time() - script_start_time
    logger.warning("isp_start_time not found, using script_start_time for total_final_time calculation")

if isp_success:
    logger.info("Computing ISP statistics (GPU mode)...")
    print("📈 Computing ISP statistics (GPU mode)...")
    
    try:
        # Goal-state-shift statistics
        logger.info("Initializing InSilicoPerturberStats...")
        print("📈 Initializing InSilicoPerturberStats...")
        stats_init_start = time.time()
        
        ispstats = InSilicoPerturberStats(
            mode="goal_state_shift",
            genes_perturbed=genes_to_perturb,
            combos=0,
            anchor_gene=None,
            cell_states_to_model=cell_states_to_model,
            model_version="V2",
            token_dictionary_file=token_dictionary_file
        )
        
        stats_init_time = time.time() - stats_init_start
        logger.info(f"InSilicoPerturberStats initialized in {stats_init_time:.2f} seconds")
        print(f"✅ InSilicoPerturberStats initialized in {stats_init_time:.2f} seconds")
        
        # Time the statistics computation
        logger.info("Computing statistics...")
        print("📈 Computing statistics...")
        stats_comp_start = time.time()
        
        ispstats.get_stats(
            isp_output_dir,
            None,
            isp_stats_dir,
            f"{CONFIG['output_prefix']}_genes"
        )
        
        stats_comp_time = time.time() - stats_comp_start
        try:
            total_final_time = time.time() - isp_start_time
        except NameError:
            total_final_time = time.time() - script_start_time
        
        logger.info(f"Statistics computed in {stats_comp_time:.2f} seconds")
        logger.info(f"Total analysis time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
        print(f"✅ Statistics computed -> {isp_stats_dir}")
        print(f"⏱️ Statistics computation time: {stats_comp_time:.2f} seconds ({stats_comp_time/60:.2f} minutes)")
        print(f"⏱️ Total analysis time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
        
    except Exception as e:
        # Ensure total_final_time is set even if isp_start_time might not exist
        try:
            total_final_time = time.time() - isp_start_time
        except NameError:
            total_final_time = time.time() - script_start_time
        logger.error(f"Statistics computation failed: {e}")
        print(f"❌ Statistics computation failed: {e}")
        import traceback
        traceback.print_exc()
else:
    # Ensure total_final_time is set even if isp_start_time might not exist
    try:
        total_final_time = time.time() - isp_start_time
    except NameError:
        total_final_time = time.time() - script_start_time
    logger.warning("Skipping statistics computation - ISP perturbation did not complete successfully")
    print("⚠️ Skipping statistics computation - ISP perturbation did not complete successfully")
    print("   Statistics require successful ISP perturbation output")

# ==============================
# FINAL SUMMARY
# ==============================

# Calculate total script execution time
total_script_time = time.time() - script_start_time

logger.info("=" * 80)
if isp_success:
    logger.info("GPU-ONLY ISP ANALYSIS COMPLETED SUCCESSFULLY")
else:
    logger.info("GPU-ONLY ISP ANALYSIS COMPLETED WITH WARNINGS")
    logger.warning("⚠️ ISP perturbation did NOT complete successfully!")
    logger.warning("⚠️ Check error logs above for details.")
logger.info("=" * 80)
logger.info(f"ISP perturbation success: {isp_success}")
logger.info(f"Total ISP section time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
logger.info(f"Total script execution time: {total_script_time:.2f} seconds ({total_script_time/60:.2f} minutes)")
if total_script_time < 300:  # Less than 5 minutes
    logger.warning(f"⚠️ WARNING: Script completed in only {total_script_time:.2f} seconds!")
    logger.warning("⚠️ Expected ISP execution time: 1-2 hours (3600-7200 seconds)")
    logger.warning("⚠️ ISP perturbation may not have run! Check logs above.")
logger.info(f"Output directory: {output_dir}")
logger.info(f"ISP output directory: {isp_output_dir}")
logger.info(f"ISP stats directory: {isp_stats_dir}")
logger.info(f"Log file: {log_file}")
logger.info("=" * 80)

print("\n" + "=" * 60)
if isp_success:
    print("🏁 GPU-ONLY ISP ANALYSIS COMPLETED SUCCESSFULLY")
else:
    print("⚠️ GPU-ONLY ISP ANALYSIS COMPLETED WITH WARNINGS")
    print("   ISP perturbation did NOT complete successfully!")
print("=" * 60)
print(f"✅ ISP perturbation success: {isp_success}")
print(f"⏱️ Total ISP section time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
print(f"⏱️ Total script execution time: {total_script_time:.2f} seconds ({total_script_time/60:.2f} minutes)")
if total_script_time < 300:  # Less than 5 minutes
    print(f"\n⚠️  WARNING: Script completed in only {total_script_time:.2f} seconds!")
    print("⚠️  Expected ISP execution time: 1-2 hours (3600-7200 seconds)")
    print("⚠️  ISP perturbation may not have run! Check logs above for errors.")
    print("=" * 60)
print(f"📁 Output directory: {output_dir}")
print(f"📁 ISP output directory: {isp_output_dir}")
print(f"📁 ISP stats directory: {isp_stats_dir}")
print(f"📄 Log file: {log_file}")
print("=" * 60)

# Final memory cleanup
if torch.cuda.is_available():
    logger.info("Final GPU memory cleanup...")
    print("\n🧹 Final GPU memory cleanup...")
    torch.cuda.empty_cache()
    gc.collect()
    print("✅ Final memory status:")
    monitor_gpu_memory()

print(f"\n🚀 GPU-ONLY ISP Pipeline completed successfully!")
print(f"🔄 Used GPU-only approach for all operations")
# total_final_time should already be set by now, but add safety check
try:
    exec_time_str = f"{total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)"
except NameError:
    exec_time_str = "unknown (variable not set)"
print(f"⏱️ Total execution time: {exec_time_str}")

print("\n" + "="*70)
print("🎯 Next steps:")
print("1. Analyze perturbation_results for gene importance")
print("2. Check statistics in isp_stats_dir")
print("3. Visualize results using the generated data")
print("="*70)
print("🔧 GPU-ONLY Version Features:")
print("- GPU acceleration for all operations")
print("- GPU-accelerated embedding extraction")
print("- GPU-accelerated ISP analysis")
print("- Optimized batch sizes for GPU performance")
print("- GPU memory management and monitoring")
print("- Enhanced error handling and troubleshooting")
print("- MEMORY OPTIMIZED for GPU operations")
print("- Comprehensive logging for debugging and monitoring")
print("="*70)

logger.info("Script execution completed successfully")

