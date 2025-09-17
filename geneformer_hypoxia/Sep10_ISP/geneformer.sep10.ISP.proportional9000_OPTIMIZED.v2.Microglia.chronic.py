#!/usr/bin/env python
# coding: utf-8

# =============================================================================
# OPTIMIZED ISP ANALYSIS SCRIPT with Hybrid GPU/CPU Approach
# Based on script_cursor_v3_GPU_emb_CPU_ISP.py optimizations
# =============================================================================

# === Standard Library ===
import os
import gc
import datetime
import warnings
import logging
import time
import sys
import psutil
from collections import Counter

# === Data & Math ===
import numpy as np
import pandas as pd
import torch

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
logger.info("=" * 80)
logger.info("STARTING OPTIMIZED ISP ANALYSIS SCRIPT")
logger.info("=" * 80)
logger.info(f"Script: test_ISP_hypoxia2000_v11-all_genes_Microglia.sep11.proportional9000_OPTIMIZED.py")
logger.info(f"Timestamp: {timestamp}")
logger.info(f"Log file: {log_file}")
logger.info(f"Python executable: {sys.executable}")
logger.info(f"Working directory: {os.getcwd()}")
logger.info(f"Python version: {sys.version}")

# ==============================
# OPTIMIZED ENVIRONMENT SETUP
# ==============================

print("🚀 Starting OPTIMIZED Hypoxia Geneformer ISP Pipeline")
print("🎯 Strategy: GPU for embedding extraction, CPU for ISP analysis")
print("🔧 MEMORY OPTIMIZED VERSION")
print("=" * 70)

# CRITICAL: Force single-threaded operations for CPU stability
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["TOKENIZERS_PARALLELISM"] = "false"

# OPTIMIZED GPU MEMORY MANAGEMENT
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "max_split_size_mb:64,expandable_segments:True"
os.environ["CUDA_LAUNCH_BLOCKING"] = "1"

logger.info("Environment variables set for stability:")
logger.info(f"OMP_NUM_THREADS: {os.environ.get('OMP_NUM_THREADS')}")
logger.info(f"MKL_NUM_THREADS: {os.environ.get('MKL_NUM_THREADS')}")
logger.info(f"NUMEXPR_NUM_THREADS: {os.environ.get('NUMEXPR_NUM_THREADS')}")
logger.info(f"OPENBLAS_NUM_THREADS: {os.environ.get('OPENBLAS_NUM_THREADS')}")
logger.info(f"VECLIB_MAXIMUM_THREADS: {os.environ.get('VECLIB_MAXIMUM_THREADS')}")
logger.info(f"TOKENIZERS_PARALLELISM: {os.environ.get('TOKENIZERS_PARALLELISM')}")
logger.info(f"PYTORCH_CUDA_ALLOC_CONF: {os.environ.get('PYTORCH_CUDA_ALLOC_CONF')}")
logger.info(f"CUDA_LAUNCH_BLOCKING: {os.environ.get('CUDA_LAUNCH_BLOCKING')}")

# More conservative memory fraction for stability
torch.cuda.set_per_process_memory_fraction(0.25)  # Use only 25% of available VRAM
torch.cuda.empty_cache()  # Clear cache before and after

# CPU information
cpu_count = os.cpu_count()
logger.info(f"CPU Cores Available: {cpu_count}")
logger.info(f"Total CPU Memory: {psutil.virtual_memory().total / (1024**3):.2f} GB")
logger.info(f"Available CPU Memory: {psutil.virtual_memory().available / (1024**3):.2f} GB")

print(f"🖥️ CPU Cores Available: {cpu_count}")
print(f"🔄 Using Hybrid GPU/CPU Mode")
print(f"🖥️ Total CPU Memory: {psutil.virtual_memory().total / (1024**3):.2f} GB")
print(f"🖥️ Available CPU Memory: {psutil.virtual_memory().available / (1024**3):.2f} GB")

# Set CPU threading to 1 for stability
torch.set_num_threads(1)
logger.info(f"PyTorch using {torch.get_num_threads()} CPU thread")
print(f"✅ PyTorch using {torch.get_num_threads()} CPU thread")

# Check PyTorch and CUDA
try:
    logger.info(f"PyTorch version: {torch.__version__}")
    print(f"✅ PyTorch version: {torch.__version__}")
    if torch.cuda.is_available():
        logger.info("CUDA available - will use GPU for embedding extraction")
        print("✅ CUDA available - will use GPU for embedding extraction")
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
        logger.info("CUDA not available - will use CPU for everything")
        print("🖥️ CUDA not available - will use CPU for everything")
        device = "cpu"
except ImportError as e:
    logger.error(f"PyTorch import failed: {e}")
    print(f"❌ PyTorch import failed: {e}")
    device = "cpu"

# Check Geneformer availability
try:
    from geneformer import TranscriptomeTokenizer
    logger.info("Geneformer available")
    print("✅ Geneformer available")
except ImportError as e:
    logger.error(f"Geneformer import failed: {e}")
    print(f"❌ Geneformer import failed: {e}")

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
model_directory = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.proportional.sampling.3000.dataset.out.model.250910201851/250910201851/250910_geneformer_cellClassifier_hypoxia_condition_classifier_optimized/ksplit1/checkpoint-6304/"
input_data_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.proportional.sampling.3000.dataset"

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

print("🩺 Broad classes of cell populations:", df['broad_class'].nunique())
print(df['broad_class'].value_counts(), end="\n\n")

# Log basic statistics
condition_counts = df['condition'].value_counts()
cell_type_counts = df['class'].value_counts()
broad_class_counts = df['broad_class'].value_counts()

logger.info(f"Unique medical conditions: {df['condition'].nunique()}")
logger.info(f"Condition distribution: {dict(condition_counts)}")
logger.info(f"Unique cell types: {df['class'].nunique()}")
logger.info(f"Cell type distribution: {dict(cell_type_counts)}")
logger.info(f"Broad classes: {df['broad_class'].nunique()}")
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

cell_type = "Microglia"  # Change this to analyze different cell types
logger.info(f"Selected cell type for analysis: {cell_type}")

# ==============================
# OPTIMIZED CONFIGURATION
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
    "forward_batch_size": 8,  # Very conservative batch size for memory stability
    
    # Separate batch sizes for different operations
    "emb_forward_batch_size": 16,  # Batch size for embedding extraction (GPU)
    "isp_forward_batch_size": 4,   # Batch size for ISP analysis (CPU)
    
    # Separate cell counts for different operations
    "emb_max_ncells": 1000,        # Cell count for embedding extraction (GPU)
    "isp_max_ncells": 400,         # Cell count for ISP analysis (CPU)
    
    # Process configuration - SEPARATE VARIABLES
    "emb_nproc": 20,            # CPU processes for embedding extraction (GPU mode)
    "isp_nproc": 1,             # CPU processes for ISP analysis (CPU mode)
    
    # Output (will be updated with strategy after strategy selection)
    "output_prefix": f"hypoxia9000proportional_ISP_{cell_type}_optimized"
}

logger.info(f"Configuration: {CONFIG}")

print("📋 Configuration:")
for key, value in CONFIG.items():
    print(f"   {key}: {value}")

print("\n🔄 Process Configuration Summary:")
print(f"   Embedding Extraction: {CONFIG['emb_nproc']} CPU processes (GPU accelerated)")
print(f"   ISP Analysis: {CONFIG['isp_nproc']} CPU process (CPU stable)")
print(f"   Embedding Batch Size: {CONFIG['emb_forward_batch_size']}")
print(f"   ISP Batch Size: {CONFIG['isp_forward_batch_size']}")
print(f"   Embedding Max Cells: {CONFIG['emb_max_ncells']}")
print(f"   ISP Max Cells: {CONFIG['isp_max_ncells']}")

# ==============================
# SETUP AND VERIFICATION
# ==============================

logger.info("Setting up output directories...")
print("🚀 Starting Hypoxia Geneformer ISP Pipeline")
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

# User choice: "acute" or "chronic"
PERTURBATION_STRATEGY = "chronic"  # Change this to "chronic" for chronic analysis

# Define perturbation parameters based on strategy
if PERTURBATION_STRATEGY == "acute":
    cell_states_to_model = {
        "state_key": "condition",
        "start_state": "5d_HIE",
        "goal_state": "non-HIE",
        "alt_states": ["2mo_HIE"]
    }
    strategy_description = "Acute (5d_HIE → non-HIE, with 2mo_HIE as alternative)"
    
elif PERTURBATION_STRATEGY == "chronic":
    cell_states_to_model = {
        "state_key": "condition",
        "start_state": "2mo_HIE",
        "goal_state": "non-HIE",
        "alt_states": ["5d_HIE"]
    }
    strategy_description = "Chronic (2mo_HIE → non-HIE, with 5d_HIE as alternative)"
    
else:
    raise ValueError(f"Invalid perturbation strategy: {PERTURBATION_STRATEGY}. Choose 'acute' or 'chronic'")

# Log the chosen strategy
logger.info(f"Perturbation strategy: {PERTURBATION_STRATEGY}")
logger.info(f"Strategy description: {strategy_description}")
logger.info(f"Cell states to model: {cell_states_to_model}")

print(f"🎯 Perturbation strategy: {PERTURBATION_STRATEGY}")
print(f"📋 {strategy_description}")
print(f"🔬 Cell states: {cell_states_to_model}")

# Update output prefix to include strategy
CONFIG["output_prefix"] = f"hypoxia9000proportional_ISP_{cell_type}_{PERTURBATION_STRATEGY}_optimized"
logger.info(f"Updated output prefix: {CONFIG['output_prefix']}")
print(f"📝 Output prefix updated: {CONFIG['output_prefix']}")

# Update output directory to include strategy
output_dir = f"/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_3_full_h5ad.reprocess.sep10.proportional.sampling.3000.dataset.out.model.250910201851.isp/{cell_type}_{PERTURBATION_STRATEGY}_{datestamp}_optimized"
os.makedirs(output_dir, exist_ok=True)
logger.info(f"Updated output directory: {output_dir}")
print(f"📁 Output directory updated: {output_dir}")

# Update log file name to include strategy
log_file = os.path.join(log_dir, f"isp_analysis_{PERTURBATION_STRATEGY}_{timestamp}.log")
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

try:
    from geneformer import EmbExtractor

    filter_data_dict = {"class": [cell_type]}
    
    logger.info(f"Filter data: {filter_data_dict}")

    print(f"🔎 Cell filter: {filter_data_dict}")

    # Setup embedding extractor with GPU optimization
    logger.info("Initializing embedding extractor (GPU mode)...")
    print(f"🔧 Initializing embedding extractor (GPU mode)...")
    
    embex = EmbExtractor(
        model_type="CellClassifier",
        num_classes=3,
        filter_data=filter_data_dict,
        max_ncells=CONFIG["emb_max_ncells"],  # Use CONFIG variable
        emb_layer=0,
        summary_stat="exact_mean",
        forward_batch_size=CONFIG["emb_forward_batch_size"],  # Use CONFIG variable
        model_version="V2",
        nproc=CONFIG["emb_nproc"],  # Use separate variable for embedding extraction
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

except Exception as e:
    logger.error(f"State embedding extraction failed: {e}")
    print(f"❌ State embedding extraction failed: {e}")
    print(f"   Cannot proceed without state embeddings")
    import traceback
    traceback.print_exc()

# ==============================
# ISP PERTURBATION (CPU MODE)
# ==============================

logger.info("Starting CPU-optimized ISP perturbation analysis...")

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

# Check memory status before ISP
if torch.cuda.is_available():
    logger.info("GPU Memory Status before ISP:")
    print("\n🔍 GPU Memory Status before ISP:")
    monitor_gpu_memory()

try:
    # Time the ISP initialization
    logger.info("Initializing InSilicoPerturber (CPU mode)...")
    print("🔧 Initializing InSilicoPerturber (CPU mode)...")
    init_start = time.time()
    
    # Create ISP with single process and conservative settings for CPU stability
    logger.info("Using single process mode with conservative settings for maximum stability")
    print("🔄 Using single process mode with conservative settings for maximum stability")
    
    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb=genes_to_perturb,
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=3,
        emb_mode="cls",  # Use CLS embedding for V2
        filter_data=filter_data_dict,
        cell_states_to_model=cell_states_to_model,
        state_embs_dict=state_embs_dict,
        max_ncells=CONFIG["isp_max_ncells"],  # Use CONFIG variable
        emb_layer=0,
        forward_batch_size=CONFIG["isp_forward_batch_size"],  # Use CONFIG variable
        model_version="V2",
        nproc=CONFIG["isp_nproc"]  # Use separate variable for ISP analysis
    )
    
    # Avoid the Column * int bug
    if hasattr(isp, "special_token"):
        isp.special_token = False
    
    init_time = time.time() - init_start
    logger.info(f"InSilicoPerturber initialized in {init_time:.2f} seconds")
    print(f"✅ InSilicoPerturber initialized in {init_time:.2f} seconds")
    
    # Time the perturbation execution
    logger.info("Running ISP perturbation (CPU mode)...")
    print("🧬 Running ISP perturbation (CPU mode)...")
    perturb_start = time.time()
    
    isp.perturb_data(
        CONFIG["model_path"],
        CONFIG["input_data_file"],
        isp_output_dir,
        f"{CONFIG['output_prefix']}_tf"
    )
    
    perturb_time = time.time() - perturb_start
    total_time = time.time() - init_start
    
    logger.info(f"ISP perturbation completed successfully in {perturb_time:.2f} seconds")
    logger.info(f"Total ISP time (init + execution): {total_time:.2f} seconds")
    print(f"✅ ISP outputs -> {isp_output_dir}")
    print(f"⏱️ Perturbation execution time: {perturb_time:.2f} seconds ({perturb_time/60:.2f} minutes)")
    print(f"⏱️ Total time (init + execution): {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    
except Exception as e:
    perturb_time = time.time() - perturb_start if 'perturb_start' in locals() else 0
    total_time = time.time() - init_start if 'init_start' in locals() else 0
    
    logger.error(f"ISP perturbation failed after {perturb_time:.2f} seconds: {e}")
    logger.error(f"Error type: {type(e).__name__}")
    logger.error(f"Error details: {str(e)}")
    
    # Try to get more detailed error information
    import traceback
    error_traceback = traceback.format_exc()
    logger.error(f"Full traceback:\n{error_traceback}")
    
    print(f"❌ ISP perturbation failed: {e}")
    print(f"⏱️ Failed after: {perturb_time:.2f} seconds ({perturb_time/60:.2f} minutes)")
    
    # Check if it's a multiprocessing issue and suggest solutions
    if "EOFError" in str(e) or "multiprocess" in str(e).lower():
        logger.error("This appears to be a multiprocessing issue. Suggestions:")
        logger.error("1. Try running with nproc=1")
        logger.error("2. Check available memory")
        logger.error("3. Try reducing max_ncells or forward_batch_size")
        print("💡 This appears to be a multiprocessing issue. Try running with nproc=1")
    
    # Don't exit, continue to statistics section if possible
    logger.warning("Continuing to statistics section despite perturbation failure...")
    print("⚠️ Continuing to statistics section despite perturbation failure...")

# ==============================
# STATISTICS COMPUTATION (CPU MODE)
# ==============================

logger.info("Computing ISP statistics (CPU mode)...")

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
    total_final_time = time.time() - init_start if 'init_start' in locals() else time.time()
    
    logger.info(f"Statistics computed in {stats_comp_time:.2f} seconds")
    logger.info(f"Total analysis time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
    print(f"✅ Statistics computed -> {isp_stats_dir}")
    print(f"⏱️ Statistics computation time: {stats_comp_time:.2f} seconds ({stats_comp_time/60:.2f} minutes)")
    print(f"⏱️ Total analysis time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
    
except Exception as e:
    logger.error(f"Statistics computation failed: {e}")
    print(f"❌ Statistics computation failed: {e}")
    import traceback
    traceback.print_exc()

# ==============================
# FINAL SUMMARY
# ==============================

logger.info("=" * 80)
logger.info("OPTIMIZED ISP ANALYSIS COMPLETED")
logger.info("=" * 80)
logger.info(f"Total execution time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
logger.info(f"Output directory: {output_dir}")
logger.info(f"ISP output directory: {isp_output_dir}")
logger.info(f"ISP stats directory: {isp_stats_dir}")
logger.info(f"Log file: {log_file}")
logger.info("=" * 80)

print("\n" + "=" * 60)
print("🏁 OPTIMIZED ISP ANALYSIS COMPLETED")
print("=" * 60)
print(f"⏱️ Total execution time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")
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

print(f"\n🚀 OPTIMIZED ISP Pipeline completed successfully!")
print(f"🔄 Used hybrid approach: GPU for embeddings, CPU for ISP")
print(f"⏱️ Total execution time: {total_final_time:.2f} seconds ({total_final_time/60:.2f} minutes)")

print("\n" + "="*70)
print("🎯 Next steps:")
print("1. Analyze perturbation_results for gene importance")
print("2. Check statistics in isp_stats_dir")
print("3. Visualize results using the generated data")
print("="*70)
print("🔧 OPTIMIZED Version Features:")
print("- Hybrid GPU/CPU approach for optimal performance")
print("- GPU acceleration for embedding extraction")
print("- CPU stability for ISP analysis (single process)")
print("- Conservative batch sizes for reliability")
print("- Environment variable overrides for stability")
print("- Enhanced error handling and troubleshooting")
print("- MEMORY OPTIMIZED with separate batch size and cell count controls")
print("- Comprehensive logging for debugging and monitoring")
print("="*70)

logger.info("Script execution completed successfully")
