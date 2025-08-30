#!/usr/bin/env python
# coding: utf-8

# =============================================================================
# STABLE CPU-OPTIMIZED ISP PIPELINE for activated-NSCs in hypoxia
# GPU for embedding extraction, CPU for ISP analysis
# MEMORY OPTIMIZED VERSION
# =============================================================================

# In[1]: Standard Library Imports
import os
import gc
import datetime
import warnings
import time
import pickle
import glob
from collections import Counter
from pathlib import Path

# In[2]: Data & Math Libraries
import numpy as np
import pandas as pd
import torch

# In[3]: Hugging Face & Transformers
import transformers
from huggingface_hub import snapshot_download
from datasets import load_from_disk, Dataset

# In[4]: Geneformer
from geneformer import EmbExtractor, Classifier, InSilicoPerturber, InSilicoPerturberStats

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# In[5]: Environment Setup - Hybrid GPU/CPU Approach
print("🚀 Starting STABLE Hypoxia Geneformer ISP Pipeline")
print("🎯 Strategy: GPU for embedding extraction, CPU for ISP analysis")
print("🔧 MEMORY OPTIMIZED VERSION")
print("=" * 70)

#Environment Variables Control:

#Single-threaded operations within each process
#Library-level threading (NumPy, MKL, OpenBLAS)
#Prevents conflicts between different libraries
#Ensures stability of individual operations

#Process Configuration Controls:

#Multiple separate processes running in parallel
#Each process uses single-threaded operations
#Process-level parallelism for data handling
#Independent memory spaces for each process

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

# More conservative memory fraction for stability
torch.cuda.set_per_process_memory_fraction(0.25)  # Use only 25% of available VRAM
torch.cuda.empty_cache()  # Clear cache before and after

# CPU information
import psutil
cpu_count = os.cpu_count()
print(f"🖥️ CPU Cores Available: {cpu_count}")
print(f"🔄 Using Single Process Mode for ISP (multiprocessing disabled)")
print(f"🖥️ Total CPU Memory: {psutil.virtual_memory().total / (1024**3):.2f} GB")
print(f"🖥️ Available CPU Memory: {psutil.virtual_memory().available / (1024**3):.2f} GB")

# Set CPU threading to 1 for stability
torch.set_num_threads(1)
print(f"✅ PyTorch using {torch.get_num_threads()} CPU thread")

# Check PyTorch and CUDA
try:
    print(f"✅ PyTorch version: {torch.__version__}")
    if torch.cuda.is_available():
        print("✅ CUDA available - will use GPU for embedding extraction")
        device = "cuda"
        # Clear GPU cache
        torch.cuda.empty_cache()
        print(f"GPU memory allocated: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
        print(f"GPU memory reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")
        print(f"GPU memory total: {torch.cuda.get_device_properties(0).total_memory / 1024**3:.2f} GB")
        print(f"GPU memory free: {(torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / 1024**3:.2f} GB")
    else:
        print("🖥️ CUDA not available - will use CPU for everything")
        device = "cpu"
except ImportError as e:
    print(f"❌ PyTorch import failed: {e}")
    device = "cpu"

# Check Geneformer availability
try:
    from geneformer import TranscriptomeTokenizer
    print("✅ Geneformer available")
except ImportError as e:
    print(f"❌ Geneformer import failed: {e}")

# Memory cleanup
gc.collect()

# Memory monitoring function
def monitor_gpu_memory():
    if torch.cuda.is_available():
        allocated = torch.cuda.memory_allocated() / 1024**3
        reserved = torch.cuda.memory_reserved() / 1024**3
        free = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / 1024**3
        print(f"🖥️ GPU Memory - Allocated: {allocated:.2f} GB, Reserved: {reserved:.2f} GB, Free: {free:.2f} GB")
        return allocated, reserved, free
    return None, None, None

# Initial memory check
print("\n🔍 Initial GPU Memory Status:")
monitor_gpu_memory()

# In[6]: Model and Token Dictionary Setup
print("\n" + "="*50)
print("🔧 MODEL AND TOKEN SETUP")
print("="*50)

# Model paths
model_directory = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_out/250821141257/250821_geneformer_cellClassifier_hypoxia_condition_classifier_optimized/ksplit1/checkpoint-1400/"
input_data_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_oldlabel.dataset"
token_dictionary_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl"

# Load token dictionary
print("📚 Loading token dictionary...")
try:
    with open(token_dictionary_file, "rb") as f:
        token_dict = pickle.load(f)
    print(f"✅ Token dictionary loaded: {len(token_dict):,} genes")
except Exception as e:
    print(f"❌ Failed to load token dictionary: {e}")
    exit()

# Verify paths
print("\n🔍 Verifying file paths...")
required_files = [
    ("Input dataset", input_data_file),
    ("Model directory", model_directory),
    ("Token dictionary", token_dictionary_file)
]

for name, path in required_files:
    if os.path.exists(path):
        print(f"✅ {name}: {path}")
    else:
        print(f"❌ {name}: {path}")
        print(f"   File not found - please check the path")
        exit()

# In[7]: Dataset Loading and Information
print("\n" + "="*50)
print("📊 DATASET LOADING AND ANALYSIS")
print("="*50)

# Load dataset
print("📖 Loading dataset...")
try:
    dataset = load_from_disk(input_data_file)
    print(f"✅ Dataset loaded successfully!")
except Exception as e:
    print(f"❌ Failed to load dataset: {e}")
    exit()

# Basic dataset information
print(f"\n📊 BASIC DATASET INFORMATION:")
print(f"   Total cells: {len(dataset):,}")
print(f"   Dataset type: {type(dataset)}")
print(f"   Number of columns: {len(dataset.column_names)}")

# Data composition analysis
df = dataset.to_pandas()

print("\n🧍 MEDICAL CONDITIONS:")
print(f"   Unique conditions: {df['condition'].nunique()}")
print(df['condition'].value_counts())

print("\n🔬 CELL TYPES:")
print(f"   Unique cell types: {df['class'].nunique()}")
print(df['class'].value_counts())

print("\n🩺 BROAD CELL CLASSES:")
print(f"   Unique broad classes: {df['broad_class'].nunique()}")
print(df['broad_class'].value_counts())

# In[8]: Configuration
print("\n" + "="*50)
print("⚙️ CONFIGURATION")
print("="*50)

CONFIG = {
    # Dataset paths
    "input_data_file": input_data_file,
    "model_path": model_directory,
    "token_dictionary_file": token_dictionary_file,
    
    # Single process configuration for ISP
    "target_size": None,  # Use full dataset
    "random_state": 73,
    
    # Model configuration
    "model_version": "V2",
    "freeze_layers": 2,
    "forward_batch_size": 8,  # Very conservative batch size for memory stability
    
    # Separate batch sizes for different operations
    "emb_forward_batch_size": 16,  # Batch size for embedding extraction
    "isp_forward_batch_size": 4,   # Batch size for ISP analysis
    
    # Separate cell counts for different operations
    "emb_max_ncells": 1000,        # Cell count for embedding extraction
    "isp_max_ncells": 400,         # Cell count for ISP analysis
    
    # Process configuration - SEPARATE VARIABLES
    "emb_nproc": 20,            # CPU processes for embedding extraction (GPU mode)
    "isp_nproc": 1,             # CPU processes for ISP analysis (CPU mode)
    
    # Output
    "output_prefix": "hypoxia_2000cells_activated_NSC_ISP_stable"
}

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

# In[9]: Output Directory Setup
print("\n" + "="*50)
print("📁 OUTPUT DIRECTORY SETUP")
print("="*50)

# Create timestamped output directory
current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"

output_dir = f"/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_oldlabel.isp/{datestamp}"
os.makedirs(output_dir, exist_ok=True)

print(f"📁 Output directory: {output_dir}")
print(f"📅 Timestamp: {datestamp}")

# Create subdirectories
isp_output_dir = os.path.join(output_dir, "isp")
isp_stats_dir = os.path.join(output_dir, "isp_stats")
os.makedirs(isp_output_dir, exist_ok=True)
os.makedirs(isp_stats_dir, exist_ok=True)

print(f"📁 ISP output directory: {isp_output_dir}")
print(f"📁 ISP stats directory: {isp_stats_dir}")

# In[10]: Generate State Embeddings - USING GPU FOR SPEED
print("\n" + "="*50)
print("🔧 GENERATING STATE EMBEDDINGS (GPU MODE)")
print("=" * 50)

# Define perturbation parameters
cell_states_to_model = {
    "state_key": "condition",
    "start_state": "5d_HIE",
    "goal_state": "non-HIE",
    "alt_states": ["2mo_HIE"]
}

# Cell filtering
filter_data_dict = {"class": ["activated-NSCs"]}

print(f"🎯 Perturbation goal: 5d_HIE → non-HIE (with 2mo_HIE as alternative)")
print(f"🎯 Cell states: {cell_states_to_model}")
print(f"🔎 Cell filter: {filter_data_dict}")

# Setup embedding extractor with GPU optimization
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

print(f"✅ State embeddings extracted successfully!")
print(f"⏱️ Execution time: {execution_time:.2f} seconds ({execution_time/60:.2f} minutes)")

# Save embeddings
pkl_file = os.path.join(output_dir, f"{CONFIG['output_prefix']}.pkl")
with open(pkl_file, 'wb') as f:
    pickle.dump(state_embs_dict, f)
print(f"💾 State embeddings saved to: {pkl_file}")

# Clear GPU memory after embedding extraction
if torch.cuda.is_available():
    print("\n🧹 Cleaning GPU memory after embedding extraction...")
    monitor_gpu_memory()
    torch.cuda.empty_cache()
    gc.collect()
    print("✅ GPU memory cleared")
    monitor_gpu_memory()

# In[11]: In Silico Perturbation Setup
print("\n" + "="*50)
print("🧬 IN SILICO PERTURBATION SETUP (CPU MODE)")
print("=" * 50)

# Define perturbation parameters
cell_states_to_model = {
    "state_key": "condition",
    "start_state": "5d_HIE",
    "goal_state": "non-HIE",
    "alt_states": ["2mo_HIE"]
}

filter_data_dict = {"class": ["activated-NSCs"]}
genes_to_perturb = "all"

print(f"🎯 Perturbation goal: 5d_HIE → non-HIE (with 2mo_HIE as alternative)")
print(f"🎯 Cell states: {cell_states_to_model}")
print(f"🔎 Cell filter: {filter_data_dict}")
print(f"🧬 Genes to perturb: {genes_to_perturb}")
print(f"🖥️ ISP analysis will run on CPU for maximum stability")

# In[13]: Run In Silico Perturbation (CPU MODE)
print("\n" + "="*50)
print("🚀 RUNNING IN SILICO PERTURBATION (CPU MODE)")
print("=" * 50)

# Check memory status before ISP
if torch.cuda.is_available():
    print("\n🔍 GPU Memory Status before ISP:")
    monitor_gpu_memory()

try:
    # Time the ISP initialization
    print("🔧 Initializing InSilicoPerturber...")
    init_start = time.time()
    
    # Create ISP with single process and conservative settings for CPU stability
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
    print(f"✅ InSilicoPerturber initialized in {init_time:.2f} seconds")
    
    # Time the perturbation execution
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
    
    print(f"✅ ISP outputs -> {isp_output_dir}")
    print(f"⏱️ Perturbation execution time: {perturb_time:.2f} seconds ({perturb_time/60:.2f} minutes)")
    print(f"⏱️ Total time (init + execution): {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    
except Exception as e:
    print(f"❌ ISP execution failed: {e}")
    print(f"   Error details:")
    import traceback
    traceback.print_exc()
    
    print("\n💡 Troubleshooting suggestions:")
    print("1. Check if the model path is correct")
    print("2. Verify the dataset file is accessible")
    print("3. Ensure sufficient memory is available")
    print("4. Try running with even smaller batch sizes")
    exit()

# In[14]: ISP Analysis (CPU MODE)
print("\n" + "="*50)
print("📈 COMPUTING STATISTICS (CPU MODE)")
print("=" * 50)

try:
    # Goal-state-shift statistics
    print("🔧 Initializing InSilicoPerturberStats...")
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
    print(f"✅ InSilicoPerturberStats initialized in {stats_init_time:.2f} seconds")
    
    # Time the statistics computation
    print("📈 Computing statistics...")
    stats_comp_start = time.time()
    
    ispstats.get_stats(
        isp_output_dir,
        None,
        isp_stats_dir,
        f"{CONFIG['output_prefix']}_genes"
    )
    
    stats_comp_time = time.time() - stats_comp_start
    total_time_with_stats = time.time() - init_start
    
    print(f"✅ Statistics computed -> {isp_stats_dir}")
    print(f"⏱️ Statistics computation time: {stats_comp_time:.2f} seconds ({stats_comp_time/60:.2f} minutes)")
    print(f"⏱️ Total time (init + execution + stats): {total_time_with_stats:.2f} seconds ({total_time_with_stats/60:.2f} minutes)")
    
except Exception as e:
    print(f"❌ Statistics computation failed: {e}")
    import traceback
    traceback.print_exc()

# In[15]: Load and Display Results
print("\n" + "="*50)
print("📊 LOADING ISP RESULTS")
print("=" * 50)

# Look for pickle files in the isp_output_dir
pickle_files = glob.glob(os.path.join(isp_output_dir, "*_raw.pickle"))

if pickle_files:
    # Use the most recent file
    latest_pickle = max(pickle_files, key=os.path.getctime)
    print(f"✅ Found ISP results: {os.path.basename(latest_pickle)}")
    
    try:
        # Load automatically
        with open(latest_pickle, 'rb') as f:
            perturbation_results = pickle.load(f)
        
        print(f"\n📊 File contents:")
        print(f"Type: {type(perturbation_results)}")
        
        if isinstance(perturbation_results, dict):
            print(f"Keys: {list(perturbation_results.keys())}")
            
            # Explore the structure
            for key, value in perturbation_results.items():
                print(f"\n--- {key} ---")
                print(f"  Type: {type(value)}")
                
                if isinstance(value, dict):
                    print(f"  Dict keys: {list(value.keys())}")
                    for sub_key, sub_value in value.items():
                        print(f"    {sub_key}: {type(sub_value)}")
                        if hasattr(sub_value, 'shape'):
                            print(f"      Shape: {sub_value.shape}")
                        elif isinstance(sub_value, list):
                            print(f"      Length: {len(sub_value)}")
                            print(f"      First 5 values: {sub_value[:5]}")
                elif hasattr(value, 'shape'):
                    print(f"  Shape: {value.shape}")
                elif isinstance(value, list):
                    print(f"  Length: {len(value)}")
                    print(f"  First 10 values: {value[:10]}")
                    print(f"  Last 5 values: {value[-5:]}")
        else:
            print(f"Content: {perturbation_results}")
        
        print(f"\n🎯 Results loaded into 'perturbation_results' variable!")
        
    except Exception as e:
        print(f"❌ Failed to load results: {e}")
        perturbation_results = None
        
else:
    print(f"❌ No ISP results found in {isp_output_dir}")
    print(f"Available files:")
    try:
        for f in sorted(os.listdir(isp_output_dir)):
            print(f"  {f}")
    except Exception as e:
        print(f"Error listing directory: {e}")
    perturbation_results = None

# In[16]: Final Summary
print("\n" + "="*50)
print("🎉 PIPELINE COMPLETION SUMMARY")
print("=" * 50)

print(f"✅ Output directory: {output_dir}")
print(f"✅ ISP results: {isp_output_dir}")
print(f"✅ Statistics: {isp_stats_dir}")
print(f"✅ State embeddings: {len(state_embs_dict) if state_embs_dict else 0} states")
print(f"✅ Perturbation results: {'Loaded' if perturbation_results else 'Not found'}")

if perturbation_results:
    print(f"✅ Results structure: {type(perturbation_results)}")
    if isinstance(perturbation_results, dict):
        print(f"✅ Result keys: {list(perturbation_results.keys())}")

# Final memory cleanup
if torch.cuda.is_available():
    print("\n🧹 Final GPU memory cleanup...")
    torch.cuda.empty_cache()
    gc.collect()
    print("✅ Final memory status:")
    monitor_gpu_memory()

print(f"\n🚀 STABLE ISP Pipeline completed successfully!")
print(f"🔄 Used hybrid approach: GPU for embeddings, CPU for ISP")
print(f"⏱️ Total execution time: {total_time_with_stats:.2f} seconds ({total_time_with_stats/60:.2f} minutes)")

print("\n" + "="*70)
print("🎯 Next steps:")
print("1. Analyze perturbation_results for gene importance")
print("2. Check statistics in isp_stats_dir")
print("3. Visualize results using the generated data")
print("="*70)
print("🔧 STABLE Version Features:")
print("- Hybrid GPU/CPU approach for optimal performance")
print("- GPU acceleration for embedding extraction")
print("- CPU stability for ISP analysis (single process)")
print("- Conservative batch sizes for reliability")
print("- Environment variable overrides for stability")
print("- Enhanced error handling and troubleshooting")
print("- MEMORY OPTIMIZED with separate batch size and cell count controls")
print("="*70)
