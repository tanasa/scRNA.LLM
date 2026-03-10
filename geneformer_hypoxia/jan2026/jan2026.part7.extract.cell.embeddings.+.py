#!/usr/bin/env python
# coding: utf-8

# === Standard Library ===
import os
import sys

# Set CUDA memory allocator config BEFORE importing torch
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

import gc
import datetime
import warnings
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
from geneformer import EmbExtractor, Classifier, TranscriptomeTokenizer
import pickle

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

FORWARD_BATCH_SIZE = 12  # Minimum batch size - DO NOT increase
NPROC = 8  # Minimum parallelism - DO NOT increase
MAX_NCELLS =  1000  # ⚠️ ULTRA-LOW: Start with 100 cells, increase gradually if successful
SKIP_PLOTTING = False  # Keep plotting enabled

print("=" * 60)
print("MEMORY OPTIMIZATION SETTINGS:")
print(f"  FORWARD_BATCH_SIZE: {FORWARD_BATCH_SIZE}")
print(f"  NPROC: {NPROC}")
print(f"  MAX_NCELLS: {MAX_NCELLS if MAX_NCELLS else 'All cells'}")
print(f"  SKIP_PLOTTING: {SKIP_PLOTTING}")
print("=" * 60)

# ==============================
# Environment setup
# ==============================

env_path = "/home/btanasa/miniconda3/envs/torch-cu"
os.environ["CONDA_PREFIX"] = env_path
os.environ["CONDA_DEFAULT_ENV"] = "torch-cu"
os.environ["RETICULATE_PYTHON"] = f"{env_path}/bin/python"
print("✅ Environment set to torch-cu")
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

# Clear GPU cache
torch.cuda.empty_cache()
gc.collect()

# Check GPU memory
print(f"GPU memory allocated: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"GPU memory reserved: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")
print(f"CUDA memory after empty_cache: {torch.cuda.memory_allocated() / (1024**3):.2f} GB")

# Verify environment variable was set (should be set at top of file)
print("✅ Environment variable set!")
print(f"PYTORCH_CUDA_ALLOC_CONF: {os.environ.get('PYTORCH_CUDA_ALLOC_CONF')}")


# the model was saved in the location :
# /mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_out/250821141257/

token_dictionary_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl"

# Check file existence before loading
if not os.path.exists(token_dictionary_file):
    raise FileNotFoundError(f"Token dictionary file not found: {token_dictionary_file}")

try:
    with open(token_dictionary_file, "rb") as f:
        token_dict = pickle.load(f)
except Exception as e:
    raise RuntimeError(f"Failed to load token dictionary from {token_dictionary_file}: {e}")

print(type(token_dict))
print(len(token_dict))

# Show first 20 entries
for i, (gene, token_id) in enumerate(token_dict.items()):
    print(f"{gene} → {token_id}")
    if i >= 6:
        break

# Convert to DataFrame
token_df = pd.DataFrame(list(token_dict.items()), columns=["Ensembl_ID", "Token_ID"])
print(token_df)

# =============================================================================
# FOLDERS and INFORMATION
# =============================================================================

model_directory = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.model/260122101411_08epochs/260122_geneformer_cellClassifier_hypoxia_18K/ksplit1/checkpoint-210"
input_data_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.dataset"

# Verify the paths

if os.path.exists(token_dictionary_file):
    print(f"✅ File exists: {token_dictionary_file}")
else:
    print(f"❌ File does NOT exist: {token_dictionary_file}")

if os.path.exists(model_directory):
    print(f"✅ File exists: {model_directory}")
else:
    print(f"❌ File does NOT exist: {model_directory}")

if os.path.exists(input_data_file):
    print(f"✅ File exists: {input_data_file}")
else:
    print(f"❌ File does NOT exist: {input_data_file}")

# Load the dataset
print(f"\n📦 Loading dataset...")
try:
    dataset = load_from_disk(input_data_file)
    print(f"✅ Dataset loaded successfully!")
except Exception as e:
    print(f"❌ Failed to load dataset: {e}")
    exit()


# during the Tokenization, I have use the following labels
# cols = ['sample', 
# 'condition', 
# 'leiden', 
# 'orig.ident', 
# 'class', 
# 'pre_ingest_annots', 
# 'orig_ident', 
# 'labels_scanvi', 
# 'predictions_scanvi', 
# 'celltypes', 
# 'dataset', 
# 'fine_annotation', 
# 'barcode', 
# 'individual', 
# 'batch', 
# 'donor', 
# 'post_filter_annots', 
# 'post_filter_annots_number',
# 'cell_id',
# ]

# =============================================================================
# BASIC INFORMATION
# =============================================================================

print(f"\n📊 BASIC DATASET INFORMATION:")
print(f"   Total cells: {len(dataset):,}")
print(f"   Dataset type: {type(dataset)}")
print(f"   Number of columns: {len(dataset.column_names)}")

print(f"\n📋 METADATA :")
print(f"   Available metadata : {dataset.column_names}")

print(f"\n📋 FEATURE SCHEMA:")
for col_name in dataset.column_names:
    col_type = dataset.features[col_name]
    print(f"   {col_name}: {col_type}")


# Check number of rows and columns
print(dataset.num_rows)
print(dataset.column_names)

# View a sample row
# print(dataset[0])
# print(dataset[1])
# print(dataset[2])

# Data composition
# Sample dataset to avoid loading entire dataset into pandas (memory efficient)
sample_size = min(5000, len(dataset))
df = dataset.select(range(sample_size)).to_pandas()

# Unique medical conditions
print("🧍 Unique medical conditions:", df['condition'].nunique())
print(df['condition'].value_counts(), end="\n\n")

# Broad cell types
print("🩺 post_filter_annots:", df['post_filter_annots'].nunique())
print(df['post_filter_annots'].value_counts(), end="\n\n")

CONFIG = {

    # Dataset paths
    "input_data_file": "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.dataset",

    # Model paths

    "model_directory" : "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.model/260122101411_08epochs/260122_geneformer_cellClassifier_hypoxia_18K/ksplit1/checkpoint-210",

    "token_dictionary_file": "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl",

    # Dataset subsampling
    "target_size": 10000,  # Set to None to use full dataset
    "random_state": 73,

    # Model configuration
    "model_version": "V2",
    "freeze_layers": 2,

    # Note: forward_batch_size is set directly in EmbExtractor calls (4) for memory efficiency
    "nproc": 10,

    # Output
    "output_prefix": "jan2026.hypoxia_model_18K_extract_embs"
}

# =============================================================================
# SETUP AND VERIFICATION
# =============================================================================

print("🚀 Starting Hypoxia Geneformer Embedding Extraction")
print("=" * 60)

# Create timestamped output directory
current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
datestamp_min = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"

output_dir = f"/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.model.emb/{datestamp}"
os.makedirs(output_dir, exist_ok=True)

print(f"📁 Output directory: {output_dir}")
print(f"📅 Timestamp: {datestamp}")

# Verify all required files exist
print(f"\n🔍 Verifying files...")
required_files = [
    ("Input dataset", CONFIG["input_data_file"]),
    ("Model directory", CONFIG["model_directory"]),
    ("Token dictionary", CONFIG["token_dictionary_file"])
]

for name, path in required_files:
    if os.path.exists(path):
        print(f"✅ {name}: {path}")
    else:
        print(f"❌ {name}: {path}")
        print(f"   File not found - please check the path")
        # exit()

# === condition ===
# cell_states_to_model = {
#    "state_key": "condition",
#    "start_state": "2dHIE",
#    "goal_state": "non-hie",
#    "alt_states": []  # No alternative states with only 2 conditions
# }

print(output_dir)

# 🔬 Unique cell types: 

# Ast_Rgl
# Endo
# Epd
# Gli_IPC
# Imm_CGE
# Imm_CGE_LGE
# Imm_Exc_1
# Imm_Exc_2
# Imm_LGE
# Imm_MGE
# Imm_MGE_MAF_CRABP1
# Imm_MGE_MAF_NPY
# Inh_IPC
# Mcg
# Mcg_Neuronal
# MSN
# Myel_Olig
# OPC
# Prem_Olig
# Rbc

### checking the labels in the dataset : 

# Print column names
print("=" * 50)
print("AVAILABLE COLUMNS:")
print("=" * 50)
for col in dataset.column_names:
    print(f"  • {col}")

# Convert first few rows to DataFrame for better visualization
print("\n" + "=" * 50)
print("FIRST 5 SAMPLES:")
print("=" * 50)
df_samples = pd.DataFrame(dataset[:5])
print(df_samples.head())

# Show info about each column
print("\n" + "=" * 50)
print("COLUMN INFO:")
print("=" * 50)
for col in dataset.column_names:
    sample_value = dataset[0][col]
    print(f"{col}:")
    print(f"  Type: {type(sample_value)}")
    if isinstance(sample_value, (list, tuple)):
        print(f"  Length: {len(sample_value)}")
    print(f"  Sample: {sample_value if not isinstance(sample_value, list) or len(sample_value) < 10 else sample_value[:10]}")
    print()


#######################################################################################################
#######################################################################################################
# extract embedding in CELL model : # Layer 0
#######################################################################################################
####################################################################################################### 

print("EMBEDDING LAYER 0 ")
print("=" * 60)
print("⚠️  MEMORY TIP: If you get OOM errors, process ONE layer at a time")
print("⚠️  Comment out the 'EMBEDDING LAYER -1' section above to process only layer 0")
print("=" * 60)
import time

# Clear GPU cache before starting layer 0
print(f"🧹 Clearing GPU cache before layer 0 extraction...")
torch.cuda.empty_cache()
gc.collect()
print(f"GPU memory before extraction: {torch.cuda.memory_allocated() / 1024**3:.2f} GB / {torch.cuda.memory_reserved() / 1024**3:.2f} GB reserved")

try:
    from geneformer import EmbExtractor

    # Define perturbation parameters
    # cell_states_to_model = {
    #    "state_key": "condition",
    #    "start_state": "5d_HIE",
    #    "goal_state": "non-HIE",
    #    "alt_states": ["2mo_HIE"]
    # }


    cell_states_to_model = {
        "state_key": "condition",
        "start_state": "2dHIE",
        "goal_state": "non-hie",
        "alt_states": []  # No alternative states with only 2 conditions
    }

    # Use all disease states (no filtering)
    
    filter_data_dict = None
    # filter_data_dict = {"class":["Imm-MGE"]}
    # filter_data_dict={"class":["Astroglia","Imm-MGE","Imm-CGE"]} 

    # Setup embedding extractor for gene embeddings
    print(f"🔧 Initializing embedding extractor for cell embeddings...")
    
    # ⚠️ WARNING: Gene embeddings are EXTREMELY memory-intensive!
    
    # Gene embeddings = (num_cells × num_genes × embedding_dim)
    # Example: 1000 cells × 20,000 genes × 768 dims × 4 bytes = ~61 GB!

    if MAX_NCELLS is None:
        print("⚠️  CRITICAL WARNING: MAX_NCELLS is None - this WILL cause OOM for gene embeddings!")
        print("⚠️  Gene embeddings require VERY limited cell count. Setting to 100 as safety measure.")
        print("⚠️  To process all cells, run the script multiple times with different MAX_NCELLS ranges.")
        max_ncells_to_use = 1000
    else:
        max_ncells_to_use = MAX_NCELLS
        print(f"📊 Processing up to {max_ncells_to_use:,} cells for cell embeddings")
        print(f"⚠️  Memory estimate: ~{max_ncells_to_use * 20000 * 768 * 4 / (1024**3):.1f} GB for embeddings alone")
        print(f"⚠️  If this fails, reduce MAX_NCELLS further (try 50 or even 25)")
    
    embex = EmbExtractor(
        
        model_type = "CellClassifier", 
        num_classes = 2,
        
        filter_data = filter_data_dict,  # Filter by cell type
        # filter_data = None,  # Uncomment to use all cells

        # Memory optimization: Limit cells if needed (set to None for all cells)
        # For very large datasets, consider processing in chunks
        max_ncells = max_ncells_to_use,  # Use global setting from top of file
        
        emb_layer = 0,      # emb_layer{-1, 0}
                            # Layer 0 is more general representation
                            # The last layer (-1) is most specifically weighted to optimize the given learning objective.

        # emb_label = ["obs_names", "leiden", "class", "condition", "broad_class", "length"],  # All available labels
        # labels_to_plot = ["condition", "class", "broad_class", "leiden"],
        
        emb_label = ["cell_id", "sample", "condition", "leiden", "orig.ident", "class", "pre_ingest_annots", "labels_scanvi", "predictions_scanvi", "celltypes", "fine_annotation", "individual", "batch", "donor", "post_filter_annots", "post_filter_annots_number"], 
        # All available labels
        labels_to_plot = ["sample","condition","leiden","orig.ident","class","pre_ingest_annots","labels_scanvi","predictions_scanvi","celltypes","fine_annotation","individual","batch","donor","post_filter_annots","post_filter_annots_number"],


        forward_batch_size = FORWARD_BATCH_SIZE,  # Use global setting from top of file
        model_version = "V2",
        nproc = NPROC,  # Use global setting from top of file         

        # cell_emb_style = "mean_pool" # Method for summarizing gene embeddings if not using CLS token.
        # summary_stat = "exact_mean",  # Not needed for cell mode
        # summary_stat = "exact_median",
        # gene_emb_style = "mean_pool" # Method for summarizing gene embeddings.

        emb_mode = "cls"   
        # emb_mode = "gene",    
        # emb_mode = "cell" 
    )

    # Get gene embeddings with timing
    print(f"📊 Extracting cell embeddings...")
    
    # Final memory check and aggressive clearing before extraction
    print(f"🧹 Final GPU memory check before extraction...")
    torch.cuda.empty_cache()
    gc.collect()
    free_memory = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / (1024**3)
    print(f"   Free GPU memory: {free_memory:.2f} GB")
    if free_memory < 5.0:
        print(f"⚠️  WARNING: Only {free_memory:.2f} GB free - extraction may fail!")
        print(f"⚠️  Consider reducing MAX_NCELLS further or restarting Python to clear memory")
    
    start_time = time.time()
    
    # Extract gene embeddings (extract_embs doesn't take cell_states_to_model)
    # This is the memory-intensive operation - monitor closely
    try:
        cell_embs = embex.extract_embs(
            CONFIG["model_directory"],
            CONFIG["input_data_file"],
            output_dir,
            CONFIG["output_prefix"]
        )
    except torch.cuda.OutOfMemoryError as oom_error:
        print(f"❌ CUDA OOM during extraction!")
        print(f"   Current MAX_NCELLS: {max_ncells_to_use}")
        print(f"   Try reducing MAX_NCELLS to {max_ncells_to_use // 2} or even {max_ncells_to_use // 4}")
        print(f"   Or process one layer at a time by commenting out the other layer section")
        raise
    print(f"✅ Cell embeddings extracted")
    
    # Debug: Print structure of cell_embs to understand its format
    print(f"\n🔍 Inspecting cell_embs structure...")
    print(f"cell_embs type: {type(cell_embs)}")
    if isinstance(cell_embs, dict):
        print(f"cell_embs keys: {list(cell_embs.keys())}")
        for k, v in cell_embs.items():
            if hasattr(v, "shape"):
                print(f"  {k} shape: {v.shape}")
            elif isinstance(v, (list, tuple)):
                print(f"  {k} type: {type(v)}, length: {len(v)}")
            else:
                print(f"  {k} type: {type(v)}")
    elif hasattr(cell_embs, "shape"):
        print(f"cell_embs shape: {cell_embs.shape}")
    
    # Clear GPU cache immediately after extraction
    print(f"🧹 Clearing GPU cache after extraction...")
    torch.cuda.empty_cache()
    gc.collect()
    
    # Convert cell_embs to DataFrame and save
    print(f"\n💾 Converting embeddings to DataFrame and saving...")
    try:
        # Handle different return types from extract_embs()
        if isinstance(cell_embs, pd.DataFrame):
            df_cell_embs = cell_embs
        elif isinstance(cell_embs, dict):
            # If it's a dictionary, try to extract the embeddings
            # Prefer known keys in order: 'embs', 'gene', 'cell'
            if 'embs' in cell_embs:
                embs_data = cell_embs['embs']
                print(f"  Using key 'embs' for embeddings")
            elif 'gene' in cell_embs:
                embs_data = cell_embs['gene']
                print(f"  Using key 'gene' for embeddings")
            elif 'cell' in cell_embs:
                embs_data = cell_embs['cell']
                print(f"  Using key 'cell' for embeddings")
            elif len(cell_embs) > 0:
                # Fallback: Get first value if it's a dict of arrays
                # This is risky - prefer explicit keys above
                first_key = list(cell_embs.keys())[0]
                embs_data = cell_embs[first_key]
                print(f"  ⚠️  Warning: Using first key '{first_key}' - verify this is correct!")
            else:
                raise ValueError("cell_embs is empty")
            
            # Convert to numpy array if needed - ensure we move to CPU first
            if torch.is_tensor(embs_data):
                embs_data = embs_data.detach().cpu().numpy()
            elif hasattr(embs_data, 'detach'):
                embs_data = embs_data.detach().cpu().numpy()
            elif hasattr(embs_data, 'cpu'):
                embs_data = embs_data.cpu().numpy()
            elif not isinstance(embs_data, np.ndarray):
                embs_data = np.array(embs_data)
            
            # Create DataFrame
            if len(embs_data.shape) == 2:
                df_cell_embs = pd.DataFrame(embs_data)
            else:
                df_cell_embs = pd.DataFrame(embs_data.reshape(-1, embs_data.shape[-1]))
        else:
            # Assume it's an array-like object
            embs_data = np.array(cell_embs) if not isinstance(cell_embs, np.ndarray) else cell_embs
            if torch.is_tensor(embs_data):
                embs_data = embs_data.detach().cpu().numpy()
            elif hasattr(embs_data, 'detach'):
                embs_data = embs_data.detach().cpu().numpy()
            df_cell_embs = pd.DataFrame(embs_data)
        
        # Save to CSV
        output_filename = f"{CONFIG['output_prefix']}_cell_embeddings_emb_layer_0_emb_mode_cls.csv"
        output_filepath = os.path.join(output_dir, output_filename)
        df_cell_embs.to_csv(output_filepath, index=False)
        print(f"✅ Cell embeddings saved to: {output_filepath}")
        print(f"   Shape: {df_cell_embs.shape} (rows: cells, columns: embedding_dim)")
        
        # Clear DataFrame immediately after saving to free memory
        del df_cell_embs
        del embs_data
        gc.collect()
        
    except Exception as save_error:
        print(f"⚠️ Warning: Could not save embeddings to DataFrame/CSV: {save_error}")
        import traceback
        traceback.print_exc()

    end_time = time.time()
    execution_time = end_time - start_time
    
    print(f"✅ Cell embeddings extracted successfully")
    print(f"⏱️ Execution time: {execution_time:.2f} seconds ({execution_time/60:.2f} minutes)")
    
    # Clear GPU cache before plotting
    torch.cuda.empty_cache()
    gc.collect()
    
    # =========================================================================
    # PLOT UMAP (optional - skip if memory is tight)
    # =========================================================================
    # SKIP_PLOTTING is set at top of script - set to True to skip visualization and save memory
    
    if not SKIP_PLOTTING:
        print(f"\n📊 Generating UMAP plot...")
        embex.max_ncells = 1000  # Set max_ncells for plotting
        
        try:
            embex.plot_embs(
                embs=cell_embs,
                plot_style="umap",
                output_directory=output_dir,
                output_prefix=f"{CONFIG['output_prefix']}_cell_embeddings_emb_layer_0_emb_mode_cls_umap_2K_cells",
                max_ncells_to_plot=1000
            )
            print(f"✅ UMAP plot completed")
            torch.cuda.empty_cache()
            gc.collect()
        except Exception as plot_error:
            print(f"⚠️ Warning: UMAP plotting failed (may be memory issue): {plot_error}")
            torch.cuda.empty_cache()
            gc.collect()
    else:
        print(f"\n⏭️  Skipping UMAP plot (SKIP_PLOTTING=True)")
    
    # =========================================================================
    # PLOT HEATMAP (optional - skip if memory is tight)
    # =========================================================================
    if not SKIP_PLOTTING:
        print(f"\n📊 Generating Heatmap...")
        embex.max_ncells = 1000  # Update max_ncells for heatmap
        
        try:
            embex.plot_embs(
                embs=cell_embs,
                plot_style="heatmap",
                output_directory=output_dir,
                output_prefix=f"{CONFIG['output_prefix']}_cell_embeddings_emb_layer_0_emb_mode_cls_heatmap_2Kcells",
                max_ncells_to_plot=1000
            )
            print(f"✅ Heatmap completed")
            torch.cuda.empty_cache()
            gc.collect()
        except Exception as plot_error:
            print(f"⚠️ Warning: Heatmap plotting failed (may be memory issue): {plot_error}")
            torch.cuda.empty_cache()
            gc.collect()
    else:
        print(f"\n⏭️  Skipping Heatmap (SKIP_PLOTTING=True)")
    
    # Clean up after layer 0 processing
    print(f"🧹 Cleaning up after layer 0 processing...")
    del cell_embs
    del embex
    torch.cuda.empty_cache()
    gc.collect()
    print(f"GPU memory after cleanup: {torch.cuda.memory_allocated() / 1024**3:.2f} GB / {torch.cuda.memory_reserved() / 1024**3:.2f} GB reserved")

except Exception as e:
    print(f"❌ Cell embedding extraction failed: {e}")
    import traceback
    traceback.print_exc()
    # Clean up on error
    torch.cuda.empty_cache()
    gc.collect()
    # exit()



