
#!/usr/bin/env python
# coding: utf-8
# =============================================================================
# Geneformer Cell Embedding Extraction - v3
# =============================================================================

# === Standard Library ===
import os
import sys
import gc
import time
import datetime
import warnings
import pickle
import traceback
from collections import Counter

# Set CUDA memory allocator config BEFORE importing torch
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# === Data & Math ===
import numpy as np
import pandas as pd

# === PyTorch ===
import torch

# === Hugging Face & Transformers ===
import transformers
from huggingface_hub import snapshot_download
from datasets import load_from_disk, Dataset

# === Geneformer ===
from geneformer import EmbExtractor, Classifier, TranscriptomeTokenizer

# =============================================================================
# GLOBAL SETTINGS  (edit these to change behavior)
# =============================================================================
FORWARD_BATCH_SIZE = 12   # Reduce if GPU OOM
NPROC              = 8    # Parallelism
MAX_NCELLS         = 1000 # Set to None to use all cells
SKIP_PLOTTING      = False

print("=" * 60)
print("MEMORY OPTIMIZATION SETTINGS:")
print(f"  FORWARD_BATCH_SIZE : {FORWARD_BATCH_SIZE}")
print(f"  NPROC              : {NPROC}")
print(f"  MAX_NCELLS         : {MAX_NCELLS if MAX_NCELLS else 'All cells'}")
print(f"  SKIP_PLOTTING      : {SKIP_PLOTTING}")
print("=" * 60)

# =============================================================================
# ENVIRONMENT SETUP
# =============================================================================
env_path = "/home/btanasa/miniconda3/envs/torch-cu"
os.environ["CONDA_PREFIX"]       = env_path
os.environ["CONDA_DEFAULT_ENV"]  = "torch-cu"
os.environ["RETICULATE_PYTHON"]  = f"{env_path}/bin/python"

print("✅ Environment set to torch-cu")
print("   Python executable :", sys.executable)
print("   CONDA_PREFIX      :", os.environ.get("CONDA_PREFIX"))
print("   PYTORCH_CUDA_ALLOC_CONF:", os.environ.get('PYTORCH_CUDA_ALLOC_CONF'))
print("✅ PyTorch available :", torch.__version__)
print("✅ CUDA available    :", torch.cuda.is_available())

# Clear GPU cache
torch.cuda.empty_cache()
gc.collect()
print(f"   GPU memory allocated : {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
print(f"   GPU memory reserved  : {torch.cuda.memory_reserved()  / 1024**3:.2f} GB")

# =============================================================================
# TIMESTAMP  (must be defined BEFORE filename / output_dir)
# =============================================================================
current_date  = datetime.datetime.now()
datestamp     = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
datestamp_min = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"

Cell_Type = "All"
filename  = f"{datestamp}_{Cell_Type}"

print(f"\n📅 Timestamp : {datestamp}")
print(f"📁 Filename  : {filename}")

# =============================================================================
# CONFIGURATION
# =============================================================================
CONFIG = {
    # Paths
    "input_data_file"      : "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.dataset",
    "model_directory"      : "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.model/260122101411_08epochs/260122_geneformer_cellClassifier_hypoxia_18K/ksplit1/checkpoint-210",
    "token_dictionary_file": "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl",
#    "token_dictionary_file": "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.redownloaded.pkl",
    # Dataset subsampling
    "target_size"          : 10000,   # None = use full dataset
    "random_state"         : 73,
    # Model configuration
    "model_version"        : "V2",
    "freeze_layers"        : 2,
    "nproc"                : 10,
    # Output
    "output_prefix"        : "jan2026.hypoxia_model_18K_extract_embs",
}

# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================
output_dir = f"/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia_newobj_jan2026/SCANVI_INTEGRATED_FILTERED.geneformer_compatible.dataset.downsample.18K.model.emb/{filename}"
os.makedirs(output_dir, exist_ok=True)
print(f"\n📁 Output directory: {output_dir}")

# =============================================================================
# VERIFY REQUIRED FILES
# =============================================================================
print(f"\n🔍 Verifying files...")
required_files = [
    ("Input dataset"    , CONFIG["input_data_file"]),
    ("Model directory"  , CONFIG["model_directory"]),
    ("Token dictionary" , CONFIG["token_dictionary_file"]),
]
all_ok = True
for name, path in required_files:
    if os.path.exists(path):
        print(f"✅ {name}: {path}")
    else:
        print(f"❌ {name} NOT FOUND: {path}")
        all_ok = False

if not all_ok:
    print("❌ One or more required files are missing. Please check paths.")
    sys.exit(1)

# =============================================================================
# LOAD TOKEN DICTIONARY
# =============================================================================
print(f"\n📖 Loading token dictionary...")
try:
    with open(CONFIG["token_dictionary_file"], "rb") as f:
        token_dict = pickle.load(f)
    print(f"✅ Token dictionary loaded: {len(token_dict):,} entries")
    token_df = pd.DataFrame(list(token_dict.items()), columns=["Ensembl_ID", "Token_ID"])
    print(token_df.head(7))
except Exception as e:
    raise RuntimeError(f"Failed to load token dictionary: {e}")

# =============================================================================
# LOAD DATASET
# =============================================================================
print(f"\n📦 Loading dataset from: {CONFIG['input_data_file']}")
try:
    dataset = load_from_disk(CONFIG["input_data_file"])
    print(f"✅ Dataset loaded: {len(dataset):,} cells, {len(dataset.column_names)} columns")
except Exception as e:
    print(f"❌ Failed to load dataset: {e}")
    sys.exit(1)

# =============================================================================
# DATASET INFORMATION
# =============================================================================
print(f"\n📊 DATASET INFORMATION:")
print(f"   Total cells    : {len(dataset):,}")
print(f"   Columns        : {dataset.column_names}")

print("\n📋 FEATURE SCHEMA:")
for col_name in dataset.column_names:
    print(f"   {col_name}: {dataset.features[col_name]}")

# Sample for quick inspection
sample_size = min(5000, len(dataset))
df = dataset.select(range(sample_size)).to_pandas()
print(f"\n🧍 Conditions ({df['condition'].nunique()}):")
print(df['condition'].value_counts())
print(f"\n🩺 post_filter_annots ({df['post_filter_annots'].nunique()}):")
print(df['post_filter_annots'].value_counts())

print("\n" + "=" * 50)
print("FIRST 5 SAMPLES:")
print("=" * 50)
print(pd.DataFrame(dataset[:5]).head())

print("\n" + "=" * 50)
print("COLUMN INFO:")
print("=" * 50)
for col in dataset.column_names:
    sample_value = dataset[0][col]
    print(f"{col}:")
    print(f"  Type  : {type(sample_value)}")
    if isinstance(sample_value, (list, tuple)):
        print(f"  Length: {len(sample_value)}")
    print(f"  Sample: {sample_value if not isinstance(sample_value, list) or len(sample_value) < 10 else sample_value[:10]}")
    print()

# =============================================================================
# EMBEDDING EXTRACTION - LAYER 0 (CLS mode)
# =============================================================================
print("=" * 60)
print("EMBEDDING EXTRACTION - LAYER 0 (CLS mode)")
print("=" * 60)

print(f"🧹 Clearing GPU cache before extraction...")
torch.cuda.empty_cache()
gc.collect()
print(f"   GPU memory: {torch.cuda.memory_allocated() / 1024**3:.2f} GB allocated / {torch.cuda.memory_reserved() / 1024**3:.2f} GB reserved")

try:
    from geneformer import EmbExtractor

    cell_states_to_model = {
        "state_key"  : "condition",
        "start_state": "2dHIE",
        "goal_state" : "non-hie",
        "alt_states" : []
    }

    filter_data_dict = None
    # filter_data_dict = {"class": ["Imm-MGE"]}

    max_ncells_to_use = MAX_NCELLS if MAX_NCELLS is not None else 1000
    print(f"📊 Processing up to {max_ncells_to_use:,} cells")

    print(f"🔧 Initializing EmbExtractor...")
    embex = EmbExtractor(
        model_type        = "CellClassifier",
        num_classes       = 2,
        filter_data       = filter_data_dict,
        max_ncells        = max_ncells_to_use,
        emb_layer         = 0,
        emb_label         = [
            "cell_id", "sample", "condition", "leiden", "orig.ident",
            "class", "pre_ingest_annots", "labels_scanvi", "predictions_scanvi",
            "celltypes", "fine_annotation", "individual", "batch", "donor",
            "post_filter_annots", "post_filter_annots_number"
        ],
        labels_to_plot    = [
            "sample", "condition", "leiden", "orig.ident", "class",
            "pre_ingest_annots", "labels_scanvi", "predictions_scanvi",
            "celltypes", "fine_annotation", "individual", "batch", "donor",
            "post_filter_annots", "post_filter_annots_number"
        ],
        forward_batch_size = FORWARD_BATCH_SIZE,
        model_version      = "V2",
        nproc              = NPROC,
        emb_mode           = "cls",
    )

    # Final memory check before extraction
    torch.cuda.empty_cache()
    gc.collect()
    free_memory = (torch.cuda.get_device_properties(0).total_memory - torch.cuda.memory_reserved()) / (1024**3)
    print(f"   Free GPU memory before extraction: {free_memory:.2f} GB")
    if free_memory < 5.0:
        print(f"⚠️  WARNING: Only {free_memory:.2f} GB free - extraction may fail!")

    # Extract embeddings
    print(f"📊 Extracting cell embeddings...")
    start_time = time.time()

    try:
        cell_embs = embex.extract_embs(
            CONFIG["model_directory"],
            CONFIG["input_data_file"],
            output_dir,
            CONFIG["output_prefix"]
        )
    except torch.cuda.OutOfMemoryError:
        print(f"❌ CUDA OOM! Try reducing MAX_NCELLS to {max_ncells_to_use // 2}")
        torch.cuda.empty_cache()
        gc.collect()
        raise

    end_time       = time.time()
    execution_time = end_time - start_time
    print(f"✅ Cell embeddings extracted in {execution_time:.2f}s ({execution_time/60:.2f} min)")

    # Inspect structure
    print(f"\n🔍 cell_embs type: {type(cell_embs)}")
    if isinstance(cell_embs, dict):
        for k, v in cell_embs.items():
            shape = v.shape if hasattr(v, "shape") else (f"len={len(v)}" if isinstance(v, (list,tuple)) else type(v))
            print(f"   {k}: {shape}")
    elif hasattr(cell_embs, "shape"):
        print(f"   shape: {cell_embs.shape}")

    # Clear GPU cache after extraction
    torch.cuda.empty_cache()
    gc.collect()

    # ─────────────────────────────────────────────
    # SAVE EMBEDDINGS TO CSV
    # ─────────────────────────────────────────────
    print(f"\n💾 Saving embeddings...")
    try:
        if isinstance(cell_embs, pd.DataFrame):
            df_cell_embs = cell_embs

        elif isinstance(cell_embs, dict):
            for key in ['embs', 'gene', 'cell']:
                if key in cell_embs:
                    embs_data = cell_embs[key]
                    print(f"   Using key '{key}'")
                    break
            else:
                first_key = list(cell_embs.keys())[0]
                embs_data = cell_embs[first_key]
                print(f"⚠️  Fallback: using first key '{first_key}'")

            if torch.is_tensor(embs_data) or hasattr(embs_data, 'detach'):
                embs_data = embs_data.detach().cpu().numpy()
            elif hasattr(embs_data, 'cpu'):
                embs_data = embs_data.cpu().numpy()
            elif not isinstance(embs_data, np.ndarray):
                embs_data = np.array(embs_data)

            df_cell_embs = pd.DataFrame(
                embs_data if len(embs_data.shape) == 2
                else embs_data.reshape(-1, embs_data.shape[-1])
            )
        else:
            embs_data = np.array(cell_embs) if not isinstance(cell_embs, np.ndarray) else cell_embs
            if torch.is_tensor(embs_data) or hasattr(embs_data, 'detach'):
                embs_data = embs_data.detach().cpu().numpy()
            df_cell_embs = pd.DataFrame(embs_data)

        output_filename = f"{CONFIG['output_prefix']}_cell_embeddings_emb_layer_0_emb_mode_cls.csv"
        output_filepath = os.path.join(output_dir, output_filename)
        df_cell_embs.to_csv(output_filepath, index=False)
        print(f"✅ Saved: {output_filepath}")
        print(f"   Shape: {df_cell_embs.shape} (cells × embedding_dim)")

    except Exception as save_error:
        print(f"⚠️  Could not save embeddings: {save_error}")
        traceback.print_exc()
    finally:
        for var in ['df_cell_embs', 'embs_data']:
            if var in locals():
                del locals()[var]
        gc.collect()

    # ─────────────────────────────────────────────
    # UMAP PLOT
    # ─────────────────────────────────────────────
    if not SKIP_PLOTTING:
        print(f"\n📊 Generating UMAP plot...")
        embex.max_ncells = 1000
        try:
            embex.plot_embs(
                embs              = cell_embs,
                plot_style        = "umap",
                output_directory  = output_dir,
                output_prefix     = f"{CONFIG['output_prefix']}_cell_embeddings_emb_layer_0_emb_mode_cls_umap_1K_cells",
                max_ncells_to_plot= 1000,
            )
            print(f"✅ UMAP plot completed")
        except Exception as e:
            print(f"⚠️  UMAP plotting failed: {e}")
        finally:
            torch.cuda.empty_cache()
            gc.collect()

        # ─────────────────────────────────────────
        # HEATMAP PLOT
        # ─────────────────────────────────────────
        print(f"\n📊 Generating Heatmap...")
        embex.max_ncells = 1000
        try:
            embex.plot_embs(
                embs              = cell_embs,
                plot_style        = "heatmap",
                output_directory  = output_dir,
                output_prefix     = f"{CONFIG['output_prefix']}_cell_embeddings_emb_layer_0_emb_mode_cls_heatmap_1K_cells",
                max_ncells_to_plot= 1000,
            )
            print(f"✅ Heatmap completed")
        except Exception as e:
            print(f"⚠️  Heatmap plotting failed: {e}")
        finally:
            torch.cuda.empty_cache()
            gc.collect()
    else:
        print(f"\n⏭️  Skipping plots (SKIP_PLOTTING=True)")

    # ─────────────────────────────────────────────
    # CLEANUP
    # ─────────────────────────────────────────────
    print(f"\n🧹 Cleaning up...")
    del cell_embs, embex
    torch.cuda.empty_cache()
    gc.collect()
    print(f"   GPU memory after cleanup: {torch.cuda.memory_allocated() / 1024**3:.2f} GB / {torch.cuda.memory_reserved() / 1024**3:.2f} GB reserved")

except Exception as e:
    print(f"❌ Embedding extraction failed: {e}")
    traceback.print_exc()
    torch.cuda.empty_cache()
    gc.collect()

print("\n" + "=" * 60)
print("✅ SCRIPT COMPLETE")
print(f"📁 Results saved to: {output_dir}")
print("=" * 60)



