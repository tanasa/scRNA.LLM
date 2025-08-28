#!/usr/bin/env python
# coding: utf-8

# In[1]:


# === Standard Library ===
import os
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
from geneformer import EmbExtractor, Classifier
import pickle

# Suppress warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


# In[2]:


# ==============================
# Environment setup
# ==============================
import os
import sys
import torch
import gc
import warnings 

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

torch.cuda.empty_cache()
print(f"CUDA memory after empty_cache: {torch.cuda.memory_allocated() / (1024**3):.2f} GB")

print("✅ Environment variable set!")
print(f"PYTORCH_CUDA_ALLOC_CONF: {os.environ.get('PYTORCH_CUDA_ALLOC_CONF')}")
os.environ['PYTORCH_CUDA_ALLOC_CONF'] = 'expandable_segments:True'

print("✅ Environment variable set!")
print(f"PYTORCH_CUDA_ALLOC_CONF: {os.environ.get('PYTORCH_CUDA_ALLOC_CONF')}")


# In[3]:


# the model was saved in the location :
# /mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_out/250821141257/

# model_version1="V1"
model_version2="V2"

token_dictionary_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl"  # or your V1 or V2 token dict
token_dict_path = token_dictionary_file

with open(token_dictionary_file, "rb") as f:
    token_dict = pickle.load(f)

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


# In[4]:


# =============================================================================
# FOLDERS and INFORMATION
# =============================================================================

model_directory = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_out/250821141257/250821_geneformer_cellClassifier_hypoxia_condition_classifier_optimized/ksplit1/checkpoint-1400/"      # pretrained or fine-tuned model folder
input_data_file = "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_oldlabel.dataset"

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
print(f"\n�� Loading dataset...")
try:
    dataset = load_from_disk(input_data_file)
    print(f"✅ Dataset loaded successfully!")
except Exception as e:
    print(f"❌ Failed to load dataset: {e}")
    exit()


# In[5]:


# =============================================================================
# BASIC INFORMATION
# =============================================================================

print(f"\n📊 BASIC DATASET INFORMATION:")
print(f"   Total cells: {len(dataset):,}")
print(f"   Dataset type: {type(dataset)}")
print(f"   Number of columns: {len(dataset.column_names)}")

print(f"\n📋 METADATA :")
print(f"   Available metadata : {dataset.column_names}")

print(f"\n�� FEATURE SCHEMA:")
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

# Data composition :

import pandas as pd

# Convert to DataFrame (if small enough)
df = dataset.to_pandas()

# Unique medical conditions
print("🧍 Unique medical conditions:", df['condition'].nunique())
print(df['condition'].value_counts(), end="\n\n")

# Unique cell types
print("🔬 Unique cell types:", df['class'].nunique())
print(df['class'].value_counts(), end="\n\n")

# Broad cell types
print("🩺 Broad classes of cell populations:", df['broad_class'].nunique())
print(df['broad_class'].value_counts(), end="\n\n")


# In[6]:


# =============================================================================
# CONFIGURATION
# =============================================================================

CONFIG = {

    # Dataset paths
    "input_data_file": "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_oldlabel.dataset",

    # Model paths
    "model_path": "/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_out/250821141257/250821_geneformer_cellClassifier_hypoxia_condition_classifier_optimized/ksplit1/checkpoint-1400/",
    "token_dictionary_file": "/mnt/nfs/CX000008_DS1/projects/btanasa/Geneformer/geneformer/token_dictionary_gc104M.pkl",

    "target_size": 10000,  # Set to None to use full dataset
    "random_state": 73,

    # Model configuration
    "model_version": "V2",
    "freeze_layers": 2,
    "forward_batch_size": 200,
    "nproc": 10, #

    # Output
    "output_prefix": "hypoxia_2000cells_OPC_ISP"
}


# In[7]:


# =============================================================================
# SETUP AND VERIFICATION
# =============================================================================

print("🚀 Starting Hypoxia Geneformer ISP Pipeline")
print("=" * 60)

# Create timestamped output directory
current_date = datetime.datetime.now()
datestamp = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}{current_date.hour:02d}{current_date.minute:02d}{current_date.second:02d}"
datestamp_min = f"{str(current_date.year)[-2:]}{current_date.month:02d}{current_date.day:02d}"

output_dir = f"/mnt/nfs/CX000008_DS1/projects/btanasa/hypoxia/test_geneformer_5_full_h5ad_aug20_subsample2000cells_dataset_oldlabel.isp/{datestamp}"
os.makedirs(output_dir, exist_ok=True)

print(f"📁 Output directory: {output_dir}")
print(f"📅 Timestamp: {datestamp}")

# Verify all required files exist
print(f"\n🔍 Verifying files...")
required_files = [
    ("Input dataset", CONFIG["input_data_file"]),
    ("Model directory", CONFIG["model_path"]),
    ("Token dictionary", CONFIG["token_dictionary_file"])
]

for name, path in required_files:
    if os.path.exists(path):
        print(f"✅ {name}: {path}")
    else:
        print(f"❌ {name}: {path}")
        print(f"   File not found - please check the path")
        # exit()


# In[8]:


# === condition ===
# Unique (3): ['2mo_HIE', '5d_HIE', 'non-HIE']

# Counts: condition
# 5d_HIE     60633
# non-HIE    25118
# 2mo_HIE    18687

print(output_dir)

# 🔬 Unique cell types: 8
# class
# Astroglia                615
# Imm-MGE                  272
# Inh-Injured              241
# Medium-Spiny-Neuron_?    207
# activated-NSCs           178
# Imm-CGE                  165
# Myelinating-Olig         164
# OPC                      158


# In[9]:


import time

try:
    from geneformer import EmbExtractor

    # Define perturbation parameters
    cell_states_to_model = {
        "state_key": "condition",
        "start_state": "5d_HIE",
        "goal_state": "non-HIE",
        "alt_states": ["2mo_HIE"]
    }

    # Use all disease states (no filtering)
    
    # filter_data_dict = None
    filter_data_dict = {"class":["OPC"]}
    # filter_data_dict={"class":["Astroglia","Imm-MGE","Imm-CGE"]} 

    print(f"Perturbation goal: 5d_HIE → non-HIE (with 2mo_HIE as alternative)")
    print(f"Cell states: {cell_states_to_model}")
 
    # Setup embedding extractor
    print(f"🔧 Initializing embedding extractor...")
    embex = EmbExtractor(
        model_type = "CellClassifier", 
        num_classes = 3,
        filter_data = filter_data_dict,
        max_ncells = 1000,  # maximum number of cells to use 
        emb_layer = 0,
        summary_stat = "exact_mean",
        forward_batch_size = 16, # 256 in the tutorial
        model_version = "V2",
        nproc = 20, # 1
        emb_mode = "cls"             # <-- ADD for V2
        # cell_emb_style="mean_pool" # <-- REMOVE in CLS mode
    )

     # Get state embeddings with timing
    print(f"📊 Extracting state embeddings...")
    start_time = time.time()
    
    # Get state embeddings
    print(f"📊 Extracting state embeddings...")
    state_embs_dict = embex.get_state_embs(
        cell_states_to_model,
        CONFIG["model_path"],
        CONFIG["input_data_file"],
        output_dir,
        CONFIG["output_prefix"]
    )
    print(f"✅ State embeddings extracted")

    end_time = time.time()
    execution_time = end_time - start_time
    
    print(f"✅ State embeddings extracted")
    print(f"⏱️ Execution time: {execution_time:.2f} seconds ({execution_time/60:.2f} minutes)")

except Exception as e:
    print(f"❌ State embedding extraction failed: {e}")
    print(f"   Cannot proceed without state embeddings")
    import traceback
    traceback.print_exc()
    # exit()


# In[ ]:





# In[10]:


print(output_dir)
print(CONFIG["output_prefix"])
print(CONFIG["input_data_file"])


# In[ ]:





# In[11]:


import pickle
import numpy as np

# Reconstruct the pickle filename
pkl_file = os.path.join(output_dir, f"{CONFIG['output_prefix']}.pkl")

print(f"📁 Output directory: {output_dir}")
print(f"🏷️ Output prefix: {CONFIG['output_prefix']}")
print(f"📄 Reconstructed pickle file: {pkl_file}")

# Check if the file exists
if os.path.exists(pkl_file):
    print(f"✅ File exists: {pkl_file}")
else:
    print(f"❌ File does NOT exist: {pkl_file}")


# Load the file
with open(pkl_file, 'rb') as f:
    data = pickle.load(f)

# Show what's inside
print("Keys:", list(data.keys()))
print("Number of states:", len(data))

# Show details for each state
for state, emb in data.items():
    print(f"\n{state}:")
    print(f"  Type: {type(emb)}")
    
    # Convert to numpy if it's a tensor
    if hasattr(emb, 'detach'):
        emb_array = emb.detach().cpu().numpy()
    else:
        emb_array = np.array(emb)
    
    print(f"  Shape: {emb_array.shape}")
    print(f"  First 5 values: {emb_array.flatten()[:10]}")

# Show size of each vector
for state, emb in data.items():
    print(f"{state}: {emb.shape}")

# So when you see:
emb_array.flatten()[:5]


# In[12]:


# 


# In[13]:


print("The analysis on the state embeddings:")

# STEP 1 — Imports

import os, json, time
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# STEP 2 — Prepare paths/names

prefix = f"{CONFIG['output_prefix']}_states"
emb_dir = Path(output_dir) / "state_embeddings"
emb_dir.mkdir(parents=True, exist_ok=True)

# Normalize embeddings to CPU NumPy (2D arrays)
import torch
import numpy as np
import pandas as pd
from pathlib import Path

prefix = f"{CONFIG['output_prefix']}_states"
emb_dir = Path(output_dir) / "state_embeddings"
emb_dir.mkdir(parents=True, exist_ok=True)

state_np = {}
for state, arr in state_embs_dict.items():
    # Convert possible Torch tensors (GPU or CPU) to NumPy
    if isinstance(arr, torch.Tensor):
        A = arr.detach().cpu().numpy()
    elif isinstance(arr, (list, tuple)) and len(arr) > 0 and isinstance(arr[0], torch.Tensor):
        A = torch.stack(arr).detach().cpu().numpy()
    else:
        A = np.asarray(arr)

    # Ensure 2D shape
    if A.ndim == 1:
        A = A.reshape(1, -1)

    state_np[state] = A

# STEP 3 — Save each state's embeddings to .npy and .csv 

for state, A in state_np.items():
    np.save(emb_dir / f"{prefix}_{state}.npy", A)
    pd.DataFrame(A).to_csv(emb_dir / f"{prefix}_{state}.csv", index=False)
print("Saved per-state embeddings.")

# STEP 4 — Build combined matrix + labels; save
X_list, labels = [], []
for state, A in state_np.items():
    X_list.append(A)
    labels += [state] * A.shape[0]

X = np.vstack(X_list)
combined = pd.DataFrame(X)
combined["state"] = labels
np.save(emb_dir / f"{prefix}_all.npy", X)
combined.to_csv(emb_dir / f"{prefix}_all.csv", index=False)
print("Saved combined embeddings.")

# STEP 5 — Cosine-similarity heatmap of STATE MEANS
import matplotlib.pyplot as plt

states = list(state_np.keys())
means = np.stack([A.mean(axis=0) for A in state_np.values()], axis=0)

# cosine similarity without external deps
norms = np.linalg.norm(means, axis=1, keepdims=True) + 1e-12
C = (means / norms) @ (means / norms).T

plt.figure(figsize=(4, 3))
im = plt.imshow(C, vmin=0, vmax=1)
plt.xticks(range(len(states)), states, rotation=45, ha="right")
plt.yticks(range(len(states)), states)
plt.title("Cosine similarity of state means")
plt.colorbar(im, fraction=0.046, pad=0.04)
plt.tight_layout()
plt.savefig(emb_dir / f"{prefix}_cosine_heatmap.png", dpi=200)
plt.close()
print("Saved cosine-similarity heatmap.")

# STEP 6 — PCA scatter 

from IPython.display import display
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if X.shape[0] >= 3:
    pca = PCA(n_components=2, random_state=73)
    Xp = pca.fit_transform(X)

    # Plot (SHOW + save)
    plt.figure(figsize=(5, 4))
    for s in sorted(set(labels)):
        m = (np.array(labels) == s)
        plt.scatter(Xp[m, 0], Xp[m, 1], s=8, alpha=0.7, label=s)
    plt.legend()
    plt.title("PCA of embeddings")
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% var)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% var)")
    plt.tight_layout()
    plt.show()  # ⬅️ display inline
    plt.savefig(emb_dir / f"{prefix}_pca.png", dpi=200)
    plt.close()

    # Save coordinates + preview
    pca_df = pd.DataFrame({"PC1": Xp[:, 0], "PC2": Xp[:, 1], "state": labels})
    pca_df.to_csv(emb_dir / f"{prefix}_pca.csv", index=False)
    np.save(emb_dir / f"{prefix}_pca.npy", Xp)
    print("PCA coords saved. Preview:")
    display(pca_df.head())
else:
    print("Skipped PCA: need ≥3 rows in X.")


# In[ ]:





# In[14]:


# create the output files 

try:
    with open(CONFIG["token_dictionary_file"], "rb") as f:
        tok = pickle.load(f)
    keys = set(tok.keys()) if isinstance(tok, dict) else set()
    vals = set(tok.values()) if isinstance(tok, dict) else set()
    missing = [g for g in genes_to_perturb if g not in keys and g not in vals]
    if missing:
        print("⚠️ WARNING: not found in token dictionary:", missing)
except Exception as e:
    print("ℹ️ Skipping token-dict validation:", e)

# Output folders
isp_output_dir = os.path.join(output_dir, "isp")
isp_stats_dir  = os.path.join(output_dir, "isp_stats")
os.makedirs(isp_output_dir, exist_ok=True)
os.makedirs(isp_stats_dir, exist_ok=True)


# In[15]:


# === Geneformer ===
from geneformer import EmbExtractor, Classifier, InSilicoPerturber, InSilicoPerturberStats
import pickle


# In[16]:


print(output_dir)
print(isp_output_dir)
print(isp_stats_dir)


# In[ ]:


genes_to_perturb = "all"

import time

# Time the ISP initialization
print("�� Initializing InSilicoPerturber...")
init_start = time.time()

isp = InSilicoPerturber(
    perturb_type = "delete",
    perturb_rank_shift = None,
    genes_to_perturb = genes_to_perturb,
    # genes_to_perturb = "all",
    combos = 0,
    anchor_gene = None,
    model_type = "CellClassifier",
    num_classes = 3,
    emb_mode = "cls",                     # ✅ Use CLS embedding for V2
    # emb_mode = "cell",                    # SET TO "CELL" FOR V1 MODEL. FOR V2, SHOULD BE "CLS" (current default).
    # emb_mode = "cls_and_gene",
    filter_data = filter_data_dict,
    cell_states_to_model = cell_states_to_model,
    state_embs_dict = state_embs_dict,
    max_ncells = 1200,
    emb_layer = 0,
    forward_batch_size = 32,
    model_version = "V2",                 # ✅ Matches CLS usage
    nproc = 130
)

# --- avoid the Column * int bug by skipping the special-token code path ---
if hasattr(isp, "special_token"):
    isp.special_token = False

# Time the perturbation execution
print("🧬 Running ISP perturbation...")
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

# Goal-state-shift statistics (Chronic → Normoxia, Acute as alt)
print("�� Initializing InSilicoPerturberStats...")
stats_init_start = time.time()

ispstats = InSilicoPerturberStats(
    mode="goal_state_shift",
    genes_perturbed = genes_to_perturb,
    combos = 0,
    anchor_gene = None,
    cell_states_to_model = cell_states_to_model,
    model_version = "V2",
    token_dictionary_file = token_dictionary_file
)

stats_init_time = time.time() - stats_init_start
total_time_with_stats = time.time() - init_start

print(f"✅ InSilicoPerturberStats initialized in {stats_init_time:.2f} seconds")
print(f"⏱️ Total time (init + execution + stats_init): {total_time_with_stats:.2f} seconds ({total_time_with_stats/60:.2f} minutes)")

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

print(f"✅ Statistics computed -> {isp_stats_dir}")
print(f"⏱️ Statistics computation time: {stats_comp_time:.2f} seconds ({stats_comp_time/60:.2f} minutes)")


# In[ ]:


# After ISP execution, automatically load results from isp_output_dir
print(f"📁 Loading ISP results from: {isp_output_dir}")

# Look for pickle files in the isp_output_dir
import glob
pickle_files = glob.glob(os.path.join(isp_output_dir, "*_raw.pickle"))

if pickle_files:
    # Use the most recent file (in case there are multiple)
    latest_pickle = max(pickle_files, key=os.path.getctime)
    print(f"✅ Found ISP results: {os.path.basename(latest_pickle)}")
    
    # Load automatically
    with open(latest_pickle, 'rb') as f:
        data = pickle.load(f)
    
    print(f"\n📊 File contents:")
    print(f"Type: {type(data)}")
    
    if isinstance(data, dict):
        print(f"Keys: {list(data.keys())}")
        
        # Explore the structure
        for key, value in data.items():
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
        print(f"Content: {data}")
    
    # Store in perturbation_results for later use
    perturbation_results = data
    print(f"\n🎯 Results loaded into 'perturbation_results' variable!")
    
else:
    print(f"❌ No ISP results found in {isp_output_dir}")
    print(f"Available files:")
    for f in sorted(os.listdir(isp_output_dir)):
        print(f"  {f}")
    perturbation_results = None


# In[ ]:





# In[ ]:




