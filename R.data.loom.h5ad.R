
############################################################################
############################################################################

# It creates loom and h5ad files for Geneformer Input

# Key Requirements for Geneformer Input
# File Format: .loom or .h5ad (both HDF5-based).
# Gene Annotation:
# Ensembl IDs: Genes must be labeled with Ensembl IDs (e.g., ENSG00000123456).

# Location:
# .loom: Stored as a row attribute named "ensembl_id".
# .h5ad: Stored as a column attribute due to transposed formatting.

# Cell Metadata:
# Total Read Counts: Column attribute "n_counts" (required for normalization).

# Optional Filtering: Binary column attribute "filter_pass" for user-defined cell filtering.

# Data Content:
# Raw Counts: Unnormalized, unfiltered UMI counts without feature selection.

# Recommended Approach: Convert Public Datasets
# The scRNAseq R package provides publicly available datasets that can be converted to 
# Geneformer-compatible HDF5 files. Below, we use the Darmanis brain dataset as an example:

############################################################################
############################################################################

# Public Datasets Compatible with Geneformer

# Dataset	Format	Access Method	Notes
# DarmanisBrain	SingleCellExperiment	scRNAseq::DarmanisBrainData()	Human brain scRNA-seq.
# GrunPancreas	SingleCellExperiment	scRNAseq::GrunPancreasData()	Human pancreas scRNA-seq.
# HeOrganAtlas	SingleCellExperiment	scRNAseq::HeOrganAtlasData()	Human cortex snRNA-seq.

############################################################################
############################################################################

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scRNAseq")
BiocManager::install("LoomExperiment")

# Install sceasy if you haven't already
install.packages("remotes")
remotes::install_github("cellgeni/sceasy")
remotes::install_github("cellgeni/reticulate") # for Python interoperability

############################################################################

library(reticulate)
use_python("/home/tanasa/anaconda3/bin/python") # or use_condaenv("your_conda_env")

library(scRNAseq)
library(LoomExperiment)
library(sceasy)

############################################################################

# Load dataset with Ensembl IDs
sce <- DarmanisBrainData(ensembl = TRUE)

# Add total read counts per cell (required by Geneformer)
colData(sce)$n_counts <- colSums(counts(sce))

# Convert to .loom format
scle <- as(sce, "SingleCellLoomExperiment")
export(scle, "darmanis_brain.loom")

# Convert loom to h5ad
sceasy::convertFormat(
  "darmanis_brain.loom",
  from = "loom",
  to = "anndata",
  outFile = "darmanis_brain.h5ad"
)

############################################################################
############################################################################

# Load dataset with Ensembl IDs
sce <- DarmanisBrainData(ensembl = TRUE)

# Add total read counts per cell (required by Geneformer)
colData(sce)$n_counts <- colSums(counts(sce))

# Convert to .loom format
scle <- as(sce, "SingleCellLoomExperiment")
export(scle, "darmanis_brain.loom")

# Convert loom to h5ad
sceasy::convertFormat(
  "darmanis_brain.loom",
  from = "loom",
  to = "anndata",
  outFile = "darmanis_brain.h5ad"
)

############################################################################
############################################################################

# Load dataset with Ensembl IDs
sce <- scRNAseq::GrunPancreasData(ensembl = TRUE)

# Add total read counts per cell (required by Geneformer)
colData(sce)$n_counts <- colSums(counts(sce))

# Convert to .loom format
scle <- as(sce, "SingleCellLoomExperiment")
export(scle, "grun_pancreas.loom")

# Convert loom to h5ad
sceasy::convertFormat(
  "grun_pancreas.loom",
  from = "loom",
  to = "anndata",
  outFile = "grun_pancreas.h5ad"
)

############################################################################
############################################################################

# Load dataset with Ensembl IDs
sce <- scRNAseq::HeOrganAtlasData(ensembl = TRUE)

# Add total read counts per cell (required by Geneformer)
colData(sce)$n_counts <- colSums(counts(sce))

# Convert to .loom format
scle <- as(sce, "SingleCellLoomExperiment")
export(scle, "he_organ_atlas.loom")

# Convert loom to h5ad
sceasy::convertFormat(
  "he_organ_atlas.loom",
  from = "loom",
  to = "anndata",
  outFile = "he_organ_atlas.h5ad"
)

############################################################################
############################################################################
