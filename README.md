# CATFISHR

**CATFISHR**: CNV And Transcriptome Framework for Identification of Similarity to High-confidence References

CATFISHR identifies malignant cells by integrating RNA expression similarity and inferred CNV similarity to high-confidence malignant reference clusters.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("pshukla63/CATFISHR",
  build_vignettes = TRUE,
  dependencies = TRUE)
```

## Overview

CATFISHR identifies malignant cells by comparing each cell to user-selected high-confidence malignant reference clusters in both RNA and CNV space.

The workflow has three main steps:

1. Calculate Mahalanobis distances from each cell to malignant reference clusters in RNA and CNV PCA space.
2. Combine RNA and CNV distance outputs.
3. Run mean shift clustering in the integrated RNA-CNV distance space.

Users are responsible for:

- Generating PCA embeddings in RNA and CNV space
- Clustering cells in each space
- Selecting high-confidence malignant reference clusters
- Interpreting mean shift cluster output

## Input requirements

CATFISHR requires three inputs for each data space, RNA and CNV:

1. **PCA matrix**  
   A cells-by-PCs matrix. For RNA, this can be generated from scRNA-seq expression data. For CNV, this can be generated from residual expression obtained from inferred CNV methods such as inferCNV.

2. **Cluster membership**  
   A vector of cluster assignments for each cell.

3. **Reference clusters**  
   One or more high-confidence malignant clusters. All query cells are compared against these reference clusters.

For example:

```r
library(CATFISHR)

# Example: extract from Seurat objects
rna_pca <- Seurat::Embeddings(rna_seurat, "pca")
cnv_pca <- Seurat::Embeddings(cnv_seurat, "pca")

rna_clusters <- setNames(rna_seurat$seurat_clusters, colnames(rna_seurat))
cnv_clusters <- setNames(cnv_seurat$seurat_clusters, colnames(cnv_seurat))

# User-defined malignant reference clusters
rna_ref_clusters <- c("0", "2", "5")
cnv_ref_clusters <- c("1", "3")
```

## Toy dataset

CATFISHR includes a small toy dataset for demonstrating the workflow.

```r
library(CATFISHR)

data("RNA_catfishr_data", package = "CATFISHR")
data("CNV_catfishr_data", package = "CATFISHR")

names(RNA_catfishr_data)
names(CNV_catfishr_data)
```

Each toy dataset contains the required inputs:

```r
RNA_catfishr_data$pca_matrix
RNA_catfishr_data$clusters
RNA_catfishr_data$ref_clusters

CNV_catfishr_data$pca_matrix
CNV_catfishr_data$clusters
CNV_catfishr_data$ref_clusters
```

## Usage

### Step 1: Calculate Mahalanobis distances

Calculate distances separately in RNA and CNV space.

```r
mahal_RNA <- calc_mahalanobis(
  pca_matrix = RNA_catfishr_data$pca_matrix,
  clusters = RNA_catfishr_data$clusters,
  ref_clusters = RNA_catfishr_data$ref_clusters,
  n_pcs = 3
)

mahal_CNV <- calc_mahalanobis(
  pca_matrix = CNV_catfishr_data$pca_matrix,
  clusters = CNV_catfishr_data$clusters,
  ref_clusters = CNV_catfishr_data$ref_clusters,
  n_pcs = 3
)
```

### Step 2: Combine RNA and CNV outputs

```r
cm_output <- format_mahal_output(
  mahal_RNA = mahal_RNA,
  mahal_CNV = mahal_CNV
)
```

### Step 3: Run mean shift clustering

```r
ms_output <- run_mean_shift(
  mahal_df = cm_output,
  bandwidths = c(0.4, 0.5, 0.6, 0.7, 0.8),
  max_clusters = 7,
  iterations = 500,
  sample_col = "sample_barcode",
  rna_dist_col = "Mahal_Dist_RNA",
  cnv_dist_col = "Mahal_Dist_CNV"
)
```

### Step 4: Plot mean shift assignments

```r
# Extracting data for plotting
ms_plotting <- ms_output$data

ggplot(ms_plotting, 
  aes(x = log2(Mahal_Dist_CNV), y = log2(Mahal_Dist_RNA))) +
  geom_point(aes(color = factor(Assignment))) +
  labs(color = "Mean Shift\nAssignment")
```

Mean shift clusters with low Mahalanobis distance to the malignant reference clusters in both RNA and CNV space are interpreted as malignant candidates. The final malignant/non-malignant call is made by the user based on the mean shift output.


```r
vignette("CATFISHR", package = "CATFISHR")
```


## Citation

Shukla et al., in preparation.
