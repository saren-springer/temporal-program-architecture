# WATER: Windowed Assignment of Temporal Expression of RNA
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE.txt)

WATER reconstructs temporal expression trajectories independently within each biological condition, audits assignment confidence, refines weak or ambiguous assignments, and compares trajectory organization across conditions.

## Table of Contents
- [Source Publication](#source-publication)
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Requirements](#requirements)
- [Input Files](#input-files)
- [Running WATER](#running-water)
- [Worked Example](#worked-example)
- [Comparing Conditions](#comparing-conditions)
- [Output Structure](#output-structure)
- [Citation](#citation)

## Source Publication
### ***WATER: A Framework for Reconstructing Temporal Gene Expression Trajectories and Detecting Molecular Heterochrony Across Development and Disease***

*Saren M. Springer<sup>1,6</sup>, Kevon O. Afriyie<sup>1</sup>, Abigail R Boria<sup>1</sup>, Isabella Stevens<sup>1</sup>, Katarina L Micin<sup>2,3</sup>, Vanessa Aguiar-Pulido<sup>3,4,5</sup>, Rahul N Kanadia<sup>1,6,\*</sup>*

*<sup>1</sup>Department of Physiology and Neurobiology, University of Connecticut, Storrs, CT 06269, USA  
<sup>2</sup>Department of Human Genetics and John P. Hussman Institute of Human Genomics, University of Miami Miller School of Medicine, Miami FL, USA  
<sup>3</sup>Sylvester Comprehensive Cancer Center at the University of Miami Miller School of Medicine, Miami, FL  
<sup>4</sup>Department of Computer Science, University of Miami, Coral Gables, FL  
<sup>5</sup>Department of Informatics and Health Data Science, University of Miami Miller School of Medicine, Miami, FL  
<sup>6</sup>Institute for Systems Genomics, University of Connecticut, Storrs, CT 06269, USA*

 
<sup>*</sup>To whom correspondence should be addressed.

<br>

## Overview

WATER is an R-based workflow for reconstructing and comparing temporal gene-expression trajectories. The core pipeline:

1. filters features by expression,
2. classifies features as differentially abundant or non-differentially abundant using supplied log2 fold-change values,
3. quantile-normalizes replicate-level expression data,
4. averages normalized replicates by time point,
5. performs hierarchical clustering with silhouette-guided selection of the cluster number,
6. collapses clusters into temporal trajectories,
7. audits assignments using correlation and leave-one-out stability analyses,
8. reclusters ambiguous and weakly assigned features, and
9. exports final trajectories, confidence classifications, plots, and analysis-ready tables.

`WATER_condition_comparison.R` compares two completed WATER runs. Condition 1 is treated as the baseline or control condition. The comparison workflow identifies shared and condition-specific trajectories, maps trajectory transitions, generates transition plots, and optionally performs Gene Ontology enrichment.

## Repository Structure

```plaintext
temporal-program-architecture/
|-- README.md
|-- LICENSE.txt
|-- vignettes/
|   `-- water_5xFAD_female_centroid_vignette.Rmd
|-- docs/
|   `-- index.html
|-- WATER/
  |-- Running_WATER.R
  |-- WATER_condition_comparison.R
  |-- filter.txt  
  |-- test_data
  |   `-- 5xFAD_female_Data
  |       |-- WT_female_cortex
  |       |  |-- WT_female_cortex_exp_data.csv 
  |       |  `-- WT_female_cortex_log2fc_data.csv  
  |       `-- 5xFAD_female_cortex
  |          |-- 5xFAD_female_cortex_exp_data.csv 
  |          `-- 5xFAD_female_cortex_log2fc_data.csv    
  |-- preprocess/
  |   |-- water_preprocess.R
  |   |-- water_initialize.R
  |   |-- water_filter_expression.R
  |   |-- water_filter_features.R
  |   |-- water_normalize.R
  |   |-- water_PCA.R
  |   |-- water_get_output_path.R
  |   `-- water_save_output.R
  |-- clustering/
  |   |-- water_cluster_classic.R
  |   |-- water_cluster_initial.R
  |   |-- water_cut.R
  |   `-- water_label.R
  |-- audit/
  |   |-- water_audit_classic.R
  |   |-- water_initial_LOO.R
  |   `-- water_correlate_initial.R
  `-- refine/
      |-- water_refine_classic.R
      |-- water_recluster_ambiguous.R
      |-- water_recluster_weak.R
      |-- water_correlate_secondary.R
      `-- water_secondary_LOO.R
```

### Main scripts

- `Running_WATER.R` contains `run_water()`, the primary entry point for a single condition.
- `WATER_condition_comparison.R` contains `run_condition_comparison()`, which compares two completed WATER analyses.

Because the scripts use relative `source()` calls, run WATER from the repository root unless the paths are modified.

## Requirements

WATER requires R and the following packages.

### CRAN packages

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "readr",
  "readxl",
  "tibble",
  "ggplot2",
  "patchwork",
  "pheatmap",
  "factoextra",
  "mclust",
  "clue",
  "magrittr",
  "openxlsx",
  "svglite"
))
```

### Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "preprocessCore",
  "clusterProfiler",
  "org.Mm.eg.db"
))
```

`clusterProfiler` and `org.Mm.eg.db` are only required for Gene Ontology enrichment during condition comparison. The current implementation assumes mouse gene symbols. For another organism, replace `org.Mm.eg.db` and update the identifier mapping accordingly.

## Input Files

### Expression matrix

`data` may be a data frame or matrix.

- For a data frame, the first column must contain feature identifiers and the remaining columns must contain samples.
- For a matrix, feature identifiers must be stored as row names.
- Sample columns must be ordered by timepoint and replicate in the same order supplied through `replicates`.

Example:

```plaintext
Gene_ID,T1_rep1,T1_rep2,T1_rep3,T2_rep1,T2_rep2,T2_rep3
GeneA,10.2,11.1,9.8,14.3,13.9,15.0
GeneB,2.1,2.4,2.0,1.8,1.9,1.7
```

### Replicate structure

The replicate structure is supplied as a named list:

```r
replicates = list(T1 = 3, T2 = 3, T3 = 3, T4 = 3)
```

Unequal replicate numbers are supported as long as the sample columns are ordered consistently.

### Log2 fold-change matrix

`log2fc` is required for differential-abundance classification.

- The first column must contain feature identifiers, or identifiers must be stored as row names.
- Remaining columns contain pairwise log2 fold-change values.
- A feature is classified as differentially abundant when at least one absolute log2 fold-change meets `log2fc_threshold`.

### Optional feature list

`feature_list` may point to a plain-text file containing one feature identifier per line. When provided, WATER restricts the analysis to matching features before trajectory reconstruction.

## Running WATER

Download or clone this repository, then set the R working directory to the local `WATER` folder.

```r
# Set working directory to the WATER folder
setwd("/WATER")

# Read input data
data <- read.csv("path/to/expression_data.csv")
log2fc <- read.csv("path/to/log2fc_data.csv")

# Load WATER
source("Running_WATER.R")

# Run one condition
water <- run_water(
  data = data,
  replicates = list(T1 = 3, T2 = 3, T3 = 3, T4 = 3),
  log2fc = log2fc,
  preprocess_params = list(
    modality = "RNA",
    output_dir = "WT_water_output",
    save_outputs = TRUE,
    tpm_threshold = 1,
    feature_list = "filter.txt",
    log2fc_threshold = 1
  ),
  r_threshold = 0.8,
  run_refine = TRUE,
  run_loo = TRUE
)
```

During clustering, WATER reports recommended values of `k` and prompts the user to select the final number of clusters. It then prompts the user to choose a trajectory-labeling approach:

1. centroid-based trajectory collapsing,
2. manual cluster-to-trajectory mapping, or
3. retention of each cluster as a separate trajectory.

Run the pipeline independently for every condition that will be compared.

## Worked Example

A complete worked example is provided using a female wildtype 5xFAD aging time course with four ages and three biological replicates per age.

The vignette demonstrates:

- input formatting,
- preprocessing and differential-abundance filtering,
- silhouette-guided cluster selection,
- centroid-based trajectory collapse,
- trajectory renaming,
- correlation-based confidence auditing,
- leave-one-out stability analysis,
- reclustering of ambiguous and weak features, and
- interpretation of the final WATER outputs.

See the [worked WATER vignette](https://saren-springer.github.io/temporal-program-architecture/).

In this example, WATER was run using the default centroid-based collapse method with a correlation threshold of 0.8.

## Comparing Conditions

First complete `run_water()` for both conditions. Then source the comparison module:

```r
source("WATER_condition_comparison.R")
```

Condition 1 is treated as the baseline or control condition.

```r
comparison <- run_condition_comparison(
  cond1_name = "WT",
  cond2_name = "MUT",
  cond1_dir = "WT_water_output",
  cond2_dir = "MUT_water_output",
  out_dir = "WT_MUT_comparison",
  mut_colors = mut_colors,
  run_GO = TRUE
)
```

An optional color map may be supplied for condition-2 trajectories:

```r
mut_colors <- tibble::tribble(
  ~Mut_Traj,    ~spaghetti, ~mean,
  "Early High", "#489744",  "#306f37",
  "Late High",  "#58acbf",  "#31798f",
  "default",    "grey60",   "grey20"
)
```

The comparison module:

- merges feature assignments across conditions,
- labels missing assignments as `NonDA` or `Not_Exp`,
- quantifies trajectory overlap and transition counts,
- produces per-transition gene lists and trajectory plots,
- produces baseline-only plots and dominant-shift plots,
- runs BP, MF, and CC Gene Ontology enrichment when requested, and
- writes combined GO enrichment tables.

## Output Structure

Each WATER run creates step-specific output directories beneath `output_dir`.

```plaintext
WT_water_output/
|-- initialize/
|-- filter_expression/
|-- filter/
|-- normalize/
|-- PCA/
|-- cluster/
|-- cut_classic/
|-- label_classic/
|-- loo_initial/
|-- correlate_initial/
|-- recluster_ambiguous/
|-- recluster_weak/
|-- correlate_secondary/
|-- loo_secondary/
`-- final/
    |-- tables/
    |   `-- water_final_gene_table.csv
    |-- plots/
    |   |-- water_final_trajectories.png
    |   `-- water_final_trajectories.pdf
    `-- objects/
        `-- water_final_results.rds
```

Intermediate folders contain the tables, plots, and objects produced at each stage. Key normalized expression outputs include:

```plaintext
normalize/tables/water_data_normalized_replicates.xlsx
normalize/tables/water_data_normalized_avg.xlsx
```

The condition-comparison output includes:

```plaintext
WT_MUT_comparison/
|-- Trajectory_overlap_between_conditions.csv
|-- Trajectory_overlap_between_conditions.pdf
|-- trajectory_transition_counts.csv
|-- trajectory_comparison_table.csv
|-- ALL_WT_MUT_GO_ENRICHMENT.csv
|-- ALL_WT_ONLY_GO_ENRICHMENT.csv
`-- <baseline_trajectory>/
    |-- WT_only_genes.csv
    |-- WT_only_plot.pdf
    |-- WT_only_dominant_shift_plot.pdf
    |-- WT_only_dominant_shift_genes.csv
    `-- <baseline_to_condition2_trajectory>/
        |-- genes.csv
        |-- trajectory_plot.pdf
        |-- *_BP_GO.csv
        |-- *_MF_GO.csv
        `-- *_CC_GO.csv
```

GO output files are created only when enrichment results are available.

## Citation

If you use WATER, please cite:

**Springer S.M., et al.**  
*WATER: A Framework for Reconstructing Temporal Gene Expression Trajectories and Detecting Molecular Heterochrony Across Development and Disease*  
Manuscript in submission.

**Repository:**  
Springer S.M., et al.  
*WATER: Windowed Assignment of Temporal Expression of RNA.*  
GitHub: [saren-springer/temporal-program-architecture](https://github.com/saren-springer/temporal-program-architecture)

<br><br>
