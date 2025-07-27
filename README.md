# GSABenchmark
Gene Set Analysis Benchmark (GSABenchmark) evaluates the ability of gene set analysis methods to identify cells based on ground truth annotation labels.

## Installation

To install GSABenchmark, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install('andrei-stoica26/GSABenchmark')
```

## Contents

Three types of correctness metrics are implemented in GSABenchmark:

1. **Class boundary determination metrics (6)** : `sensitivity`, `specificity`, `precision`, `accuracy`, `size proximity` and `score coverage`.
2. **Matthews correlation coefficient-based metrics (2)**: `Threshold-based MCC`, `Comprehensive MCC`
3. **Global evaluation metrics (5)**: `AUROC`, `PRAUC`, `label rank alignment`, `silhouette rank alignment` and `centrality`.

Class boundary metrics provide values for each cell in the input dataset, representing values obtained if the class identification cutoff were set at the gene set analysis score recorded by the corresponding cell. Subsequently, the results of the six metrics are averaged and the gene set analysis method score corresponding to the highest average is taken to determine the class boundary cutoff. Global evaluation metrics provide threshold-independent assessments of prediction quality. As the gold standard in binary classification, the Matthews correlation coefficient is used separetely both as a class boundary metric and as a global evaluation metrics.

In addition, GSABenchmark uses [peakRAM](https://cran.r-project.org/web//packages//peakRAM/index.html) to assess the efficiency of the gene set analysis methods in terms of both running time and peak memory usage.

## Methods

Currently, GSABenchmark supports 15 gene set analysis methods for both `Seurat` 
and `SingleCellExperiment` objects:

- `AUCell`
- `CSOA`
- `GSVA`
- `MDT`
- `MLM`
- `ORA`
- `Pagoda2`
- `PLAGE`
- `Singscore`
- `SiPSiC`
- `ssGSEA`
- `UCell`
- `UDT`
- `VAM`
- `Zscore`

Users can also define their own gene set analysis methods to be assessed with GSABenchmark. User-defined methods must meet the following criteria:

- The name of the method starts with `run`.
- The first parameter is a Seurat object.
- The second parameter is a character vector of genes.
- The third parameter is a character.
- The method returns a Seurat object with cell scores in the metadata column denoted by the third parameter. 

Additional parameters are allowed. For instance, this function heading would be recognized by GSABenchmark's relevant functions, `runGSAMethods` and `efficiencyBenchmark`: 

`runCustomMethod(seuratObj, genes, colStr, tol = 0.2, verbose = TRUE)`. 

**! Note** : When calling `runGSAMethods` or `efficiencyBenchmark`, the method needs to be referred to as `CustomMethod` (no `run`).

GSABenchmark min-max-normalizes the per-cell gene set analysis scores of all supported methods, returning scores between 0 and 1. It is recommended that user-defined methods follow the same convention.
