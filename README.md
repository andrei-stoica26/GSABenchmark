# GSABenchmark
Gene Set Analysis Benchmark (GSABenchmark) evaluates the ability of gene set analysis methods to identify cells based on ground truth annotation labels.

## Installation

To install GSABenchmark, run the following R code:

```BiocManager::install('andrei-stoica26/GSABenchmark')```

## Contents

The package provides two types of correctness metrics:

1. **Class boundary metrics (8)** : `sensitivity`, `specificity`, `precision`, `accuracy`, `size proximity`, `score specificity`, `silhouette coverage` and `centrality`.
2. **Distribution metrics (4)**: `AUROC`, `Gini coefficient`, `Kolmogorov-Smirnov statistic` and `PRAUC`.

Class boundary metrics provide values for each cell in the input dataset. These signify the metric results obtained if the class identification cutoff were set at the score recorded for the cell by the cell set analysis method. Subsequently, the gene set analysis method score of the cell with the highest average over all class boundary metrics is identified as the class boundary score, and this average is used for comparisons with other methods.

Distribution metrics directly provide single-value scores that can be compared with the scores of other gene set analysis methods.

In addition, GSABenchmark uses [peakRAM](https://cran.r-project.org/web//packages//peakRAM/index.html) to assess the efficiency of the gene set analysis methods in terms of both running time and peak memory usage.

## Methods

Currently, GSABenchmark supports 15 gene set analysis methods:

- `AddModuleScore`
- `AUCell`
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

