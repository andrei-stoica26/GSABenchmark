# GSABenchmark
GSABenchmark is a package designed for benchmarking scRNA-seq gene set analysis 
(scGSA) methods. It provides both traditional and novel benchmark metrics, 
as well as visualization tools. Currently, GSABenchmark supports 17 scGSA 
methods.

## Installation

To install GSABenchmark, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install('GSABenchmark')
```

## Benchmark metrics

Three types of correctness metrics are implemented in GSABenchmark:

1. **Class boundary determination metrics (6)** : `sensitivity`, `specificity`, `precision`, `accuracy`, `size proximity` and `score coverage`.
2. **Matthews correlation coefficient-based metrics (2)**: `Threshold-based MCC` and `Comprehensive MCC`.
3. **Global evaluation metrics (7)**: `AUROC`, `PRAUC`, `label rank alignment`, `silhouette rank alignment`, `centrality`, `label Jaccard score` and  `label cosine score`.

A brief explanation of these metric classes is as follows:
- Class boundary metrics are evaluated at all distinct scores obtained by a 
scGSA method. They represent values obtained if the class identification 
cutoff were set at the corresponding score. Subsequently, the results of the 
six metrics are averaged and the scGSA method score corresponding to the 
highest average is taken to determine the class boundary cutoff. 
- Global evaluation metrics provide threshold-independent assessments of 
prediction quality. 
- As the gold standard in binary classification, the Matthews correlation 
coefficient (MCC) is used separately both as a class boundary metric and as a 
global evaluation metric.

Additionally, GSABenchmark uses [peakRAM](https://cran.r-project.org/web//packages//peakRAM/index.html) to assess the efficiency of the gene set analysis methods in terms of both running time and peak memory usage.

## Methods

Currently, GSABenchmark supports 17 scGSA methods: 

- `AddModuleScore` (reimplemented from [Seurat](https://github.com/satijalab))
- `AUCell`
- `CSOA`
- `GSVA`
- `JASMINE`(reimplemented from [JASMINE](https://github.com/NNoureen/JASMINE))
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

All methods can run on both `Seurat` and `SingleCellExperiment` objects.

## Assessing new methods

Users can also define their own gene set analysis methods to be assessed with 
GSABenchmark. User-defined methods must meet the following criteria:

- The name of the method starts with `run`.
- The first parameter is a `Seurat` or `SingleCellExperiment` object.
- The second parameter is a list of character vectors representing gene sets.
- The method returns an object of the same type as the first parameter, having
cell scores in the metadata column denoted by the third parameter. 

Additional parameters are allowed. For instance, this function heading would 
be recognized by GSABenchmark: 

`runCustomMethod(seuratObj, geneSets, colStr, tol = 0.2, verbose = TRUE)`. 

**Note**: When passing the method name to a GSABenchmark function as part of 
the `gsaMethods` argument, the method needs to be referred to as `CustomMethod` 
(no `run`).

GSABenchmark min-max-normalizes the per-cell gene set analysis scores of all 
supported methods, returning scores between 0 and 1. It is recommended that 
user-defined methods return scores in the same range.
