#' Perform accuracy benchmarks for multiple sets of GSA method scores and class
#' labels
#'
#' This function performs accuracy benchmarks for multiple sets of GSA method scores
#' and class labels.
#'
#' @inheritParams accuracyBenchmark
#' @inheritParams runGSAMethods
#' @param labels The assessed labels. Must be a subset of the values in labelCol.
#' @param benchmarkFun Benchmark function. Must be either boundaryBenchmark or
#' distributionBenchmark.
#' @param verbose Whether to output a message whenever the benchmarking starts on
#' a gene set
#' @param ... Additional paramters passed to benchmarkFun.
#'
#' @return A list of lists of benchmark data frames.
#'
#'
accuracyBenchmarkMultiple <- function(seuratObj, labelCol, geneSetNames, gsaMethods,
                                      labels = geneSetNames, benchmarkFun, verbose = FALSE, ...){
  res <- lapply(seq_along(geneSetNames), function(i) {
    setName <- geneSetNames[i]
    if (verbose)
      message('Running benchmark on gene set: ', geneSetNames[i], '...')
    setRes <- lapply(seq_along(gsaMethods), function(j){
      method <- gsaMethods[j]
      scoreCol <- paste0(method, '_', setName)
      df <- benchmarkFun(seuratObj, labelCol, scoreCol, labels[i], ...)
      return(df)
    })
    names(setRes) <- gsaMethods
    return(setRes)
  })
  names(res) <- labels
  return(res)
}

#' Perform boundary accuracy benchmarks for multiple sets of GSA method scores
#' and class labels
#'
#' This function performs boundary accuracy benchmarks for multiple sets of GSA
#' method scores and class labels.
#'
#' @inheritParams accuracyBenchmarkMultiple
#'
#' @return A list of lists of boundary benchmark data frames
#'
#' @export
#'
boundaryBenchmarkMultiple <- function(seuratObj, labelCol, geneSetNames, gsaMethods,
                                      labels = geneSetNames, verbose = FALSE)
  return(accuracyBenchmarkMultiple(seuratObj, labelCol, geneSetNames, gsaMethods, labels,
                                   boundaryBenchmark))

#' Perform direct MCC benchmarks for multiple sets of GSA method scores
#' and class labels
#'
#' This function performs direct MCC benchmarks for multiple sets of GSA
#' method scores and class labels.
#'
#' @inheritParams accuracyBenchmarkMultiple
#'
#' @return A list of lists of boundary benchmark data frames
#'
#' @export
#'
directMCCBenchmarkMultiple <- function(seuratObj, labelCol, geneSetNames, gsaMethods,
                                      labels = geneSetNames)
  return(accuracyBenchmarkMultiple(seuratObj, labelCol, geneSetNames, gsaMethods, labels,
                                   directMCCBenchmark, verbose = FALSE))

#' Perform distribution accuracy benchmarks for multiple sets of GSA method scores
#' and class labels
#'
#' This function performs distribution accuracy benchmarks for multiple sets of GSA
#' method scores and class labels.
#'
#' @inheritParams accuracyBenchmarkMultiple
#' @inheritParams globalBenchmark
#'
#' @return A list of lists of distribution benchmark data frames
#'
#' @export
#'
globalBenchmarkMultiple <- function(seuratObj, labelCol, geneSetNames, gsaMethods,
                                          labels = geneSetNames, normSilDF = NULL,
                                          dimMat = NULL, maxDist = NULL)
  return(accuracyBenchmarkMultiple(seuratObj, labelCol, geneSetNames, gsaMethods, labels,
                                   globalBenchmark, normSilDF, dimMat, maxDist, verbose = FALSE))
