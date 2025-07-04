#' Perform accuracy benchmarks for multiple sets of GSA method scores and class
#' labels
#'
#' This function performs accuracy benchmarks for multiple sets of GSA method scores
#' and class labels
#'
#' @inheritParams accuracyBenchmark
#' @inheritParams runGSAMethods
#' @param labels The assessed labels. Must be a subset of the values in labelCol.
#' @param benchmarkFun Benchmark function. Must be either boundaryBenchmark or
#' distributionBenchmark
#' @param ... Additional paramters passed to benchmarkFun
#'
#' @return A list of lists of benchmark data frames
#'
#'
accuracyBenchmarkMultiple <- function(seuratObj, labelCol, geneSetNames, gsaMethods,
                                      labels = geneSetNames, benchmarkFun, ...){
  res <- lapply(seq_along(geneSetNames), function(i) {
    setName <- geneSetNames[i]
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
#' @inheritParams boundaryBenchmark
#'
#' @return A list of lists of boundary benchmark data frames
#'
#' @export
#'
boundaryBenchmarkMultiple <- function(seuratObj, labelCol, geneSetNames, gsaMethods,
                                      labels = geneSetNames, normSilDF = NULL,
                                      dimMat = NULL, maxDist = NULL,
                                      metrics = c('sensitivity', 'specificity', 'precision',
                                                  'accuracy', 'sizeProximity', 'scoreSpecificity'))
  return(accuracyBenchmarkMultiple(seuratObj, labelCol, geneSetNames, gsaMethods, labels,
                                   boundaryBenchmark, normSilDF, dimMat, maxDist, metrics))

#' Perform distribution accuracy benchmarks for multiple sets of GSA method scores
#' and class labels
#'
#' This function performs distribution accuracy benchmarks for multiple sets of GSA
#' method scores and class labels.
#'
#' @inheritParams accuracyBenchmarkMultiple
#' @inheritParams distributionBenchmark
#'
#' @return A list of lists of distribution benchmark data frames
#'
#' @export
#'
distributionBenchmarkMultiple <- function(seuratObj, labelCol, geneSetNames, gsaMethods,
                                          labels = geneSetNames,
                                          metrics = c('AUC', 'Gini', 'KS_Stat', 'PRAUC'))
  return(accuracyBenchmarkMultiple(seuratObj, labelCol, geneSetNames, gsaMethods, labels,
                                   distributionBenchmark, metrics))
