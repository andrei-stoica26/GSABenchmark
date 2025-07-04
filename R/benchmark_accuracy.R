#' Perform an accuracy benchmark on a set of GSA method scores
#'
#' This function perform an accuracy benchmark on a set of GSA method scores
#'
#' @inheritParams extractCellScores
#' @param computeMetricsFun Function used to compute metrics. Options are
#' computeBoundaryMetrics and computeDistributionMetrics
#' @inheritParams computeOveralls
#' @param ... Additional paramters passed to fun.
#'
#' @return A benchmark data frame
#'
#'
accuracyBenchmark <- function(seuratObj, labelCol, scoreCol, label, computeMetricsFun, ...){
  df <- extractCellScores(seuratObj, labelCol, scoreCol, label)
  df <- computeMetricsFun(df, ...)
  return(df)
}

#' Perform a class boundary accuracy benchmark on a set of GSA method scores
#'
#' This function performs a class boundary accuracy benchmark on a set of GSA
#' method scores
#'
#' @inheritParams accuracyBenchmark
#' @inheritParams computeBoundaryMetrics
#'
#' @return A benchmark data frame with cells as row names, and labels, GSA scores,
#' and class boundary metric scores as columns
#'
#' @export
#'
boundaryBenchmark <- function(seuratObj, labelCol, scoreCol, label,
                              normSilDF = NULL, dimMat = NULL, maxDist = NULL,
                              metrics = c('sensitivity', 'specificity', 'precision',
                                          'accuracy', 'sizeProximity', 'scoreSpecificity')
                                          )
  return(accuracyBenchmark(seuratObj, labelCol, scoreCol, label, computeBoundaryMetrics,
                           normSilDF, dimMat, maxDist, metrics))

#' Perform a distribution benchmark on a set of GSA method scores
#'
#' This function performs a distribution benchmark on a set of GSA method scores.
#'
#' @inheritParams accuracyBenchmark
#' @inheritParams computeDistributionMetrics
#'
#' @return A one-row benchmark data frame with distribution metric scores as columns.
#'
#' @export
#'
distributionBenchmark <- function(seuratObj, labelCol, scoreCol, label, metrics = c('AUC', 'Gini', 'KS_Stat', 'PRAUC'))
  return(accuracyBenchmark(seuratObj, labelCol, scoreCol, label, computeDistributionMetrics, metrics))
