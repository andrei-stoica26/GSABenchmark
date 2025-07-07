#' Perform an accuracy benchmark on a set of gene set analysis method scores
#'
#' This function perform an accuracy benchmark on a set of gene set analysis method scores
#'
#' @inheritParams extractCellScores
#' @param computeMetricsFun Function used to compute metrics. Options are
#' computeBoundaryMetrics and computeDistributionMetrics.
#' @inheritParams computeOveralls
#' @param ... Additional parameters passed to computeMetricsFun.
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
#'
#' @return A benchmark data frame with cells as row names, and labels, GSA scores,
#' and class boundary metric scores as columns
#'
#' @export
#'
boundaryBenchmark <- function(seuratObj, labelCol, scoreCol, label)
  return(accuracyBenchmark(seuratObj, labelCol, scoreCol, label, computeBoundaryMetrics))

#' Perform a global evaluation benchmark on a set of GSA method scores
#'
#' This function performs a global evaluation benchmark on a set of GSA method scores.
#'
#' @inheritParams accuracyBenchmark
#' @inheritParams computeGlobalMetrics
#'
#' @return A one-row benchmark data frame with global metric scores as columns.
#'
#' @export
#'
globalBenchmark <- function(seuratObj, labelCol, scoreCol, label,
                            normSilDF = NULL, dimMat = NULL, maxDist = NULL)
  return(accuracyBenchmark(seuratObj, labelCol, scoreCol, label, computeGlobalMetrics,
                           normSilDF, dimMat, maxDist))

#' Compute the MCC based on the previously identified boundary cutoff
#'
#' This function computes the MCC based on the previously identified boundary cutoff
#'
#' @param boundaryBenchmarkLL A list of lists of data frames generated with boundaryBenchmarkMultiple
#'
#' @return A benchmark data frame with cells as row names, and labels, GSA scores,
#' and MCC scores as columns
#'
#' @export
#'
boundaryMCCBenchmark <- function(boundaryBenchmarkLL){
  mccValues <- c()
  geneSetNames <- names(boundaryBenchmarkLL)
  gsaMethods <- names(boundaryBenchmarkLL[[1]])
  for (setName in geneSetNames)
    mccValues <- c(mccValues, sapply(boundaryBenchmarkLL[[setName]], computeMCC))
  mccValues <- data.frame(matrix(mccValues, length(gsaMethods), length(geneSetNames)))
  colnames(mccValues) <- geneSetNames
  mccValues <- computeMethodMeans(mccValues, gsaMethods)
  return(mccValues)
}

#' Perform a direct MCC benchmark on a set of GSA method scores
#'
#' This function performs a direct MCC benchmark on a set of GSA
#' method scores
#'
#' @inheritParams accuracyBenchmark
#'
#' @return A benchmark data frame with cells as row names, and labels, GSA scores,
#' and MCC scores as columns
#'
#' @export
#'
directMCCBenchmark <- function(seuratObj, labelCol, scoreCol, label)
  return(accuracyBenchmark(seuratObj, labelCol, scoreCol, label, computeMCCMetric))
