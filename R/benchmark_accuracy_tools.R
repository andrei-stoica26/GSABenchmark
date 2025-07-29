#' @importFrom hammers centerOfMass computeSilhouette joinCharCombs metadataDF metadataNames normalizeSilhouette scCol scPCAMat scExpMat
#' @importFrom plyr count
#'
NULL

#' Extract GSA method scores and truth labels from single-cell
#' expression object.
#'
#' This function extracts GSA method scores and truth labels from a single-cell
#' expression object.
#'
#' @param scObj A \code{Seurat} or \code{SingleCellExperiment} object.
#' @param labelCol The metadata column containing the ground truth annotation.
#' @param scoreCol The metadata column containing the gene set analysis method
#' score.
#' @param label The identity assessed from the labelCol column.
#'
#' @return A data frame with two columns: truth labels (1 or 0) and
#' GSA method scores,
#'
#' @export
#'
extractCellScores <- function(scObj, labelCol, scoreCol, label){
  df <- metadataDF(scObj)[, c(labelCol, scoreCol)]
  colnames(df)[1] <- 'label'
  df$label <- as.integer(df[, 1] %in% label)
  df <- df[order(df[, scoreCol], decreasing=TRUE),]
  return(df)
}

#' Condense repeated GSA method scores
#'
#' This function condenses repeated GSA method scores, summing labels and recording
#' frequencies in the process.
#'
#' @param df A data frame with two columns: truth labels (1 or 0) and GSA method scores.
#'
#' @return A condensed scores data frame
#'
condenseRepeatedScores <- function(df){
  denseDF <- plyr::count(df[, 2])
  denseDF$label <- sapply(denseDF[, 1], function(score) sum(df[, 1][which(df[, 2] == score)]))
  denseDF <- denseDF[, c(3, 1, 2)]
  denseDF <- denseDF[order(denseDF[, 2], decreasing=TRUE), ]
  colnames(denseDF)[2] <- c(colnames(df)[2])
  return(denseDF)
}

#' Add overall scores to a benchmark data frame
#'
#' This function adds overlap scores to a benchmark data frame.
#'
#' @param df Benchmark data frame.
#' @param startCol Column where metric scores start.
#'
#' @return The input data frame sorted decreasingly by the newly added overall column.
#'
#'
computeMetricMeans <- function(df, startCol){
  df$avg <- rowMeans(df[, seq(startCol, ncol(df))])
  df <- df[order(df$avg, decreasing=TRUE),]
  return(df)
}

#' Extends summary by addding overall results for each metric
#'
#' This function extends summary by addding overall results for each metric.
#'
#' @param smr List of result data frames for each metrics. Each data frame contains
#' the results for each tested gene set analysis method for each gene set for the
#' corresponding gmethod.
#' @param metrics Metrics.
#' @inheritParams computeMethodMeans
#'
#' @return Extended summary list with an additional data frame showing the average
#' results obtained for each metric.
#'
addMetricSummary <- function(smr, metrics, gsaMethods){
  names(smr) <- metrics
  df <- data.frame(lapply(smr, function(x) x[gsaMethods, ]$avg))
  rownames(df) <- gsaMethods
  df <- df[order(df$avg, decreasing=TRUE), ]
  smr <- c(smr, list(metricSummary=df))
  return(smr)
}

#' Add means for metric results data frame
#'
#' This function adds means for a data frame of metric results.
#'
#' @param df A data frame where the values represent the scores obtained by a gene
#' set analysis method (row) on a gene set (column) for a metric.
#' @param gsaMethods Gene set analysis methods.
#'
#' @return A metric results data frame with added means, sorted decreasingly by
#' these means.
#'
#' @export
#'
computeMethodMeans <- function(df, gsaMethods){
  rownames(df) <- gsaMethods
  df$avg <- rowMeans(df)
  df <- df[order(df$avg, decreasing=TRUE), ]
  return(df)
}

#' Extract gene set results from a list of summary data frames
#'
#' This function extracts gene set results from a list of summary data frames.
#'
#' @param smr List of summary data frames.
#' @inheritParams mdsScoreSummary
#' @param gsName Gene set name.
#' @param nAggCols Number of columns of aggregate results in the list of
#' summary data frames.
#'
#' @return A matrix of gene set summary results.
#'
#' @noRd
#'
geneSetMetrics <- function(smr, gsaMethods, gsName, nAggCols)
    return(vapply(smr[seq(length(smr) - nAggCols)],
                  function(x) x[gsaMethods, gsName],
                  numeric(length(gsaMethods))))
