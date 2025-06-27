
#' Add the score-based metrics to the benchmark
#'
#' This function adds the score-based metrics to the benchmark
#'
#' @param seuratObj A Seurat object
#' @param labelCol The Seurat metadata column containing the ground truth annotation
#' @param scoreCol The Seurat metadata column containing the gene set analysis method score
#' @param label The identity assessed from the labelCol column
#' @param nMetrics Number of metrics
#' @inheritParams addCentralityMetrics
#'
#' @return A benchmark data frame
#'
#' @export
#'
identityClassMatch <- function(seuratObj, labelCol, scoreCol, normSilDF, dimMat, label, nMetrics){
  df <- seuratObj@meta.data[, c(labelCol, scoreCol)]
  df$label <- as.integer(df[, 1] %in% label)
  df <- df[order(df[, scoreCol], decreasing=TRUE),]

  df <- addCellCountMetrics(df)
  df <- addScoreMetric(df)
  df <- addCentralityMetrics(df, normSilDF, dimMat)

  df$accuracy <- rowMeans(df[, 2 + seq_len(nMetrics)])
  df <- df[order(df$accuracy, decreasing=TRUE),]
  return(df)
}

#' Performs a comparison of different gene set analysis methods using annotations
#'
#' This function compares the ability of different gene set analysis methods to
#' correctly identity cells with a certain label.
#'
#' @param seuratObj A Seurat object
#' @param labelCol The column of labels
#' @param markerNames The names of marker sets
#' @param gsaMethods The assessed gene set analysis methods
#' @param labels The assessed labels. Must be a subset of the values in labelCol
#' @inheritParams identityClassMatch
#'
#' @return A benchmark data frame
#'
#' @export
#'
identityClassBenchmark <- function(seuratObj, labelCol, markerNames, gsaMethods, normSilDF, dimMat, nMetrics = 6, labels = markerNames){
  res <- lapply(seq_along(markerNames), function(i) {
    setName <- markerNames[i]
    setRes <- lapply(seq_along(gsaMethods), function(j){
      method <- gsaMethods[j]
      scoreCol <- paste0(method, '_', setName)
      df <- identityClassMatch(seuratObj, labelCol, scoreCol, normSilDF, dimMat, labels[i], nMetrics)
      return(df)
    })
    names(setRes) <- gsaMethods
    return(setRes)
  })
  names(res) <- labels
  return(res)
}

#' Compares the identification of different labels through different gene set
#' analysis methods
#'
#' This function compares the identification of different labels through different gene set
#' analysis methods.
#'
#' @param icBenchmark A list of benchmark data frames
#' @inheritParams identityClassMatch
#'
#' @return Summary data frames
#'
#' @export
#'
identityClassBenchmarkSummary <- function(icBenchmark, nMetrics = 6){
  markerNames <- names(icBenchmark)
  gsaMethods <- names(icBenchmark[[1]])
  metrics <- colnames(icBenchmark[[1]][[1]])[2 + seq_len(nMetrics + 1)]
  smr <- lapply(metrics, function(metric){
    df <- data.frame(Reduce(rbind, lapply(gsaMethods, function(method)
      sapply(markerNames, function(setName)
        icBenchmark[[setName]][[method]][[metric]][1]))))
    rownames(df) <- gsaMethods
    df$avg <- rowMeans(df)
    df <- df[order(df$avg, decreasing=TRUE), ]
    return(df)
  })
  names(smr) <- metrics
  df <- data.frame(lapply(smr, function(x) x[gsaMethods, ]$avg))
  rownames(df) <- gsaMethods
  df <- df[order(df$accuracy, decreasing=TRUE), ]
  smr <- c(smr, list(total = df))
  return(smr)
}

#' Calculate the summary of the results based on the minmax normalization
#'
#' This function strips the last column of the results (accuracy), minmax-normalizes
#' the results of each metric, then defines the accuracy as the minmax-normalized
#' mean of the results
#'
#' @param df A data frame with metric results and accuracy
#'
#' @return A minmax-normalized summary data frame
#'
#' @export
#'
minmaxSummary <- function(df){
  df <- df[ , seq_len(ncol(df) - 1)]
  df <- data.frame(apply(df, 2, liver::minmax))
  df$accuracy <- rowSums(df)
  df$accuracy <- liver::minmax(df$accuracy)
  df <- df[order(df$accuracy, decreasing = TRUE), ]
  return(df)
}
