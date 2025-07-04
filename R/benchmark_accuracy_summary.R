
#' Compares the identification of different labels through different gene set
#' analysis methods
#'
#' This function compares the identification of different labels through different gene set
#' analysis methods.
#'
#' @param benchmarkLL A list of lists of benchmark data frames
#'
#' @return Summary data frames
#'
#' @export
#'
benchmarkSummary <- function(benchmarkLL){
  markerNames <- names(benchmarkLL)
  gsaMethods <- names(benchmarkLL[[1]])
  tableCols <- colnames(benchmarkLL[[1]][[1]])
  if (length(intersect(tableCols, c('AUC', 'Gini', 'KS_Stat', 'PRAUC'))) > 0)
    metrics <- tableCols else
      metrics <- tableCols[seq(3, length(tableCols))]
  smr <- lapply(metrics, function(metric){
    df <- data.frame(Reduce(rbind, lapply(gsaMethods, function(method)
      sapply(markerNames, function(setName)
        benchmarkLL[[setName]][[method]][[metric]][1]))))
    df <- computeMethodMeans(df, gsaMethods)
    return(df)
  })
  smr <- addMetricSummary(smr, metrics, gsaMethods)
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
  df$score <- rowSums(df)
  df$score <- liver::minmax(df$score)
  df <- df[order(df$score, decreasing = TRUE), ]
  return(df)
}

