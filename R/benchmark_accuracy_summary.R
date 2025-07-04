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
