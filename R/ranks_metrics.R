#' Extract metric ranks list from a list of summary data frames
#'
#' This function extracts metric ranks from a list of summary data frames.
#'
#' @inheritParams bindSummary
#'
#' @return List of metric ranks.
#'
#' @noRd
#'
metricRanks <- function(smr, nAggDFs=0, nAvgCols=0)
    return(lapply(smr[seq(length(smr) - nAggDFs)],
                  function(df) {
                      df <- df[, seq(ncol(df) - nAvgCols)]
                      df <- apply(df, 2, function(x) rank(-x))
                      df <- rankSummary(df)
                      return(df)
                  }))

#' Extract all metric ranks list of summary data frames
#'
#' This function extracts metric ranks from a list of summary data frames.
#'
#' @inheritParams bindSummary
#'
#' @return List of metric ranks.
#'
#' @noRd
#'
allMetricRanks <- function(smr){
    df <- c(metricRanks(smr$boundary, 2, 1),
            metricRanks(smr$MCC, 0, 1),
            metricRanks(smr$global, 2, 1))
    return(df)
}
