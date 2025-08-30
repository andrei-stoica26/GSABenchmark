#' Extract metric ranks list from a list of summary data frames
#'
#' This function extracts metric ranks from a list of summary data frames.
#'
#' @inheritParams geneSetRanks
#'
#' @return List of metric ranks.
#'
#' @keywords internal
#'
metricRanks <- function(smr, nAggDFs=0, nAvgCols=0, rankMethod='min')
    return(lapply(smr[seq(length(smr) - nAggDFs)],
                  function(df) {
                      df <- df[, seq(ncol(df) - nAvgCols)]
                      df <- apply(df, 2, function(x) rank(-x, ties.method=rankMethod))
                      df <- rankSummary(df)
                      return(df)
                  }))

#' Extract all metric ranks list of summary data frames
#'
#' This function extracts metric ranks from a list of summary data frames.
#'
#' @inheritParams metricRanks
#'
#' @return List of metric ranks.
#'
#' @keywords internal
#'
allMetricRanks <- function(smr, rankMethod='min'){
    df <- c(metricRanks(smr$boundary, 2, 1, rankMethod),
            metricRanks(smr$MCC, 0, 1, rankMethod),
            metricRanks(smr$global, 2, 1, rankMethod))
    return(df)
}
