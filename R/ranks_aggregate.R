#' Extract results list of summary data frames
#'
#' This function extracts results from a list of summary data frames.
#'
#' @param smr List of summary data frames.
#' @param nAggDFs Number of data frames of aggregate results in the list of
#' summary data frames.
#' @param nAvgCols Number of average columns for each data frame in the list.
#'
#' @return A data frame of aggregate summary results.
#'
#' @keywords internal
#'
bindSummary <- function(smr, nAggDFs=0, nAvgCols=0)
    return(do.call(cbind,
                   lapply(smr[seq(length(smr) - nAggDFs)],
                          function(x) x[, seq(ncol(x) - nAvgCols)])))

#' Compute aggregate ranks
#'
#' This function computes aggregate ranks from a list of summary data frames.
#'
#' @inheritParams mdsScoreSummary
#' @param nAggCols Number of columns of aggregate results in the list of
#' summary data frames.
#' @param nAvgCols Number of average columns for each data frame in the list.
#'
#' @return A data frame of aggregate summary results.
#'
#' @export
#'
aggregateRanks <- function(smr){
    df <- do.call(cbind, list(bindSummary(smr$boundary, 2, 1),
                              bindSummary(smr$MCC, 0, 1),
                              bindSummary(smr$global, 2, 1)))
    df <- apply(df, 2, function(x) rank(-x))
    df <- rankSummary(df)
    return(df)
}


