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
bindSummary <- function(smr, nAggDFs=0, nAvgCols=0){
    gsaMethods <- sort(rownames(smr[[1]]))
    return(do.call(cbind,
                   lapply(smr[seq(length(smr) - nAggDFs)],
                          function(x) x[gsaMethods, seq(ncol(x) - nAvgCols)])))
}

#' Compute aggregate ranks
#'
#' This function computes aggregate ranks from a list of summary data frames.
#'
#' @inheritParams mdsScoreSummary
#' @inheritParams geneSetRanks
#'
#' @return A data frame of aggregate summary results.
#'
#' @keywords internal
#'
aggregateRanks <- function(smr,
                           rankMethod = c('min', 'average', 'first',
                                          'last', 'random', 'max')){
    rankMethod <- match.arg(rankMethod)
    df <- do.call(cbind, list(bindSummary(smr$boundary, 2, 1),
                              bindSummary(smr$MCC, 0, 1),
                              bindSummary(smr$global, 2, 1)))
    df <- apply(df, 2, function(x) rank(-x, ties.method=rankMethod))
    df <- rankSummary(df)
    return(df)
}
