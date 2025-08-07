#' Extract ratio ranks list from a list of summary data frames
#'
#' This function extracts ratio (maximum over mean) ranks list from a list of
#' summary data frames.
#'
#' @inheritParams bindSummary
#'
#' @return List of ratio ranks.
#'
#' @noRd
#'
ratioRanks <- function(smr, nAggDFs=0, nAvgCols=0){
    metricNames <- names(smr)
    metricNames <- metricNames[seq(length(metricNames) - nAggDFs)]

    geneSetNames <- colnames(smr[[1]])
    geneSetNames <- geneSetNames[seq(length(geneSetNames) - nAvgCols)]

    titles <- metricTitles()

    ratioDF <- do.call(rbind, lapply(metricNames, function(metricName){
        df <- smr[[metricName]]
        do.call(rbind, lapply(geneSetNames, function(gsName) {
            colVals <- df[[gsName]]
            maxVal <- max(colVals)
            colSum <- sum(colVals)
            complMean <- (colSum - maxVal) / nrow(df)
            ratio <- maxVal / complMean
            res <- data.frame(Method = rownames(df)[which(colVals == maxVal)],
                              Metric = paste0(titles[metricName], ' (',
                                              str_to_title(gsName), ')'),
                              Ratio = ratio)
            return(res)
        }))
    }))
    ratioDF <- ratioDF[order(ratioDF$Ratio, decreasing=TRUE),]
    return(ratioDF)
}

#' Extract all ratio ranks list of summary data frames
#'
#' This function extracts all ratio ranks list from a list
#' of summary data frames.
#'
#' @inheritParams bindSummary
#'
#' @return List of gene set ranks.
#'
#' @noRd
#'
allRatioRanks <- function(smr){
    allRatioDF <- do.call(rbind, list(ratioRanks(smr$boundary, 2, 1),
                                      ratioRanks(smr$MCC, 0, 1),
                                      ratioRanks(smr$global, 2, 1)))
    allRatioDF <- allRatioDF[order(allRatioDF$Ratio, decreasing=TRUE),]
    return(allRatioDF)
}
