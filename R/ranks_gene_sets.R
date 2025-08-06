#' Extract gene set ranks list from a list of summary data frames
#'
#' This function extracts gene set ranks list from a list of summary data frames.
#'
#' @inheritParams bindSummary
#'
#' @return List of gene set ranks.
#'
#' @noRd
#'
geneSetRanks <- function(smr, nAggDFs=0, nAvgCols=0){
    geneSetNames <- colnames(smr[[1]])
    geneSetNames <- geneSetNames[seq(length(geneSetNames) - nAvgCols)]

    res <- lapply(geneSetNames, function(gsName){
        metrics <- names(smr)
        metrics <- metrics[seq(length(metrics) - nAggDFs)]
        df <- do.call(cbind, lapply(metrics, function(metric)
            setNames(smr[[metric]][, gsName, drop=FALSE],
                     metric)))
        df <- apply(df, 2, function(x) rank(-x))
        return(df)
    })

    names(res) <- geneSetNames
    return(res)
}

#' Extract all gene set ranks list of summary data frames
#'
#' This function extracts all gene set ranks list from a list
#' of summary data frames.
#'
#' @inheritParams bindSummary
#'
#' @return List of gene set ranks.
#'
#' @noRd
#'
allGeneSetRanks <- function(smr){
    boundaryRanks <- collectGeneSetRanks(smr$boundary, 2, 1)
    mccRanks <- collectGeneSetRanks(smr$MCC, 0, 1)
    globalRanks <- collectGeneSetRanks(smr$global, 2, 1)
    geneSetNames <- names(boundaryRanks)
    res <- lapply(geneSetNames, function(gsName){
        df <- do.call(cbind, list(boundaryRanks[[gsName]],
                                  mccRanks[[gsName]],
                                  globalRanks[[gsName]]))
        df <- rankSummary(df)
        return(df)
    })
    names(res) <- geneSetNames
    return(res)
}
