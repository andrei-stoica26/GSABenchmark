#' Extract gene set ranks list from a list of summary data frames
#'
#' This function extracts gene set ranks list from a list of summary data
#' frames.
#'
#' @inheritParams bindSummary
#' @param rankMethod Rank method. Choose between 'min', 'average', 'first',
#' 'last', 'random' and 'max'. Default in 'min'.
#'
#' @return List of gene set ranks.
#'
#' @keywords internal
#'
geneSetRanks <- function(smr,
                         nAggDFs = 0,
                         nAvgCols = 0,
                         rankMethod = c('min', 'average', 'first',
                                        'last', 'random', 'max')){
    rankMethod <- match.arg(rankMethod, c('min', 'average', 'first',
                                          'last', 'random', 'max'))

    geneSetNames <- colnames(smr[[1]])
    geneSetNames <- geneSetNames[seq(length(geneSetNames) - nAvgCols)]
    gsaMethods <- sort(rownames(smr[[1]]))

    res <- setNames(lapply(geneSetNames, function(gsName){
        metrics <- names(smr)
        metrics <- metrics[seq(length(metrics) - nAggDFs)]
        df <- do.call(cbind, lapply(metrics, function(metric)
            setNames(smr[[metric]][gsaMethods, gsName, drop=FALSE],
                     metric)))
        df <- apply(df, 2, function(x) rank(-x, ties.method=rankMethod))
        return(df)
    }), geneSetNames)

    return(res)
}

#' Extract all gene set ranks list of summary data frames
#'
#' This function extracts all gene set ranks list from a list
#' of summary data frames.
#'
#' @inheritParams geneSetRanks
#'
#' @return List of gene set ranks.
#'
#' @keywords internal
#'
allGeneSetRanks <- function(smr, rankMethod = 'min'){
    boundaryRanks <- geneSetRanks(smr$boundary, 2, 1, rankMethod)
    mccRanks <- geneSetRanks(smr$MCC, 0, 1, rankMethod)
    globalRanks <- geneSetRanks(smr$global, 2, 1, rankMethod)
    geneSetNames <- names(boundaryRanks)

    res <- setNames(lapply(geneSetNames, function(gsName){
        df <- do.call(cbind, list(boundaryRanks[[gsName]],
                                  mccRanks[[gsName]],
                                  globalRanks[[gsName]]))
        df <- rankSummary(df)
        return(df)
    }), geneSetNames)

    return(res)
}
