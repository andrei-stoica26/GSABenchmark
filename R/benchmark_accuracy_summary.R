#' Compares the identification of different labels through different gene set
#' analysis methods
#'
#' This function compares the identification of different labels through
#' different gene set analysis methods.
#'
#' @param benchmarkLL A list of lists of benchmark data frames.
#' @param summarizeMetrics Whether to add a metric summary. Must be set to
#' FALSE when summarizing a MCC benchmark list of lists.
#'
#' @return Summary data frames.
#'
#' @export
#'
benchmarkSummary <- function(benchmarkLL, summarizeMetrics=TRUE){
    geneSetNames <- names(benchmarkLL)
    gsaMethods <- names(benchmarkLL[[1]])
    tableCols <- colnames(benchmarkLL[[1]][[1]])
    if('MCC' %in% tableCols & summarizeMetrics)
        stop('You are summarizing an MCC benchmark, which only has one ',
             'metric. Set summarizeMetrics to FALSE.')
    if (length(intersect(tableCols, c('AUC', 'PRAUC'))) > 0)
        metrics <- tableCols else
            metrics <- tableCols[seq(3, length(tableCols))]
    smr <- lapply(metrics, function(metric){
        df <- data.frame(Reduce(rbind, lapply(gsaMethods, function(method)
            sapply(geneSetNames, function(setName)
                benchmarkLL[[setName]][[method]][[metric]][1]))))
        df <- computeMethodMeans(df, gsaMethods)
        return(df)
    })
    if(summarizeMetrics)
        return(addMetricSummary(smr, metrics, gsaMethods))
    return(smr[[1]])
}

#' Compute the scored MDS summary
#'
#' This function computes the score MDS summary. For each gene set, MDS is
#' performed for the cell-level method predicted scores. The resulting
#' two-column data frame is appended the average accuracy benchmark score
#' of each method, obtained on the gene set.
#'
#' @inheritParams runBenchmark
#' @param smr List containing boundary, MCC and global summary data frames.
#'
#' @return A list containing a data frame with three columns (two MDS
#' coordinates and method score) for each gene set.
#'
#' @export
#'
mdsScoreSummary <- function(scObj,
                            gsaMethods,
                            geneSetNames,
                            smr){

    mdsScoreSmr <- lapply(geneSetNames, function(gsName){
    setDF <- metadataDF(scObj)[, paste0(gsaMethods, '_', gsName)]
    colnames(setDF) <- gsaMethods

    distMat <- (1 - cor(setDF, method='spearman')) / 2
    mdsMat <- as.data.frame(cmdscale(distMat))

    colnames(mdsMat) <- c('x', 'y')

    summaryMat <- do.call(cbind, list(
        geneSetMetrics(smr$boundary, gsaMethods, gsName, 2),
        geneSetMetrics(smr$MCC, gsaMethods, gsName, 0),
        geneSetMetrics(smr$global, gsaMethods, gsName, 2)
        ))

    mdsMat$Score <- rowMeans(summaryMat)
    mdsMat$nn <- nearestNeighbors(distMat)
    mdsMat <- mdsMat[order(mdsMat$Score, decreasing=TRUE), ]
    return(mdsMat)
    })

    names(mdsScoreSmr) <- geneSetNames
    return(mdsScoreSmr)
}

#' Compute benchmark ranks
#'
#' This function computes benchmark ranks from a list of summary data frames.
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
benchmarkRanks <- function(smr){
    df <- do.call(cbind, list(bindSummary(smr$boundary, 2, 1),
                              bindSummary(smr$MCC),
                              bindSummary(smr$global, 2, 1)))
    df <- apply(df, 2, function(x) rank(-x, ties.method='min'))
    df <- rankSummary(df)
    return(df)
}
