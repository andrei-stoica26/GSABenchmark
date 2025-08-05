#' Extract gene set results from a list of summary data frames
#'
#' This function extracts gene set results from a list of summary data frames.
#'
#' @param smr List of summary data frames.
#' @inheritParams mdsScoreSummary
#' @param gsName Gene set name.
#' @param nAggCols Number of columns of aggregate results in the list of
#' summary data frames.
#'
#' @return A matrix of gene set summary results.
#'
#' @noRd
#'
geneSetMetrics <- function(smr, gsaMethods, gsName, nAggCols)
    return(vapply(smr[seq(length(smr) - nAggCols)],
                  function(x) x[gsaMethods, gsName],
                  numeric(length(gsaMethods))))

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
