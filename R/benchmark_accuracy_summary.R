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
    'metric. Set summarizeMetrics to FALSE')
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
#' @param boundarySmr Class boundary determination summary.
#' @param mccSmr MCC summary.
#' @param globalSmr Global evaluation summary.
#'
#' @return A list containing a data frame with three columns (two MDS
#' coordinates and method score) for each gene set.
#'
#' @export
#'
mdsScoreSummary <- function(scObj,
                            gsaMethods,
                            geneSetNames,
                            boundarySmr,
                            mccSmr,
                            globalSmr){

    mdsScoreSmr <- lapply(geneSetNames, function(gsName){
    setDF <- metadataDF(scObj)[, paste0(gsaMethods, '_', gsName)]
    colnames(setDF) <- gsaMethods

    distMat <- as.matrix(stats::dist(base::t(setDF)))
    mdsMat <- as.data.frame(cmdscale(distMat))

    colnames(mdsMat) <- c('x', 'y')

    summaryMat <- do.call(cbind, list(
        geneSetMetrics(boundarySmr, gsaMethods, gsName, 2),
        geneSetMetrics(mccSmr, gsaMethods, gsName, 0),
        geneSetMetrics(globalSmr, gsaMethods, gsName, 2)
        ))

    mdsMat$Score <- rowMeans(summaryMat)
    mdsMat$nn <- nearestNeighbors(distMat)
    mdsMat <- mdsMat[order(mdsMat$Score, decreasing=TRUE), ]
    return(mdsMat)
    })

    names(mdsScoreSmr) <- geneSetNames
    return(mdsScoreSmr)
}
