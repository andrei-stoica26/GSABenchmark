#' @importFrom jaccard jaccard
#'
NULL

#' Compares the identification of different labels through different gene set
#' analysis methods
#'
#' This function compares the identification of different labels through
#' different gene set analysis methods.
#'
#' @param benchmarkLL A list of lists of benchmark data frames.
#' @inheritParams addMethodInfo
#' @inheritParams addMetricSummary
#'
#' @return A list of summary data frames.
#'
#' @keywords internal
#'
benchmarkSummary <- function(benchmarkLL,
                             doAverage = TRUE,
                             doSummarize = TRUE){
    geneSetNames <- names(benchmarkLL)
    gsaMethods <- names(benchmarkLL[[1]])
    tableCols <- colnames(benchmarkLL[[1]][[1]])
    if('MCC' %in% tableCols & doSummarize)
        stop('You are summarizing an MCC benchmark, which only has one ',
             'metric. Set doSummarize to FALSE.')
    if (length(intersect(tableCols, c('AUC', 'PRAUC'))) > 0)
        metrics <- tableCols else
            metrics <- tableCols[seq(3, length(tableCols))]
    smr <- lapply(metrics, function(metric){
        df <- data.frame(Reduce(rbind, lapply(gsaMethods, function(method)
            vapply(geneSetNames, function(setName)
                benchmarkLL[[setName]][[method]][[metric]][1], numeric(1)))))
        df <- addMethodInfo(df, gsaMethods, doAverage)
        return(df)
    })
    smr <- addMetricSummary(smr, metrics, gsaMethods, doSummarize)
    return(smr)
}

#' Compares the identification of different labels through different gene set
#' analysis methods
#'
#' This function compares the identification of different labels through
#' different gene set analysis methods.
#'
#' @param globalBenchmarkLL A list of lists of benchmark data frames.
#'
#' @return A list of summary data frames.
#'
#' @keywords internal
#'
globalBenchmarkSummary <- function(globalBenchmarkLL,
                                   geneSetNames,
                                   gsaMethods,
                                   predList){
    globalSmr <- benchmarkSummary(globalBenchmarkLL, TRUE, FALSE)
    globalSmr <- binaryPredMetric(globalSmr, geneSetNames, gsaMethods,
                                  predList, 'labJaccard', jaccard)
    globalSmr <- binaryPredMetric(globalSmr, geneSetNames, gsaMethods,
                                  predList, 'labCosine', numCosine)
    globalSmr <- addGlobalAverage(globalSmr, gsaMethods)
    globalSmr <- addMetricSummary(globalSmr, names(globalSmr), gsaMethods)
    return(globalSmr)
}
