#' Compares the identification of different labels through different gene set
#' analysis methods
#'
#' This function compares the identification of different labels through
#' different gene set analysis methods.
#'
#' @param benchmarkLL A list of lists of benchmark data frames.
#' @param summarizeMetrics Whether to add a metric summary. Must be set to
#' \code{FALSE} when summarizing a MCC benchmark list of lists.
#'
#' @return A list of summary data frames.
#'
#' @keywords internal
#'
benchmarkSummary <- function(benchmarkLL, summarizeMetrics = TRUE){
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
            vapply(geneSetNames, function(setName)
                benchmarkLL[[setName]][[method]][[metric]][1], numeric(1)))))
        df <- computeMethodMeans(df, gsaMethods)
        return(df)
    })
    if(summarizeMetrics)
        return(addMetricSummary(smr, metrics, gsaMethods))
    return(smr[[1]])
}
