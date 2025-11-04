#' Construct binary prediction data frames for the gene set analysis methods
#'
#' This function used the boundary benchmark results to construct binary
#' prediction data frames for the gene set analysis methods.
#'
#' @inheritParams extractCellScores
#' @inheritParams boundaryMCCBenchmark
#'
#' @return A list of binary data frames.
#'
#' @noRd
#'
binaryPred <- function(scObj, labelCol, boundaryBenchmarkLL){
    geneSetNames <- names(boundaryBenchmarkLL)
    gsaMethods <- names(boundaryBenchmarkLL[[1]])

    predList <- setNames(lapply(geneSetNames, function(gsName){
        res <- data.frame(label = as.integer(scCol(scObj, labelCol) == gsName))
        rownames(res) <- colnames(scObj)
        benchmarkRes <- boundaryBenchmarkLL[[gsName]]
        for (method in gsaMethods){
            cutoff <- benchmarkRes[[method]][1, 2]
            res[[method]] <- as.integer(scCol(scObj,
                                              paste0(method, '_', gsName)) >= cutoff)
        }
        return(res)
    }), geneSetNames)

    return(predList)
}

#' Generate metric results for a gene set using binary predictions
#'
#' This function generates metric results for a gene set using
#' binary predictions.
#'
#' @inheritParams runGSAMethods
#' @param gsName Gene set name.
#' @param predList A list of data frames of the same length as the number of
#' gene sets, in which predictions made by each gene set analysis method are
#' represented as binary columns.
#' @param metricFun Metric function. Choose between \code{jaccard}
#' and \code{numCosine}.
#'
#' @return A numeric vector.
#'
#' @keywords internal
#'
binaryPredMetricCore <- function(gsName, gsaMethods, predList, metricFun){
    predDF <- predList[[gsName]]
    return(vapply(gsaMethods, function(method)
        metricFun(predDF[, 1], predDF[, method]), numeric(1)))
}

#' Generate metric results using binary predictions
#'
#' This function generates metric results using binary predictions.
#'
#' @inheritParams addGlobalAverage
#' @inheritParams binaryPredMetricCore
#' @param smrCol Summary column to be added
#'
#' @return A global summary object with an added binary prediction data frame.
#'
#' @keywords internal
#'
binaryPredMetric <- function(globalSmr, geneSetNames, gsaMethods,
                             predList, smrCol, metricFun){
    gsResList <- setNames(lapply(geneSetNames, function(gsName)
        binaryPredMetricCore(gsName, gsaMethods, predList, metricFun)),
        geneSetNames)
    df <- data.frame(do.call(cbind, gsResList))
    globalSmr[[smrCol]] <- addMethodInfo(df, gsaMethods)
    return(globalSmr)
}
