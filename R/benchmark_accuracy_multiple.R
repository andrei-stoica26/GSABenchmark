#' Perform accuracy benchmarks for multiple sets of GSA method scores and class
#' labels
#'
#' This function performs accuracy benchmarks for multiple sets of GSA method
#' scores and class labels.
#'
#' @inheritParams accuracyBenchmark
#' @inheritParams runGSAMethods
#' @param geneSetNames Gene set names.
#' @param checkLabels Whether to check that \code{geneSetNames} is a subset of
#' the values in \code{labelCol}.
#' @param benchmarkFun Benchmark function. Must be either boundaryBenchmark or
#' distributionBenchmark.
#' @param verbose Whether to output a message whenever the benchmarking
#' begins for a gene set.
#' @param ... Additional paramters passed to benchmarkFun.
#'
#' @return A list of lists of benchmark data frames.
#'
#' @keywords internal
#'
accuracyBenchmarkMultiple <- function(scObj,
                                      labelCol,
                                      geneSetNames,
                                      gsaMethods,
                                      checkLabels = TRUE,
                                      verbose = FALSE,
                                      benchmarkFun,
                                      ...){

    if(checkLabels)
        checkSetNames(scObj, labelCol, geneSetNames)
    res <- lapply(seq_along(geneSetNames), function(i) {
        setName <- geneSetNames[i]
        if (verbose)
            message('Running benchmark on gene set: ',
                    geneSetNames[i], '...')
        setRes <- lapply(seq_along(gsaMethods), function(j){
            method <- gsaMethods[j]
            scoreCol <- paste0(method, '_', setName)
            df <- benchmarkFun(scObj,
                               labelCol,
                               scoreCol,
                               geneSetNames[i],
                               checkLabel=FALSE,
                               ...)
            return(df)
            })
        names(setRes) <- gsaMethods
        return(setRes)
    })
    names(res) <- geneSetNames
    return(res)
}

#' Perform boundary accuracy benchmarks for multiple sets of gene set analysis
#' method scores and class labels
#'
#' This function performs boundary accuracy benchmarks for multiple sets of
#' gene set analysis method scores and class labels.
#'
#' @inheritParams accuracyBenchmarkMultiple
#'
#' @return A list of lists of boundary benchmark data frames.
#'
#' @keywords internal
#'
boundaryBenchmarkMultiple <- function(scObj,
                                      labelCol,
                                      geneSetNames,
                                      gsaMethods,
                                      checkLabels = TRUE,
                                      verbose = FALSE)
    return(accuracyBenchmarkMultiple(scObj,
                                     labelCol,
                                     geneSetNames,
                                     gsaMethods,
                                     checkLabels,
                                     verbose,
                                     boundaryBenchmark))

#' Perform direct MCC benchmarks for multiple sets of GSA method scores
#' and class labels
#'
#' This function performs direct MCC benchmarks for multiple sets of GSA
#' method scores and class labels.
#'
#' @inheritParams accuracyBenchmarkMultiple
#'
#' @return A list of lists of boundary benchmark data frames.
#'
#' @keywords internal
#'
directMCCBenchmarkMultiple <- function(scObj,
                                       labelCol,
                                       geneSetNames,
                                       gsaMethods,
                                       checkLabels = TRUE,
                                       verbose = FALSE)
    return(accuracyBenchmarkMultiple(scObj,
                                     labelCol,
                                     geneSetNames,
                                     gsaMethods,
                                     checkLabels,
                                     verbose,
                                     directMCCBenchmark
                                     ))

#' Perform distribution accuracy benchmarks for multiple sets of GSA method
#' scores and class labels
#'
#' This function performs distribution accuracy benchmarks for multiple sets
#' of GSA method scores and class labels.
#'
#' @inheritParams accuracyBenchmarkMultiple
#' @inheritParams globalBenchmark
#'
#' @return A list of lists of distribution benchmark data frames.
#'
#' @keywords internal
#'
globalBenchmarkMultiple <- function(scObj,
                                    labelCol,
                                    geneSetNames,
                                    gsaMethods,
                                    checkLabels = TRUE,
                                    verbose = FALSE,
                                    normSilDF = NULL,
                                    dimMat = NULL,
                                    maxDist = NULL)
    return(accuracyBenchmarkMultiple(scObj,
                                     labelCol,
                                     geneSetNames,
                                     gsaMethods,
                                     checkLabels,
                                     verbose,
                                     globalBenchmark,
                                     normSilDF,
                                     dimMat,
                                     maxDist
                                     ))
