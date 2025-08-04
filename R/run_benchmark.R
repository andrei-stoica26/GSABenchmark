#' Generate all benchmark results with some precomputed arguments.
#
#' This function generates all benchmark results by using precomputed values
#' for \code{normSilDF}, \code{dimMat} and \code{maxDist}.
#'
#' @details This function is designed to save some computational time when
#' multiple runs are sequentially performed for the same single-cell expression
#' object and cell identity classes, but different choices of gene sets
#' (so long as they still describe the same classes) or methods. If a single
#' run is planned rather than multiple ones, \code{runBenchmark} is a more
#' straightforward choice, as it takes care of generating
#' \code{normSilDF}, \code{dimMat} and \code{dimMat}.
#'
#' @inheritParams globalBenchmarkMultiple
#' @inheritParams runGSAMethods
#' @param geneSets A list of gene sets. If not \code{NULL} while
#' \code{efBenchmark} is \code{NULL}, the efficiency benchmark will be run.
#' @param efBenchmark A list of dataframes generated with efficiencyBenchmark.
#'
#' @return A list of benchmark results.
#'
#' @export
#'
allBenchmarkResults <- function(scObj,
                                labelCol,
                                geneSets,
                                gsaMethods,
                                geneSetNames,
                                checkLabels = TRUE,
                                normSilDF,
                                dimMat,
                                maxDist,
                                efBenchmark = NULL){
    x <- Sys.time()
    if(checkLabels)
        checkSetNames(scObj, labelCol, geneSetNames)

    message('Running class determination boundary benchmark...')
    boundaryRes <- boundaryBenchmarkMultiple(scObj,
                                             labelCol,
                                             gsaMethods,
                                             geneSetNames,
                                             checkLabels=FALSE,
                                             verbose=FALSE)
    boundarySmr <- benchmarkSummary(boundaryRes)

    message('Running MCC benchmark...')
    directMCCRes <- directMCCBenchmarkMultiple(scObj,
                                               labelCol,
                                               gsaMethods,
                                               geneSetNames,
                                               checkLabels=FALSE,
                                               verbose=FALSE)

    directMCCSmr <- benchmarkSummary(directMCCRes,
                                     summarizeMetrics=FALSE)
    mccSmr <- list(boundary = boundaryMCCBenchmark(boundaryRes),
                 direct = directMCCSmr)

    message('Running global evaluation benchmark...')
    globalRes <- globalBenchmarkMultiple(scObj,
                                         labelCol,
                                         gsaMethods,
                                         geneSetNames,
                                         checkLabels=FALSE,
                                         verbose=FALSE,
                                         normSilDF,
                                         dimMat,
                                         maxDist)
    globalSmr <- benchmarkSummary(globalRes)

    if(is.null(efBenchmark) & !is.null(geneSets)){
        message('Running efficiency benchmark. This may take a long time...')
        efBenchmark <- efficiencyBenchmark(scObj,
                                           labelCol,
                                           geneSets,
                                           gsaMethods,
                                           geneSetNames,
                                           checkLabels=FALSE)
    }

    message('Computing scored MDS summary...')
    mdsScoreSmr <- mdsScoreSummary(scObj, gsaMethods, geneSetNames,
                                   boundarySmr, mccSmr, globalSmr)

    message('Collating results...')
    smr <- list(boundary = boundarySmr,
                MCC = mccSmr,
                global = globalSmr,
                MDS = mdsScoreSmr)

    if(!is.null(efBenchmark) | !is.null(geneSets))
        smr$efficiency <- efBenchmark

    y <- Sys.time()
    print(y - x)
    return(smr)
}

#' Generate all benchmark results
#
#' This function performs the entire \code{GSABenchmark} pipeline.
#'
#' @details A wrapper around \code{allBenchmarkResults}. Slower for repeated
#' runs, but it does not require the user to manually
#' generate \code{normSilDF}, \code{dimMat} and \code{dimMat}.
#'
#' @inheritParams allBenchmarkResults
#'
#' @return A list of benchmark results.
#'
#' @export
#'
runBenchmark <- function(scObj,
                         labelCol,
                         geneSets,
                         gsaMethods,
                         geneSetNames){

    checkSetNames(scObj, labelCol, geneSetNames)
    scObj <- computeSilhouette(scObj, labelCol)
    message('Computing identity class-normalized silhouette...')
    normSilDF <- normalizeSilhouette(scObj, labelCol)
    message('Computing PCA-based distance matrix...')
    dimMat <- scPCAMat(scObj)
    maxDist <- max(dist(dimMat))
    res <- allBenchmarkResults(scObj,
                               labelCol,
                               geneSets,
                               gsaMethods,
                               geneSetNames,
                               checkLabels=FALSE,
                               normSilDF,
                               dimMat,
                               maxDist)
    return(res)
}
