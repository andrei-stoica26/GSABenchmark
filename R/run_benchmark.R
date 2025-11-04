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
#' \code{normSilDF} and \code{dimMat}.
#'
#' @inheritParams globalBenchmarkMultiple
#' @inheritParams runGSAMethods
#' @param geneSets A list of gene sets. If not \code{NULL} while
#' \code{efBenchmark} is \code{NULL}, the efficiency benchmark will be run.
#' @param efBenchmark A list of data frames generated with efficiencyBenchmark.
#' @param runEFBenchmark Whether to run efficiency benchmark.
#' @param verbose Whether the output of the efficiency benchmark should be
#' verbose. Ignored if \code{runEFBenchmark} is \code{FALSE}.
#'
#' @return A list of benchmark results.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' scObj <- hammers::computeSilhouette(scObj, 'label')
#' normSilDF <- hammers::normalizeSilhouette(scObj, 'label')
#' dimMat <- hammers::scPCAMat(scObj)
#' maxDist <- max(dist(dimMat))
#' smr <- allBenchmarkResults(scObj, 'label', geneSets, c('CSOA', 'Zscore'),
#' TRUE, normSilDF, dimMat, maxDist, runEFBenchmark=FALSE)
#'
#' @export
#'
allBenchmarkResults <- function(scObj,
                                labelCol,
                                geneSets,
                                gsaMethods,
                                checkLabels = TRUE,
                                normSilDF = NULL,
                                dimMat = NULL,
                                maxDist = NULL,
                                efBenchmark = NULL,
                                runEFBenchmark = TRUE,
                                verbose = TRUE){
    geneSetNames <- names(geneSets)
    if(checkLabels)
        checkSetNames(scObj, labelCol, geneSetNames)

    message('Running class determination boundary benchmark...')
    boundaryRes <- boundaryBenchmarkMultiple(scObj,
                                             labelCol,
                                             geneSetNames,
                                             gsaMethods,
                                             checkLabels=FALSE,
                                             verbose=FALSE)
    boundarySmr <- benchmarkSummary(boundaryRes)

    message('Running MCC benchmark...')
    directMCCRes <- directMCCBenchmarkMultiple(scObj,
                                               labelCol,
                                               geneSetNames,
                                               gsaMethods,
                                               checkLabels=FALSE,
                                               verbose=FALSE)

    directMCCSmr <- benchmarkSummary(directMCCRes, TRUE, FALSE)
    mccSmr <- list(boundaryMCC = boundaryMCCBenchmark(boundaryRes),
                 directMCC = directMCCSmr)

    message('Constructing binary predictions...')
    predList <- binaryPred(scObj, labelCol, boundaryRes)

    message('Running global evaluation benchmark...')
    globalRes <- globalBenchmarkMultiple(scObj,
                                         labelCol,
                                         geneSetNames,
                                         gsaMethods,
                                         checkLabels=FALSE,
                                         verbose=FALSE,
                                         normSilDF,
                                         dimMat,
                                         maxDist)
    globalSmr <- globalBenchmarkSummary(globalRes,
                                        geneSetNames,
                                        gsaMethods,
                                        predList)

    smr <- list(boundary = boundarySmr,
                MCC = mccSmr,
                global = globalSmr,
                predictions = predList)

    if(!is.null(efBenchmark))
       smr$efficiency <- efBenchmark else{
           if(runEFBenchmark){
               message('Running efficiency benchmark.',
               ' This may take a long time...')
               smr$efficiency <- efficiencyBenchmark(scObj,
                                                  labelCol,
                                                  geneSets,
                                                  gsaMethods,
                                                  checkLabels=FALSE,
                                                  verbose=verbose)
           }}
    return(smr)
}

#' Generate all benchmark results
#
#' This function performs the entire \code{GSABenchmark} pipeline.
#'
#' @details A wrapper around \code{allBenchmarkResults}. Slower for repeated
#' runs, but it does not require users to manually generate \code{normSilDF}
#' and \code{dimMat}.
#'
#' @inheritParams allBenchmarkResults
#'
#' @return A list of benchmark results.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' smr <- runBenchmark(scObj, 'label', geneSets, c('CSOA', 'Zscore'), FALSE)
#'
#' @export
#'
runBenchmark <- function(scObj,
                         labelCol,
                         geneSets,
                         gsaMethods,
                         runEFBenchmark = TRUE,
                         verbose = TRUE){

    checkSetNames(scObj, labelCol, names(geneSets))
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
                               checkLabels=FALSE,
                               normSilDF,
                               dimMat,
                               maxDist,
                               runEFBenchmark=runEFBenchmark,
                               verbose=verbose)
    return(res)
}
