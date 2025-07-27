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
#' @param efBenchmark A list of dataframes generated with efficiencyBenchmark.
#'
#' @return A list of benchmark results.
#'
#' @export
#'
allBenchmarkResults <- function(scObj,
                                labelCol,
                                geneSets,
                                geneSetNames,
                                gsaMethods,
                                labels = geneSetNames,
                                normSilDF,
                                dimMat,
                                maxDist,
                                efBenchmark = NULL){
  x <- Sys.time()
  message('Running class determination boundary benchmark...')
  boundaryRes <- boundaryBenchmarkMultiple(scObj,
                                           labelCol,
                                           geneSetNames,
                                           gsaMethods,
                                           labels = geneSetNames)
  boundarySmr <- benchmarkSummary(boundaryRes)
  message('Running MCC benchmark...')
  directMCCRes <- directMCCBenchmarkMultiple(scObj,
                                             labelCol,
                                             geneSetNames,
                                             gsaMethods)
  directMCCSmr <- benchmarkSummary(directMCCRes, FALSE)
  mccSmr <- list(boundary = boundaryMCCBenchmark(boundaryRes),
                 direct = directMCCSmr)
  message('Running global evaluation benchmark...')
  globalRes <- globalBenchmarkMultiple(scObj,
                                       labelCol,
                                       geneSetNames,
                                       gsaMethods,
                                       labels = geneSetNames,
                                       normSilDF,
                                       dimMat,
                                       maxDist)
  globalSmr <- benchmarkSummary(globalRes)
  if(is.null(efBenchmark)){
    message('Running efficiency benchmark. This may take a long time...')
    efBenchmark <- efficiencyBenchmark(scObj, geneSets,
                                       geneSetNames, gsaMethods)
  } else message('Loading efficiency benchmark...')
  message('Collating results...')
  smr <- list(boundary = boundarySmr,
              MCC = mccSmr,
              global = globalSmr,
              efficiency = efBenchmark)
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
                         geneSetNames,
                         gsaMethods,
                         labels = geneSetNames){
    scObj <- computeSilhouette(scObj, labelCol)
    message('Computing identity class-normalized silhouette...')
    normSilDF <- normalizeSilhouette(scObj, labelCol)
    message('Computing PCA-based distance matrix...')
    dimMat <- scPCAMat(scObj)
    maxDist <- max(dist(dimMat))
    res <- allBenchmarkResults(scObj,
                               labelCol,
                               geneSets,
                               geneSetNames,
                               gsaMethods,
                               labels,
                               normSilDF,
                               dimMat,
                               maxDist)
    return(res)
}
