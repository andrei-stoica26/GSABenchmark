#' Generate all benchmark results
#
#' This function generates all benchmark results.
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
