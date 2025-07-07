#' Generate all benchmark results
#
#' This function generates all benchmark results.
#'
#' @inheritParams globalBenchmarkMultiple
#' @param efBenchmark A list of dataframes generated with efficiencyBenchmark.
#'
#' @return A list of benchmark results.
#'
#' @export
#'
allBenchmarkResults <- function(seuratObj, labelCol, gsaMethods, labels = geneSetNames,
                                normSilDF, dimMat, maxDist, efBenchmark = NULL){
  x <- Sys.time()
  message('Running class determination boundary benchmark...')
  boundaryRes <- boundaryBenchmarkMultiple(seuratObj, 'label', geneSetNames, gsaMethods, labels = geneSetNames)
  boundarySmr <- benchmarkSummary(res)
  message('Running MCC benchmark...')
  directMCCRes <- directMCCBenchmarkMultiple(seuratObj, 'label', geneSetNames, gsaMethods)
  directMCCSmr <- benchmarkSummary(directMCCRes, F)
  mccSmr <- list(boundary = boundaryMCCBenchmark(boundaryRes),
                 direct = directMCCSmr)
  message('Running global evaluation benchmark...')
  globalRes <- globalBenchmarkMultiple(seuratObj, 'label', geneSetNames, gsaMethods, labels = geneSetNames,
                                       normSilDF, dimMat, maxDist)
  globalSmr <- benchmarkSummary(globalRes)
  if(is.null(efBenchmark)){
    message('Running efficiency benchmark. This may take a long time...')
    efBenchmark <- efficiencyBenchmark(seuratObj, geneSets, geneSetNames, gsaMethods)
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
