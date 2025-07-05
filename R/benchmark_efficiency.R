#' @importFrom stringr str_remove str_split str_sub str_replace
#' @importFrom peakRAM peakRAM
#'
NULL

#' Extract running times and peak memory usage for gene set analysis methods
#'
#' This function extracts running times and peak memory usage for gene set analysis
#' methods.
#'
#' @inheritParams runGSAMethods
#'
#' @return A list of two data frames, the first comprising running times in seconds,
#' the second comprising peak memory usage in mebibytes.
#'
#' @export
#'
efficiencyBenchmark <- function(seuratObj, geneSets, geneSetNames, gsaMethods){
  elapsedTime <- c()
  peakMemory <- c()
  for (i in seq_along(geneSets))
    for (method in gsaMethods){
      methodCall <- paste0('run', method)
      message(paste0('Computing running time and peak memory usage for ', method, ' on ', geneSetNames[i], ' genes...'))
      df <- peakRAM(silently_run(do.call(methodCall, list(seuratObj, geneSets[[i]]))))
      message(paste0('Time elapsed: ', round(df$Elapsed_Time_sec, 4), ' seconds'))
      elapsedTime <- c(elapsedTime, df$Elapsed_Time_sec)
      peakMemory <- c(peakMemory, df$Peak_RAM_Used_MiB)
    }
  time <- tabulateVector(elapsedTime, gsaMethods, geneSetNames, TRUE)
  space <- tabulateVector(peakMemory, gsaMethods, geneSetNames, TRUE)
  return(list(time = time, space = space))
}
