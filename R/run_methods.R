#' @importFrom fabR silently_run
#' @importFrom logr log_close log_open log_print
NULL


#' Run gene set analysis methods
#'
#' This function runs the gene set analysis methods.
#'
#' @param seuratObj A Seurat object
#' @param markerList A list of gene sets
#' @param markerNames The names of the gene sets
#' @param gsaMethods List of gene set analysis methods
#' @param logFile Name of the log file saving running times
#'
#' @return A Seurat object with the results of the runs stored as metadata columns
#'
#' @export
#'
runGSAMethods <- function(seuratObj, markerList, markerNames, gsaMethods, logFile = 'gsa'){
  log_open(logFile, logdir=FALSE, show_notes=FALSE, compact=TRUE, header=FALSE)
  for (i in seq_along(markerList)){
    setName <- markerNames[i]
    for (j in seq_along(gsaMethods)){
      method <- gsaMethods[j]
      runningMsg <- paste0('Running ', method, ' for ', setName, ' genes...')
      message(runningMsg)
      log_print(runningMsg, console=F)
      fun <- eval(as.name(paste0('run', method)))
      colStr <- paste0(method, '_', setName)
      seuratObj <- silently_run(timeCode(fun, logPrint=TRUE, seuratObj, markerList[[i]], colStr))
    }
  }
  log_close(footer=F)
  return(seuratObj)
}
