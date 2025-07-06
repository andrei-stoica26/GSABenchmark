#' @importFrom fabR silently_run
#'
NULL


#' Run gene set analysis methods
#'
#' This function runs the gene set analysis methods.
#'
#' @param seuratObj A Seurat object.
#' @param geneSets A list of gene sets.
#' @param geneSetNames The names of the gene sets.
#' @param gsaMethods Character vector of gene set analysis methods.
#'
#' @return A Seurat object with the results of the runs stored as metadata columns.
#'
#' @export
#'
runGSAMethods <- function(seuratObj, geneSets, geneSetNames, gsaMethods){
  for (i in seq_along(geneSets)){
    setName <- geneSetNames[i]
    for (j in seq_along(gsaMethods)){
      method <- gsaMethods[j]
      message(paste0('Running ', method, ' for ', setName, ' genes...'))
      fun <- eval(as.name(paste0('run', method)))
      colStr <- paste0(method, '_', setName)
      seuratObj <- silently_run(fun(seuratObj, geneSets[[i]], colStr))
    }
  }
  return(seuratObj)
}

#' Show supported methods
#'
#' This function shows methods currently supported by GSABenchmark
#'
#' @return A character vector of supported methods
#'
#' @export
#'
supportedMethods <- function()
  return(c('AddModuleScore', 'AUCell', 'GSVA', 'MDT', 'MLM', 'ORA', 'Pagoda2', 'PLAGE', 'Singscore', 'SiPSiC', 'ssGSEA', 'UCell',
           'UDT', 'VAM', 'Zscore'))

#' Perform min-max normalization when possible; otherwise return a zero-vector
#'
#' This function min-max-normalizes a vector when possible, and otherwise returns
#' the zero vector
#'
#' @param scores Numeric vector
#'
#' @return Min-max-normalized scores or the zero vector
#'
#' @export
#'
safeMinmax <- function(scores){
  if(length(unique(scores)) < 2)
    return(rep(0, length(scores)))
  return(liver::minmax(scores))
}
