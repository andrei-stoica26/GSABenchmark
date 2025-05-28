#' @importFrom qs qsave
#'
NULL

#' Run gene set analysis methods
#'
#' This function runs the gene set analysis methods.
#'
#' @param seuratObj A Seurat object
#' @param namedGeneSets A named list of character vectors
#' @param fileName Log file name
#' @param gsaMethods List of gene set analysis methods
#' @param joinChar Character used to create column names by joining the name of
#' each method with the name of each gene set
#'
#' @return A Seurat object with the results of the runs stored as metadata columns
#'
#' @export
#'
runMethods <- function(seuratObj, namedGeneSets, fileName, joinChar = '', gsaMethods = list('RunAddModuleScore', 'RunAUCell', 'RunCSOA',
                                                                        'RunMDT', 'RunPLAGE', 'RunSargent', 'RunSiPSiC',
                                                                        'RunSingscore', 'RunscSE', 'RunssGSEA', 'RunUCell',
                                                                        'RunVAM', 'RunVIPER', 'RunZScore')){
  if(!is.null(fileName))
    sink(fileName)
  tryCatch(
    {
      for (i in seq_len(length(namedGeneSets)))
        for (fun in gsaMethods){
          messageMethod(fun, names(namedGeneSets)[i])
          if (!is.null(fileName))
            messageMethod(fun, names(namedGeneSets)[i], print)
          seuratObj <- timeCode(eval(as.name(fun)), doMessage=TRUE,
                                seuratObj, namedGeneSets[[i]], paste0(str_remove(fun, 'Run'), names(namedGeneSets)[i]))
        }
      safeSink(fileName)
    },
    error = function(e){
      safeSink(fileName)
      stop(e)
    }
  )
  return(seuratObj)
}

#' Run gene set analysis methods on top cluster markers
#'
#' This function runs the gene set analysis methods on top cluster markers.
#'
#' @inheritParams prepareSeurat
#' @inheritParams MarkerList::allMarkers
#' @param ... Additional parameters to functions
#'
#' @return A Seurat object with the results of the runs stored as metadata columns
#'
#' @export
#'
clusterRun <- function(seuratObj, resolution, objectName, ...){
  seuratObj <- prepareSeurat(seuratObj, resolution)
  namedGeneSets <- prepareMarkers(seuratObj, objectName, ...)
  fileName <- paste0(objectName, 'Log.txt')
  seuratObj <- runMethods(seuratObj, namedGeneSets, fileName, ...)
  qsave(seuratObj, paste0(objectName, 'SeuratGSA.qs'))
  return(seuratObj)
}
