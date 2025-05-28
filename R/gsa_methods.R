#' @importFrom SiPSiC getPathwayScores
#' @importFrom SeuratObject LayerData
#'
NULL

#' Run AddModuleScore using irGSEA
#'
#' This function runs AddModuleScore using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runAddModuleScore <- function(seuratObj, genes, colStr = 'AddModuleScore')
  return(runIrGSEA(seuratObj, genes, 'AddModuleScore', colStr))

#' Run AUCell using irGSEA
#'
#' This function runs AUCell using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runAUCell <- function(seuratObj, genes, colStr = 'AUCell')
  return(runIrGSEA(seuratObj, genes, 'AUCell', colStr))

#' Run MDT using irGSEA
#'
#' This function runs MDT using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runMDT <- function(seuratObj, genes, colStr = 'MDT')
  return(runIrGSEA(seuratObj, genes, 'mdt', colStr))

#' Run PLAGE using irGSEA
#'
#' This function runs PLAGE using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runPLAGE <- function(seuratObj, genes, colStr = 'PLAGE')
  return(runIrGSEA(seuratObj, genes, 'plage', colStr, invert = T))

#' Run Sargent using irGSEA
#'
#' This function runs Sargent using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runSargent <- function(seuratObj, genes, colStr = 'Sargent')
  return(runIrGSEA(seuratObj, genes, 'Sargent', colStr))

#' Run singscore using irGSEA
#'
#' This function runs singscore using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runSingscore <- function(seuratObj, genes, colStr = 'singscore')
  return(runIrGSEA(seuratObj, genes, 'singscore', colStr))

#' Run scSE using irGSEA
#'
#' This function runs scSE using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runscSE <- function(seuratObj, genes, colStr = 'scSE')
  return(runIrGSEA(seuratObj, genes, 'scSE', colStr))

#' Run UCell using irGSEA
#'
#' This function runs UCell using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runUCell <- function(seuratObj, genes, colStr = 'UCell')
  return(runIrGSEA(seuratObj, genes, 'UCell', colStr))

#' Run VAM using irGSEA
#'
#' This function runs VAM using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runVAM <- function(seuratObj, genes, colStr = 'VAM')
  return(runIrGSEA(seuratObj, genes, 'VAM', colStr))

#' Run VIPER using irGSEA
#'
#' This function runs VIPER using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runVIPER <- function(seuratObj, genes, colStr = 'VIPER')
  return(runIrGSEA(seuratObj, genes, 'viper', colStr))

#' Run Zscore using irGSEA
#'
#' This function runs Zscore using irGSEA
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runZScore <- function(seuratObj, genes, colStr = 'Zscore')
  return(runIrGSEA(seuratObj, genes, 'zscore', colStr))


#' Run SiPSiC
#'
#' This function runs SiPSiC
#'
#' @inheritParams runIrGSEA
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runSiPSiC <- function(seuratObj, genes, colStr = 'SiPSiC', slot = 'counts'){
  expMat <- LayerData(seuratObj, layer = slot)
  scores <- SiPSiC::getPathwayScores(expMat, genes)[[2]]
  seuratObj@meta.data[[colStr]] <- liver::minmax(scores)
  return(seuratObj)
}

#' Run ssGSEA
#'
#' This function runs ssGSEA using escape
#'
#' @inheritParams runEscape
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runssGSEA <- function(seuratObj, genes, colStr = 'test')
  return(runEscape(seuratObj, genes, 'ssGSEA', colStr))
