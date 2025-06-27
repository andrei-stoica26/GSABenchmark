#' @importFrom Seurat AddModuleScore
#' @importFrom SiPSiC getPathwayScores
#' @importFrom VAM vam
NULL

#' Run AddModuleScore
#'
#' This function runs AddModuleScore
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to AddModuleScore
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runAddModuleScore <- function(seuratObj, genes, colStr = 'AddModuleScore', ...){
  seuratObj <- AddModuleScore(seuratObj, features = list(genes), name = colStr, ...)
  seuratObj@meta.data[[colStr]] <- liver::minmax(seuratObj@meta.data[[paste0(colStr, 1)]])
  seuratObj@meta.data[[paste0(colStr, 1)]] <- c()
  return(seuratObj)
}

#' Run SiPSiC
#'
#' This function runs SiPSiC
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to getPathwayScores
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runSiPSiC <- function(seuratObj, genes, colStr = 'SiPSiC', ...){
  mat <- LayerData(seuratObj, layer = 'counts')
  scores <- SiPSiC::getPathwayScores(mat, genes, ...)[[2]]
  seuratObj@meta.data[[colStr]] <- liver::minmax(scores)
  return(seuratObj)
}

#' Run VAM
#'
#' This function runs VAM
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to VAM
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runVAM <- function(seuratObj, genes, colStr = 'VAM', ...){
  seuratSubs <- subset(seuratObj, features = genes)
  mat <- t(as.matrix(LayerData(seuratSubs, layer = 'data')))
  v <- vam(mat, ...)
  seuratObj@meta.data[[colStr]] <- v$cdf.value
  return(seuratObj)
}
