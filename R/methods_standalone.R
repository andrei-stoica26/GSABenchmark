#' @importFrom Matrix t
#' @importFrom pagoda2 score.cells.puram
#' @importFrom Seurat AddModuleScore
#' @importFrom singscore rankGenes simpleScore
#' @importFrom SiPSiC getPathwayScores
#' @importFrom VAM vam
#'
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

#' Run pagoda2
#'
#' This function runs pagoda2
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to pagoda2::score.cells.puram
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runPagoda2 <- function(seuratObj, genes, colStr = 'Pagoda2', ...){
  mat <- Matrix::t(LayerData(seuratObj, layer = 'data'))
  scores <- pagoda2::score.cells.puram(mat, genes, ...)
  seuratObj@meta.data[[colStr]] <- liver::minmax(scores)
  return(seuratObj)
}

#' Run Singscore
#'
#' This function runs Singscore
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to singscore::simpleScore
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runSingscore <- function(seuratObj, genes, colStr = 'Singscore', ...){
  mat <- as.matrix(LayerData(seuratObj, layer = 'data'))
  mat <- singscore::rankGenes(mat)
  scores <- singscore::simpleScore(mat, genes, ...)$TotalScore
  seuratObj@meta.data[[colStr]] <- liver::minmax(scores)
  return(seuratObj)
}

#' Run SiPSiC
#'
#' This function runs SiPSiC
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to SiPSiC::getPathwayScores
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
#' @param ... Additional arguments passed to VAM::vam
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runVAM <- function(seuratObj, genes, colStr = 'VAM', ...){
  mat <- t(as.matrix(LayerData(seuratObj, layer='data')[genes, ]))
  v <- vam(mat, ...)
  seuratObj@meta.data[[colStr]] <- v$cdf.value
  return(seuratObj)
}
