#' @importFrom escape escape.matrix
#'
NULL

#' Run a gene set analysis method using escape
#'
#' This function runs one of the gene set analysis methods supported by escape
#'
#' @inheritParams runDecoupleRMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
runEscapeMethod <- function(seuratObj, genes, method, colStr = method, ...){
  mat <- suppressWarnings(as.matrix(LayerData(seuratObj, layer = 'data')))
  scores <- escape.matrix(mat, list(set1 = genes), method = method, ...)[, 1]
  seuratObj@meta.data[[colStr]] <- safeMinmax(scores)
  return(seuratObj)
}

#' Run AUCell
#'
#' This function runs AUCell using escape
#'
#' @inheritParams runEscapeMethod
#' @param ... Additional parameters passed to runEscapeMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runAUCell <- function(seuratObj, genes, colStr = 'AUCell', ...)
  return(runEscapeMethod(seuratObj, genes, 'AUCell', colStr, ...))

#' Run ssGSEA
#'
#' This function runs ssGSEA using escape
#'
#' @inheritParams runAUCell
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runssGSEA <- function(seuratObj, genes, colStr = 'ssGSEA', ...)
  return(runEscapeMethod(seuratObj, genes, 'ssGSEA', colStr, ...))

#' Run UCell
#'
#' This function runs UCell using escape
#'
#' @inheritParams runAUCell
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runUCell <- function(seuratObj, genes, colStr = 'UCell', ...)
  return(runEscapeMethod(seuratObj, genes, 'UCell', colStr, ...))
