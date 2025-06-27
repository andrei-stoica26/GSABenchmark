#' @importFrom decoupleR run_mdt run_ora run_udt
#' @importFrom SeuratObject LayerData
#'
NULL

#' Run a gene set analysis method using decoupleR
#'
#' This function runs one of the gene set analysis methods supported by decoupleR
#'
#' @param seuratObj A Seurat object
#' @param genes A vector of genes
#' @param method Gene set analysis method
#' @param colStr Name of the results column
#'
#' @return A Seurat object with the results saved as a metadata column
#'
runDecoupleRMethod <- function(seuratObj, genes, method, colStr = method){
  mat <- suppressWarnings(as.matrix(LayerData(seuratObj, layer = 'data')))
  scores <- do.call(paste0('run_', tolower(method)), list(mat, network = data.frame(source='geneSet', target=genes, mor=1)))$score
  seuratObj@meta.data[[colStr]] <- liver::minmax(scores)
  return(seuratObj)
}

#' Run MDT using decoupleR
#'
#' This function runs MDT using decoupleR
#'
#' @inheritParams runDecoupleRMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runMDT <- function(seuratObj, genes, colStr = 'MDT')
  return(runDecoupleRMethod(seuratObj, genes, 'MDT', colStr))

#' Run ORA using decoupleR
#'
#' This function runs ORA using decoupleR
#'
#' @inheritParams runDecoupleRMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runORA <- function(seuratObj, genes, colStr = 'ORA')
  return(runDecoupleRMethod(seuratObj, genes, 'ORA', colStr))

#' Run UDT using decoupleR
#'
#' This function runs UDT using decoupleR
#'
#' @inheritParams runDecoupleRMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runUDT <- function(seuratObj, genes, colStr = 'UDT')
  return(runDecoupleRMethod(seuratObj, genes, 'UDT', colStr))
