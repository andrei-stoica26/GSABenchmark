#' @importFrom CSOA expMat
#' @importFrom escape escape.matrix
#' @importFrom irGSEA irGSEA.score
NULL

#' Runs a gene set analysis method using escape
#'
#' This function runs one of the gene set analysis methods supported by escape
#'
#' @param seuratObj A Seurat object
#' @param genes A vector of genes
#' @param method Gene set analysis method
#' @param colStr Name of the results column
#'
#' @return A Seurat object with the results saved as a metadata column
#'
runEscape <- function(seuratObj, genes, method, colStr = method){
  expObj <- CSOA::expMat(seuratObj)
  scores <- escape.matrix(expObj, list(set1 = genes), method = method)[, 1]
  seuratObj@meta.data[[colStr]] <- liver::minmax(scores)
  return(seuratObj)
}


#' Runs a gene set analysis method using irGSEA
#'
#' This function runs one of the gene set analysis methods supported by irGSEA
#'
#' @inheritParams runEscape
#' @param slot Seurat object slot
#' @param invert Return the complement of the original scores (1 - x)
#' @param nCores Number of cores
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#'
runIrGSEA <- function(seuratObj, genes, method, colStr = method, slot = 'data', invert = FALSE, nCores = 1){
  scores <- irGSEA.score(seuratObj, geneset = list(genes), method = method, slot = slot, custom = TRUE, minGSSize = 5,
                         maxGSSize = 1000, ncores = nCores)[[method]]@data['geneset1', ]
  if (invert)
    scores <- 1 - scores
  seuratObj@meta.data[[colStr]] <- liver::minmax(scores)
  return(seuratObj)
}
