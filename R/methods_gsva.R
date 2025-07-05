#' @importFrom GSVA gsva gsvaParam plageParam zscoreParam
#'
NULL

#' Run a gene set analysis method using GSVA
#'
#' This function runs one of the gene set analysis methods supported by escape
#'
#' @inheritParams runDecoupleRMethod
#' @param invert Whether to transform the scores from x to 1 - x
#' @param filter Whether to filter the expression matrix as to contain only
#' signature genes
#'
#' @return A Seurat object with the results saved as a metadata column
#'
runGSVAMethod <- function(seuratObj, genes, method, colStr = method,
                          invert = FALSE, filter = FALSE, ...){
  if(filter)
    mat <- as.matrix(LayerData(seuratObj, layer='data')[genes, ]) else
      mat <- as.matrix(LayerData(seuratObj, layer='data'))

  geneSet <- setNames(list(genes), 'sigScore')
  gsvaPar <- do.call(paste0(tolower(method), 'Param'), list(mat, geneSet, ...))
  scores <- gsva(gsvaPar)[1, ]
  scores <- liver::minmax(scores)
  if (invert)
    scores <- 1 - scores
  seuratObj@meta.data[[colStr]] <- scores
  return(seuratObj)
}

#' Run GSVA
#'
#' This function runs GSVA using GSVA
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to runGSVAMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runGSVA <- function(seuratObj, genes, colStr = 'GSVA', ...)
  return(runGSVAMethod(seuratObj, genes, 'GSVA', colStr, ...))

#' Run PLAGE
#'
#' This function runs PLAGE using GSVA
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to runGSVAMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runPLAGE <- function(seuratObj, genes, colStr = 'PLAGE', ...)
  return(runGSVAMethod(seuratObj, genes, 'PLAGE', colStr, TRUE, TRUE, ...))

#' Run Zscore
#'
#' This function runs Zscore using GSVA
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to runGSVAMethod
#'
#' @return A Seurat object with the results saved as a metadata column
#'
#' @export
#'
runZscore <- function(seuratObj, genes, colStr = 'Zscore', ...)
  return(runGSVAMethod(seuratObj, genes, 'Zscore', colStr, FALSE, TRUE, ...))
