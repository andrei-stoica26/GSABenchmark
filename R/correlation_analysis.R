#' @importFrom stats cor dist setNames
#'
NULL


#' Calculate correlation matrix for method results
#'
#' This function calculates the correlation matrix for all the methods for an
#' identity class
#'
#' @inheritParams extractCellScores
#' @param gsaMethods Gene set analysis methods.
#' @param scoreCols Gene set analysis method score columns.
#' @param corMethod Correlation method.
#'
#' @return A correlation matrix
#'
#' @export
#'
methodsCor <- function(scObj, gsaMethods, scoreCols = gsaMethods, corMethod = 'pearson'){
  m <- as.matrix(metadataDF(scObj)[, scoreCols])
  corMat <- cor(m, method = corMethod)
  rownames(corMat) <- gsaMethods
  colnames(corMat) <- gsaMethods
  return(corMat)
}

#' Calculate correlation matrix for method results for multiple identity classes
#' and Seurat objects
#'
#' This function gets the mean correlation matrix for gene set analysis method
#' results for multiple identity classes and Seurat objects
#'
#' @inheritParams methodsCor
#' @param scoreColList List of equally-sized character vectors
#'
#' @return A correlation matrix
#'
#' @export
#'
allCor <- function(scObj, gsaMethods, scoreColList){
  corMats <- lapply(scoreColList, function(x) methodsCor(scObj, gsaMethods, x))
  return(Reduce(`+`, corMats) / length(corMats))
}
