
#' Calculate correlation matrix for method results
#'
#' This function calculates the correlation matrix for all the methods for an
#' identity class
#'
#' @inheritParams runMethods
#' @param group Identity class
#' @param corMethod Correlation method
#'
#' @return A correlation matrix
#'
#' @export
#'
methodsCor <- function(seuratObj, gsaMethods, joinChar = '', group = NULL, corMethod = 'pearson'){
  m <- as.matrix(seuratObj@meta.data[, paste0(gsaMethods, joinChar, group)])
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
#' @param seurats List of Seurat objets
#' @inheritParams clusterBenchmark
#'
#' @return A correlation matrix
#'
#' @export
#'
allCor <- function(seurats, gsaMethods, groups){
  corMats <- mapply(function(x, y){
    return(lapply(y, function(z) methodsCor(x, gsaMethods, z)))
  }, seurats, groups)
  corMats <- unlist(corMats, recursive = F)
  return(Reduce(plus, corMats) / length(corMats))
}
