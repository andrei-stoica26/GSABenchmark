#' @importClassesFrom Seurat Seurat
#' @importFrom SeuratObject Embeddings
#'
NULL

#' Calculate the coordinates of the center of mass
#'
#' This function calculates the coordinates of the center of mass based on a
#' matrix of cell embeddings and a vector of weights
#'
#' @param dimMat A matrix of cell embeddings
#' @param weights A vector of weights
#'
#' @return A vector containing the coordinates of the center of mass
#'
#' @export
#'
centerOfMass <- function(dimMat, weights){
  totalWeight <- sum(weights)
  return(apply(dimMat, 2, function(x) sum(x * weights) / totalWeight))
}
