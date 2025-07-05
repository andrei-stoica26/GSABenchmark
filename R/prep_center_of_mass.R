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

#' Calculate the coordinates of the center of mass at each position
#'
#' This function calculates the coordinates of the centers of mass obtained based
#' on a matrix of cell embeddings and a vector of weights, for each cutoff.
#'
#' @param dimMat A matrix of cell embeddings
#' @param weights A vector of weights
#'
#' @return A matrix of centers of mass
#'
#' @export
#'
centerOfMassV <- function(dimMat, weights)
  return(apply(dimMat * weights, 2, cumsum) / cumsum(weights))
