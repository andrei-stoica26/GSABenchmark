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


#' Calculate the coordinates of the silhouette center of mass
#'
#' This function calculates the coordinates of the silhouette-weight center of
#' mass of a given Seurat cluster.
#'
#' @param seuratObj A Seurat object
#' @param cluster Cluster
#'
#' @return A vector containing the coordinates of the center of mass
#'
#' @export
#'
silhouetteCM <- function(seuratObj, cluster){
  seuratSub <- subset(seuratObj, seurat_clusters == cluster)
  return(centerOfMass(Embeddings(seuratSub, reduction="umap"), seuratSub$normSilhouette))
}

#' Calculate the coordinates of the center of mass
#'
#' This function calculates the in-cluster centers of mass using each of the
#' gene set analysis method scores, as well as cluster silhouette
#'
#' @param seuratObj A Seurat object
#' @param clusters Character vector of Seurat clusters
#' @param gsaMethods Character vector of gene set analysis methods
#'
#' @return A vector containing the coordinates of the center of mass
#'
#' @export
#'
centersAnalysis <- function(seuratObj, clusters, gsaMethods){
  dimMat <- Embeddings(seuratObj, reduction="umap")
  res <- lapply(clusters, function(x) {
    df <- data.frame(t(sapply(gsaMethods, function(y){
      weights <- seuratObj@meta.data[, paste0(y, x)]
      res <- centerOfMass(dimMat, weights)
    })))
    df <- rbind(df, silhouetteCM(seuratObj, x))
    rownames(df) <- c(gsaMethods, 'normSilhouette')
    return(df)
  })
  names(res) <- clusters
  return(res)
}
