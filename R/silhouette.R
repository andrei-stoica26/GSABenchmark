#' @importFrom cluster silhouette
#' @importFrom liver minmax
#' @importFrom stylo dist.cosine
#' @importFrom stats cor dist weighted.mean
NULL

#' Compute cluster silhouette for Seurat object
#'
#' This function computes the silhouette for each cell in the Seurat object
#'
#' @param seuratObj A Seurat object
#' @param distMetric A distance metric. Must be one of "cosine", euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"
#'
#' @return A Seurat object with a silhouette column in the metadata
#'
#' @export
#'
computeSilhouette <- function(seuratObj, distMetric = 'cosine'){
  pcaMat <- as.matrix(Embeddings(seuratObj, reduction = "pca"))
  message('Computing distance matrix...')
  if (distMetric == 'cosine') distMat <- stylo::dist.cosine(pcaMat) else distMat <- stats::dist(x=pcaMat, method=distMetric)
  message('Computing silhouette...')
  clusters <- as.integer(seuratObj$seurat_clusters) - 1L
  seuratObj$silhouette <- cluster::silhouette(clusters, distMat)[, 3]
  return(seuratObj)
}

#' Normalize silhouette by cluster for Seurat object
#'
#' This function normalized the already computed silhouette for each cluster in
#' a Seurat object
#'
#' @param seuratObj A Seurat object with a silhouette metadata column
#'
#' @return A Seurat object with a normSilhouette column in the metadata
#'
#' @export
#'
normalizeSilhouette <- function(seuratObj){
  seuratObj$normSilhouette <- 0
  for (i in levels(seuratObj)){
    message(paste0('Normalizing silhouette: Cluster ', i, '...'))
    clusterSubset <- subset(seuratObj, seurat_clusters == i)
    seuratObj$normSilhouette[colnames(clusterSubset)] <- liver::minmax(clusterSubset$silhouette)
  }
  return(seuratObj)
}

#' Calculate silhouette-weighted mean for column
#'
#' This function computed the silhouette-weighted mean for the column in a
#' dataframe having normSilhouette as a column
#'
#' @param df Data frame
#' @param colStr Name of column
#'
#' @return A numeric vector
#'
#' @export
#'
silhouetteWeightedMean <- function(df, colStr)
  return(weighted.mean(df[, colStr], df$normSilhouette))
