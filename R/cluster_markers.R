#' @importFrom Seurat FindClusters
#' @importFrom MarkerList allMarkers markerListNames
#'
NULL

#' Prepare Seurat object for gene set analysis using cluster markers
#'
#' This function finds clusters at a given resolution and computes the silhouette
#' and the cluster-normalized silhouette.
#'
#' @param seuratObj A Seurat object
#' @param resolution Clustering resolution
#'
#' @return A Seurat object with a silhouette column in the metadata
#'
#' @export
#'
prepareSeurat <- function(seuratObj, resolution){
  seuratObj <- FindClusters(seuratObj, resolution = resolution)
  seuratObj <- computeSilhouette(seuratObj)
  seuratObj <- normalizeSilhouette(seuratObj)
  return(seuratObj)
}

#' Prepare gene sets for gene set analysis using cluster markers
#'
#' This function constructs a list of gene set using the top markers from each
#' cluster
#'
#' @inheritParams MarkerList::allMarkers
#' @param ... Additional parameters to functions
#'
prepareMarkers <- function(seuratObj, objectName, ...){
  allClusterMarkers <- allMarkers(seuratObj, objectName, ...)
  return(markerListNames(allClusterMarkers, ...))
}


