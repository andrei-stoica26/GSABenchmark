#' @importFrom abdiv cosine_distance
#' @importFrom cluster silhouette
#' @importFrom liver minmax
#' @importFrom stats cor dist weighted.mean
NULL

#' Compute cluster silhouette for Seurat object
#'
#' This function computes the silhouette for each cell in the Seurat object.
#'
#' @param seuratObj A Seurat object.
#' @param idClass Seurat identity class.
#' @param distMetric A distance metric. Must be one of "cosine", euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'
#' @return A Seurat object with a metadata silhouette column.
#'
#' @export
#'
computeSilhouette <- function(seuratObj, idClass = 'seurat_clusters', distMetric = 'cosine'){
  if (!idClass %in% colnames(seuratObj@meta.data))
    stop(paste0(idClass, ' not found in the metadata of the Seurat object'))
  if (!distMetric %in% c('cosine', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'))
    stop('Unknown distance metric. See ?GSABenchmark::computeSilhouette for supported metrics')
  pcaMat <- as.matrix(Embeddings(seuratObj, reduction = "pca"))
  message('Computing distance matrix...')
  if (distMetric == 'cosine') distMat <- cosine_distance(pcaMat) else
    distMat <- stats::dist(x=pcaMat, method=distMetric)
  message(paste0('Computing silhouette for identity class: ', idClass, '...'))
  groupVals <- unclass(factor(seuratObj@meta.data[[idClass]]))
  seuratObj$silhouette <- cluster::silhouette(groupVals, distMat)[, 3]
  return(seuratObj)
}

#' Normalize silhouette by identity class for Seurat object
#'
#' This function normalizes the already computed silhouette for each identity
#' class in a Seurat object
#'
#' @inheritParams computeSilhouette
#'
#' @return A data frame with normalized silhouettes for each unique element in the
#' identity class
#'
#' @export
#'
normalizeSilhouette <- function(seuratObj, idClass){
  df <- seuratObj@meta.data[, c(idClass, 'silhouette')]
  colnames(df)[1] <- 'label'
  groups <- unique(df[[idClass]])
  res <- data.frame(matrix(0, nrow(df), length(groups)))
  rownames(res) <- rownames(df)
  colnames(res) <- groups
  for (group in groups){
    dfSub <- subset(df, label == group)
    sortedSil <- dfSub$silhouette
    auxMin <- 2 * sortedSil[1] - sortedSil[2]
    normSilVals <- liver::minmax(c(dfSub$silhouette, auxMin))[seq_along(dfSub$silhouette)]
    res[rownames(dfSub), group] <- normSilVals
  }
  return(res)
}
