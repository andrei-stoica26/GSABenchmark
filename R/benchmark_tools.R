#' Calculate the percentage of in-cluster score
#'
#' This function calculates the percentage of score attained by a method inside
#' a cluster
#'
#' @param seuratObj A Seurat object
#' @param cluster The cluster whose top markers were used as input for the method
#' @param colStr The Seurat metadata column storing the method results
#' @param verbose Verbose
#'
#' @return Percentage value
#'
#' @export
#'
clusterScorePerc <- function(seuratObj, cluster, colStr, verbose = FALSE){
  df <- seuratObj@meta.data[, c('seurat_clusters', colStr)]
  clusterScore <- sum(subset(df, seurat_clusters == cluster)[, 2])
  totalScore <- sum(df[, 2])
  perc <- clusterScore / totalScore * 100
  if(verbose)
    message(paste0('Percentage of ', colStr, ' score captured by Cluster ', cluster, ': ', perc, '%.'))
  return(perc)
}

#' Calculate the percentage of high scores in cluster
#'
#' This function calculates the percentage of cells with a higher than average
#' method score in cluster. If invert is set to TRUE, it calculates the percentage
#' of cells with a lower than average method score outside of cluster.
#'
#' @inheritParams clusterScorePerc
#' @param invert Assess high scores in cluster (FALSE, default) or low scores
#' outside of cluster (TRUE)
#'
#' @return Percentage value
#'
#' @export
#'
clusterHighPerc <- function(seuratObj, cluster, colStr, invert = FALSE, verbose = FALSE){
  cutoff <- mean(seuratObj@meta.data[, colStr])
  if (invert){
    df <- subset(seuratObj@meta.data, seurat_clusters != cluster)
    topCount <- length(which(df[, colStr] < cutoff))
    prep <- 'outside'
    comp <- 'lower'
  } else {
    df <- subset(seuratObj@meta.data, seurat_clusters == cluster)
    topCount <- length(which(df[, colStr] > cutoff))
    prep <- 'in'
    comp <- 'higher'
  }
  bottomCount <- nrow(df)
  perc <- topCount / bottomCount * 100
  if(verbose)
    message(paste0('Percentage of cells ', prep, ' Cluster ', cluster, ' with a ', comp, ' than average ', colStr,
                  ' score: ', perc, '%.'))
  return(perc)
}

#' Calculate the percentage of low scores outside of cluster
#'
#' This function calculates the percentage of cells with a lower than average
#' method score outside of cluster.
#'
#' @inheritParams clusterScorePerc
#'
#' @return Percentage value
#'
#' @export
#'
outclusterLowPerc <- function(seuratObj, cluster, colStr, verbose = FALSE)
  return(clusterHighPerc(seuratObj, cluster, colStr, invert = T, verbose))


#' Calculate the ratio of silhouette-weighted mean over mean
#'
#' This function calculate the ratio of silhouette-weighted mean over mean for a
#' method inside a cluster
#'
#' @inheritParams clusterScorePerc
#'
#' @return A ratio
#'
#' @export
#'
clusterSilhouetteMean <- function(seuratObj, cluster, colStr, verbose = FALSE){
  df <- subset(seuratObj@meta.data, seurat_clusters == cluster)
  res <- silhouetteWeightedMean(df, colStr) / mean(df[, colStr])
  if(verbose)
    message(paste0('Ratio between silhouette-weighted mean ', colStr,
                  ' score and mean ', colStr, ' score in Cluster ', cluster, ': ', res, '.'))
  return(res)
}
