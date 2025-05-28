#' Calculate cluster benchmark
#'
#' This function calculates one of three supported benchmarks after the gene set
#' analysis methods have been run with cluster markers as input.
#'
#' @inheritParams runMethods
#' @param fun Cluster benchmark function. One of clusterScorePerc,
#' clusterSilhouetteMean and outclusterLowPerc
#' @param groups Vector of identity classes
#' @param groupStr Type of identity class
#'
#' @return List of two dataframes: raw benchmark scores and minmax-normalized ones
#'
#' @export
#'
clusterBenchmark <- function(seuratObj, fun, groups, gsaMethods,  groupStr = 'Cluster'){
  scores <- c()
  for (group in groups)
    for (method in gsaMethods)
      scores <- c(scores, fun(seuratObj, group, paste0(method, group)))
  df <- data.frame(pracma::Reshape(scores, length(gsaMethods), length(groups)))
  rownames(df) <- gsaMethods
  colnames(df) <- paste0(groupStr, groups)
  return(list(df, apply(df, 2, liver::minmax)))
}

#' Find center-of-mass distances between method scores and clusters
#'
#' This function construct a distance matrix and a proximity matrix based on the
#' method scores for each method and set of cluster markers and the normalized
#' silhouette score in each cluster.
#'
#' @param centersList A list of data frames generated with centersAnalysis. The
#' list contains the coordinates of the centers of mass for each method score
#' calculated by using markers of (selected) Seurat clusters as input gene sets,
#' as well as the coordinates of normalized silhouette-based cluster centers of
#' mass.
#'
#' @return A list containing a distance matrix and a proximity matrix, each having
#' method names as rows and clusters as solumns
#'
#' @export
#'
centerDistances <- function(centersList){
  distDF <- data.frame(sapply(centersList, function(x) apply(x[-nrow(x), ], 1, function(y) nnspat::euc.dist(x['normSilhouette', ], y))))
  proxDF <- data.frame(apply(distDF, 2, function(x) liver::minmax(min(x) / x)))
  colnames(distDF) <- str_replace(colnames(distDF), 'X', 'Cluster')
  colnames(proxDF) <- colnames(distDF)
  return(list(distDF, proxDF))
}


#' Calculate four cluster benchmarks
#'
#' This function calculates four benchmarks after the gene set analysis methods
#' have been run with cluster markers as input.
#'
#' @inheritParams clusterBenchmark
#' @param centersDF List of data frames produced by centerDistances
#'
#' @return List of four dataframes containing minmax-normalized ones
#'
#' @export
#'
accuracyScores <- function(seuratObj, groups, gsaMethods, centersDF){
  return(list(clusterBenchmark(seuratObj, clusterScorePerc, clusters, gsaMethods)[[2]],
              clusterBenchmark(seuratObj, clusterSilhouetteMean, clusters, gsaMethods)[[2]],
              clusterBenchmark(seuratObj, outclusterLowPerc, clusters, gsaMethods)[[2]],
              centerDistances(centersDF)[[2]]))
}
