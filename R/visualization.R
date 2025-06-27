#' @importFrom ggplot2 aes geom_point ggtitle
#' @importFrom ggrepel geom_text_repel
#' @importFrom Seurat FeaturePlot
#' @importFrom SeuratObject Embeddings
NULL

#' Plot centers of mass of all methods over the feature plot of a method of
#' choice
#'
#' This function shows the locations of the centers of mass of all methods
#' superimposed over a feature plot showing a method of choice
#'
#' @param seuratObj A Seurat object
#' @param method The name of the method. Must exist as a metadata column in the
#' Seurat object
#' @param centersDF The data frame with the UMAP coordinates for the centers of
#' mass for each method
#' @param plotTitle Plot title
#'
#' @return A ggplot object
#'
#' @export
#'
cmPlot <- function(seuratObj, method, centersDF, plotTitle = method){
  p <- FeaturePlot(seuratObj, method, label=T, label.size=2.2, repel=T) +
    geom_point(aes(umap_1, umap_2), data=centersDF, shape=4, size=2) +
    geom_text_repel(aes(umap_1, umap_2), data=centersDF, label=rownames(centersDF), size=2.5, max.overlaps=30) + ggtitle(plotTitle)
  return(p)
}

#' Call cmPlot for all methods
#'
#' This function shows the locations of the centers of mass of all methods
#' superimposed over a feature plot showing a method of choice
#'
#' @param seuratObj A Seurat object
#' @param method The name of the method. It must have been run in the Seurat
#' object using clusterRun
#' @param cluster The cluster as character
#' @param centersList The list of data frame with the UMAP coordinates for the
#' centers of mass for each method, for each cluster where the methods were run
#' @param plotTitle Plot title
#' @param suffix A character object to append at the end of the plot title, separated with a dash. Ignored if plotTitle is not NULL
#'
#' @return A ggplot object
#'
#' @export
#'
cmPlotClusters <- function(seuratObj, method, cluster, centersList, plotTitle=NULL, suffix=NULL){
  if (!is.character(cluster))
    stop('cluster must be a character object')
  df <- centersList[cluster]
  if (is.null(plotTitle)){
    plotTitle <- paste0(method, ' in cluster ')
    if (!is.null(suffix))
      plotTitle <- paste0(plotTitle, ' - ', suffix)
  }
  p <- cmPlot(seuratObj, paste0(method, cluster), df, plotTitle)
  return(p)
}
