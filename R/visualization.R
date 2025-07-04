#' @importFrom ggplot2 aes geom_point ggplot ggtitle labs scale_color_manual theme theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices rainbow
#' @importFrom reshape2 melt
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

#' Plot a data frame score
#'
#' This function plots a dataframe score with methods as rows, gene sets and the
#' average of scores across all gene sets as columns.
#'
#' @param scoreDF A score data frame
#' @param title Plot title
#'
#' @return A ggplot object
#'
#' @export
#'
scorePlot <- function(scoreDF, title){
  scoreDF <- scoreDF[order(scoreDF$avg), ]
  longDF <- reshape2::melt(as.matrix(scoreDF [, -ncol(scoreDF)]))
  pal <- rainbow(length(rownames(scoreDF )))
  p <- ggplot(data=longDF) +
    geom_point(mapping=aes(x=value, y=Var1, color=Var1)) +
    labs(y = 'Score', x = 'Method', color = 'Method', title = title) +
    theme_minimal() + scale_color_manual(values=pal, breaks = rev(rownames(scoreDF))) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}

#' Plot a data frame score
#'
#' This function plots a dataframe score with methods as rows, gene sets and the
#' average of scores across all gene sets as columns.
#'
#' @param smr Summary list
#' @param datasetName Dataset name
#'
#' @return A list of ggplot objects
#'
#' @export
#'
benchmarkPlots <- function(smr, datasetName = NULL){
  v <- c('Sensitivity', 'Specificity', 'Precision', 'Accuracy', 'Size proximity',
         'Score specificity', 'Silhouette coverage', 'Centrality', 'AUROC',
         'Gini coefficient', 'Kolgomorov-Smirnov statistics', 'PRAUC')
  if (!length(intersect(names(smr), c('AUC', 'Gini', 'KS_Stat', 'PRAUC'))))
    v <- c(v, c('Class boundary benchmark gene set summary',
           'Class boundary benchmark metric summary')) else
             v <- c(v, c('Distribution benchmark gene set summary',
                    'Distribution benchmark metric summary'))

  if(!is.null(datasetName))
    v <- paste0(v, ' - ', datasetName)
  names(v) <- c('sensitivity', 'specificity', 'precision', 'accuracy',
                'sizeProximity','scoreSpecificity', 'silhouetteCoverage',
                'centrality', 'AUC','Gini', 'KS_Stat', 'PRAUC',
                'avg', 'metricSummary')
  plots <- lapply(seq_len(length(smr)), function(i) scorePlot(smr[[i]], v[names(smr)[i]]))
  return(plots)
}
