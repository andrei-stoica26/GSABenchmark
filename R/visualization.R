#' @importFrom ggplot2 aes element_text geom_point ggplot ggtitle labs scale_color_manual theme theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices rainbow
#' @importFrom reshape2 melt
#' @importFrom Seurat FeaturePlot
#' @importFrom SeuratObject Embeddings
NULL

#' Plot a data frame consisting of gene set analysis method scores
#
#' This function plots data frame consisting of method scores with methods
#' as rows, gene sets and the gene set average as columns.
#'
#' @param scoreDF A score data frame.
#' @param title Plot title.
#' @param xLabel x axis label.
#'
#' @return A ggplot object.
#'
#' @export
#'
scorePlot <- function(scoreDF, title, xLabel = 'Score'){
  scoreDF <- scoreDF[order(scoreDF$avg), ]
  longDF <- reshape2::melt(as.matrix(scoreDF [, -ncol(scoreDF)]))
  pal <- rainbow(length(rownames(scoreDF )))
  p <- ggplot(data=longDF) +
    geom_point(mapping=aes(x=value, y=Var1, color=Var1)) +
    labs(x = xLabel, y = 'Method', color = 'Method', title = title) +
    theme_minimal() + scale_color_manual(values=pal, breaks = rev(rownames(scoreDF))) +
    theme(plot.title = element_text(hjust = 0.5))
  return(p)
}


#' Plot a data frame consisting of gene set analysis method running times
#
#' This function plots data frame consisting of method running times with methods
#' as rows, gene sets and the gene set average as columns.
#'
#' @param efBenchmark A list of dataframes generated with efficiencyBenchmark, containing an
#' element labeled "time".
#' @param datasetName Dataset name
#'
#' @return A ggplot object.
#'
#' @export
#'
timePlot <- function(efBenchmark, datasetName = NULL){
  title <- 'Running times'
  if (!is.null(dataset))
    title <- paste0(title, ' - ', datasetName, ' dataset')
  p <- scorePlot(efBenchmark$time, title, 'Running time (s)')
  return(p)
}

#' Plot a data frame consisting of gene set analysis method peak memory usage
#
#' This function plots data frame consisting of method peak memory usages with
#' methods as rows, gene sets and the gene set average as columns.
#'
#' @param efBenchmark A list of dataframes generated with efficiencyBenchmark, containing an
#' element labeled "space".
#' @param inheritParams timePlot
#'
#' @return A ggplot object.
#'
#' @export
#'
memoryPlot <- function(efBenchmark, datasetName = NULL){
  title <- 'Peak memory usage'
  if (!is.null(dataset))
    title <- paste0(title, ' - ', datasetName, ' dataset')
  p <- scorePlot(efBenchmark$space, title, 'Peak memory usage (MiB)')
  return(p)
}

#' Plot a list of data frame scores
#'
#' This function plots a list of dataframe scores with methods as rows, gene sets and the
#' average of scores across all gene sets as columns.
#'
#' @param smr Summary list
#' @inheritParams timePlot
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
    v <- paste0(v, ' - ', datasetName, ' dataset')
  names(v) <- c('sensitivity', 'specificity', 'precision', 'accuracy',
                'sizeProximity','scoreSpecificity', 'silhouetteCoverage',
                'centrality', 'AUC','Gini', 'KS_Stat', 'PRAUC',
                'avg', 'metricSummary')
  plots <- lapply(seq_len(length(smr)), function(i) scorePlot(smr[[i]], v[names(smr)[i]]))
  return(plots)
}
