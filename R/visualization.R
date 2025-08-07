#' @importFrom ggplot2 aes cut_number element_text geom_point ggplot ggtitle labs scale_color_manual theme theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices rainbow
#' @importFrom henna densityPlot rankSummary rankPlot
#' @importFrom reshape2 melt
#'
NULL

#' Plot a data frame consisting of gene set analysis method scores
#
#' This function plots data frame consisting of method scores with methods
#' as rows, gene sets and the gene set average as columns.
#'
#' @param scoreDF A score data frame.
#' @param title Plot title.
#' @param xLabel x axis label.
#' @param legendTitle Legend title.
#'
#' @return A ggplot object.
#'
#' @export
#'
scorePlot <- function(scoreDF, title, xLabel = 'Score', legendTitle = 'Gene set'){
    scoreDF <- scoreDF[order(scoreDF$avg), ]
    longDF <- reshape2::melt(as.matrix(scoreDF [, -ncol(scoreDF),
                                                drop=FALSE]))
    pal <- rainbow(length(colnames(scoreDF)))
    p <- ggplot(data=longDF) +
        geom_point(mapping=aes(x=value, y=Var1, color=Var2), size=2.5) +
        labs(x=xLabel, y ='Method', color=legendTitle, title=title) +
        theme_minimal() + scale_color_manual(values=pal, breaks=colnames(scoreDF)) +
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
#' @param titleSuffix Plot title suffix to be appended to the default title
#'
#' @return A ggplot object.
#'
#' @export
#'
timePlot <- function(efBenchmark, titleSuffix = NULL){
    title <- paste0('Running times', titleSuffix)
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
#' @inheritParams timePlot
#'
#' @return A ggplot object.
#'
#' @export
#'
memoryPlot <- function(efBenchmark, titleSuffix = NULL){
    title <- paste0('Peak memory usage', titleSuffix)
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
benchmarkPlots <- function(smr, titleSuffix = NULL){
    boundaryMetrics <- c('Sensitivity', 'Specificity', 'Precision',
                         'Accuracy', 'Size proximity','Score coverage')
    globalMetrics <- c('AUROC', 'PRAUC', 'Label rank alignment',
                       'Silhouette rank alignment', 'Centrality')

    boundarySummaries <- c('Class boundary determination gene set summary',
                           'Class boundary determination metric summary')
    globalSummaries <- c('Global evaluation gene set summary',
                         'Global evaluation metric summary')

    if ('sensitivity' %in% names(smr)){
        v <- c(boundaryMetrics, boundarySummaries)
        names(v) <- c('sensitivity', 'specificity', 'precision', 'accuracy',
                      'sizeProximity','scoreCoverage', 'avg', 'metricSummary')
    } else{
        v <- c(globalMetrics, globalSummaries)
        names(v) <- c('AUROC', 'PRAUC', 'labRankAlignment',
                      'silRankAlignment', 'centrality',
                      'avg', 'metricSummary')
    }

    if(!is.null(titleSuffix))
        v <- setNames(paste0(v, titleSuffix), names(v))

    plots <- lapply(seq_len(length(smr)), function(i) scorePlot(smr[[i]], v[names(smr)[i]]))
    plots[[length(plots)]] <- plots[[length(plots)]] + labs(color='Metric')
    return(plots)
}

#' Plot the complete list of benchmark summaries
#'
#' This function plots the complete list of benchmark summaries.
#'
#' @param smr Complete summary list generated with allBenchmarkResults.
#' @inheritParams timePlot
#' @param ... Additional arguments passed to \code{henna::densityPlot}.
#'
#' @return A list of ggplot objects
#'
#' @export
#'
allBenchmarkPlots <- function(smr, titleSuffix = NULL, ...){
    p1 <- benchmarkPlots(smr$boundary, titleSuffix)
    p2 <- list(scorePlot(smr$MCC$boundary,
                         paste0('MCC with boundary threshold', titleSuffix)),
               scorePlot(smr$MCC$direct,
                         paste0('Comprehensive MCC', titleSuffix)))
    p3 <- benchmarkPlots(smr$global, titleSuffix)
    p4 <- mapply(function(mdsDF, gsName)
        densityPlot(mdsDF, paste0('MDS plot - ', gsName, ' genes',
                                  titleSuffix), ...),
        smr$MDS, names(smr$MDS), SIMPLIFY=FALSE)

    geneSetNames <- colnames(smr$boundary[[1]])
    geneSetNames <- geneSetNames[seq(length(geneSetNames) - 1)]

    metricNames <- c('sensitivity', 'specificity', 'precision',
    'accuracy', 'size proximity','score coverage',
    'MCC with boundary threshold', 'comprehensive MCC',
    'AUROC', 'PRAUC', 'label rank alignment',
    'silhouette rank alignment', 'centrality')

    p5 <- mapply(function(gsRankDF, gsName)
        rankPlot(gsRankDF,
                 paste0('Distribution of ranks for ',
                        gsName, ' gene set',
                        titleSuffix), summarize=FALSE, xLab='Method'),
        smr$gsRanks, geneSetNames, SIMPLIFY=FALSE)

    p6 <- mapply(function(metricRankDF, metricName)
                 rankPlot(metricRankDF,
                          paste0('Distribution of ',
                                 metricName, ' ranks',
                                 titleSuffix), summarize=FALSE, xLab='Method'),
                 smr$metricRanks, metricNames, SIMPLIFY=FALSE)

    p7 <- list(rankPlot(smr$aggRanks,
                   paste0('Distribution of aggregate ranks',
                          titleSuffix),
                   summarize=FALSE,
                   xLab='Method',
                   showMeanRanks=TRUE))

    plots <- c(p1, p2, p3, p4, p5, p6, p7)
    if ('efficiency' %in% names(smr)){
        p8 <- list(timePlot(smr$efficiency, titleSuffix),
                   memoryPlot(smr$efficiency, titleSuffix))
        plots <- c(plots, p8)
    }
    return(plots)
}
