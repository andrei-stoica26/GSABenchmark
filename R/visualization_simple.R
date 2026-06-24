#' @importFrom ggplot2 aes cut_number element_text geom_point ggplot ggtitle labs scale_color_manual theme theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom henna classPlot correlationPlot densityPlot rankSummary rankPlot tilePlot
#' @importFrom paletteer paletteer_c
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#'
NULL

#' Plot a data frame consisting of gene set analysis method scores
#
#' This function plots a data frame consisting of method scores with methods
#' as rows, gene sets and the gene set average as columns.
#'
#' @param scoreDF A summary data frame.
#' @param title Plot title.
#' @param xLab x axis label.
#' @param yLab y axis label.
#' @param isDecreasing Logical; whether the methods should be displayed on the
#' plot in decreasing order of the obtained average scores. If \code{FALSE}
#' (as default), the methods will be displayed in increasing order of the
#' average scores.
#' @param palette Color palette.
#' @param legendLab Legend label.
#' @param pointSize Point size.
#' @param pointShape Point shape.
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs2', package='GSABenchmark')
#' smr <- qs2::qs_read(sPath)
#' scorePlot(smr[[1]][[1]])
#'
#' @export
#'
scorePlot <- function(scoreDF,
                      title = NULL,
                      xLab = 'Score',
                      yLab = 'Method',
                      isDecreasing = FALSE,
                      palette = 'grDevices::Dark 3',
                      legendLab = 'Gene set',
                      pointSize = 1.5,
                      pointShape = 19){
    scoreDF <- scoreDF[order(scoreDF$avg, decreasing=isDecreasing), ]
    longDF <- reshape2::melt(as.matrix(scoreDF [, -ncol(scoreDF),
                                                drop=FALSE]))
    nColors <- length(colnames(scoreDF))
    pal <- paletteer_c(palette, nColors)
    p <- ggplot(data=longDF) +
        geom_point(mapping=aes(x=.data[['value']],
                               y=.data[['Var1']],
                               color=.data[['Var2']]),
                   size=pointSize,
                   shape=pointShape) +
        labs(x=xLab, y=yLab, color=legendLab, title=title) +
        theme_minimal() + scale_color_manual(values=pal,
                                             breaks=colnames(scoreDF)) +
        theme(plot.title = element_text(hjust = 0.5))
    return(p)
}

#' Plot a data frame consisting of gene set analysis method running times
#
#' This function plots data frame consisting of method running times with
#' methods as rows, gene sets and the gene set average as columns.
#'
#' @inheritParams scorePlot
#' @param efBenchmark A list of dataframes generated
#' with \code{efficiencyBenchmark}.
#' @param ... Additional parameters to be passed to \code{scorePlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs2', package='GSABenchmark')
#' smr <- qs2::qs_read(sPath)
#' timePlot(smr[[5]])
#'
#' @export
#'
timePlot <- function(efBenchmark,
                     title = NULL,
                     xLab = 'Running time (s)',
                     pointShape=17,
                     ...)
    return(scorePlot(efBenchmark$time,
                     title,
                     xLab,
                     isDecreasing=TRUE,
                     pointShape=pointShape,
                     ...))

#' Plot a data frame consisting of gene set analysis method peak memory usage
#
#' This function plots data frame consisting of method peak memory usages with
#' methods as rows, gene sets and the gene set average as columns.
#'
#' @inheritParams timePlot
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs2', package='GSABenchmark')
#' smr <- qs2::qs_read(sPath)
#' memoryPlot(smr[[5]])
#'
#' @export
#'
memoryPlot <- function(efBenchmark,
                       title = NULL,
                       xLab = 'Peak memory usage (MiB)',
                       pointShape=17,
                       ...)
    return(scorePlot(efBenchmark$space,
                     title,
                     xLab,
                     isDecreasing=TRUE,
                     pointShape=pointShape,
                     ...))

#' Create an aggregate rank plot from a summary object
#'
#' This function creates an aggregate rank plot from a summary object.
#'
#' @inheritParams scorePlot
#' @inheritParams allBenchmarkPlots
#' @inheritParams geneSetRankPlots
#' @param sigDigits Number of significant digits used when displaying mean
#' ranks. If \code{NULL}, the mean ranks will not be displayed.
#' @inheritParams aggregateRanks
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs2', package='GSABenchmark')
#' smr <- qs2::qs_read(sPath)
#' p <- aggregateRankPlot(smr)
#'
#' @export
#'
aggregateRankPlot <- function(smr,
                              title = NULL,
                              xLab = 'Method',
                              sigDigits = 2,
                              rankMethod = c('min', 'average', 'first',
                                             'last', 'random', 'max'),
                              ...){
    aggRanks <- aggregateRanks(smr, rankMethod)
    p <- rankPlot(aggRanks,
                  title,
                  summarize=FALSE,
                  xLab=xLab,
                  sigDigits=sigDigits,
                  ...)
    return(p)
}

#' Create a ratio rank plot for the method results
#'
#' This function creates a ratio rank plot for the method results.
#'
#' @inheritParams scorePlot
#' @inheritParams allBenchmarkPlots
#' @inheritParams allTopRatios
#' @param ... Additional arguments passed to \code{henna::classRank}.
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs2', package='GSABenchmark')
#' smr <- qs2::qs_read(sPath)
#' p <- ratioPlot(smr)
#'
#' @export
#'
ratioPlot <- function(smr,
                      title = NULL,
                      nItems = 25,
                      xLab = 'Maximum over mean ratio',
                      yLab = 'Metric',
                      legendLab = 'Method',
                      ...){
    ratioDF <- allTopRatios(smr, nItems)
    p <- classPlot(ratioDF,
                   title,
                   xLab,
                   yLab,
                   legendLab,
                   ...)
    return(p)
}
