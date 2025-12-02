#' Plot a list of summary data frames
#'
#' This function plots a list of summary data frames with methods as rows,
#' gene sets and the average of scores across all gene sets as columns.
#'
#' @inheritParams plotList
#' @param smr List of summary data frames, whether boundary, MCC or global.
#' @param lastPlotLegendLab The legend label of the last plot.
#' @param ... Additional arguments passed to \code{scorePlot}.
#'
#' @return A list of ggplot objects.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' plots <- benchmarkPlots(smr[[1]])
#'
#' @export
#'
benchmarkPlots <- function(smr,
                           titleSuffix = NULL,
                           lastPlotLegendLab = 'Metric',
                           ...){

    if ('sensitivity' %in% names(smr)){
        titleInfixes <- c('Sensitivity', 'Specificity', 'Precision',
                            'Accuracy', 'Size proximity','Score coverage',
                            'Class boundary determination gene set summary',
                            'Class boundary determination metric summary')
        prefix <- 'boundary_'
    } else {
        titleInfixes <- c('AUROC', 'PRAUC', 'Label rank alignment',
                            'Silhouette rank alignment', 'Centrality',
                            'Label Jaccard score', 'Label cosine score',
                            'Global evaluation gene set summary',
                            'Global evaluation metric summary')
        prefix <- 'global_'
    }

    plots <- plotList(args=list(smr),
                      calcFun=identity,
                      plotFun=scorePlot,
                      unlistPlotArgs=FALSE,
                      titlePrefix=NULL,
                      titleInfixes=titleInfixes,
                      titleSuffix=titleSuffix,
                      ...)
    plots[[length(plots)]] <- plots[[length(plots)]] +
        labs(color=lastPlotLegendLab)
    slice <- c(length(plots) - 1, length(plots))
    names(plots)[slice] <- paste0(prefix, names(plots)[slice])
    return(plots)
}

#' Plot the complete list of benchmark summaries
#'
#' This function plots the complete list of benchmark summaries.
#'
#' @param smr Complete summary list generated with \code{allBenchmarkResults}
#' or \code{runBenchmark}.
#' @inheritParams benchmarkPlots
#' @param ... Additional parameters passed to other functions.
#'
#' @return A list of ggplot objects.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' plots <- allBenchmarkPlots(smr)
#'
#' @export
#'
allBenchmarkPlots <- function(smr, titleSuffix = NULL, ...){
    p1 <- benchmarkPlots(smr$boundary, titleSuffix, ...)
    p2 <- setNames(list(scorePlot(smr$MCC$boundary,
                                  paste0('MCC with boundary threshold - ',
                                         titleSuffix), ...),
                        scorePlot(smr$MCC$direct,
                                  paste0('Comprehensive MCC - ', titleSuffix),
                                  ...)), names(smr$MCC))
    p3 <- benchmarkPlots(smr$global, titleSuffix, ...)

    plots <- c(p1, p2, p3)
    if ('efficiency' %in% names(smr)){
        p4 <- setNames(list(timePlot(smr$efficiency, titleSuffix, ...),
                            memoryPlot(smr$efficiency, titleSuffix, ...)),
                       names(smr$efficiency))
        plots <- c(plots, p4)
    }
    return(plots)
}

#' Create gene set rank plots for the method results
#'
#' This function creates gene set rank plots for method results.
#'
#' @inheritParams allBenchmarkPlots
#' @inheritParams allGeneSetRanks
#' @inheritParams plotList
#' @inheritParams scorePlot
#'
#' @param ... Additional arguments passed to \code{henna::rankPlot}.
#'
#' @return A named list of ggplot objects.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' plots <- geneSetRankPlots(smr)
#'
#' @export
#'
geneSetRankPlots <- function(smr,
                             rankMethod = 'min',
                             titlePrefix = 'Distribution of ranks',
                             titleInfixes = NULL,
                             titleSuffix = NULL,
                             xLab = 'Method',
                             ...)
    return(plotList(args=list(smr, rankMethod),
                    calcFun=allGeneSetRanks,
                    plotFun=rankPlot,
                    unlistPlotArgs=FALSE,
                    titlePrefix=titlePrefix,
                    titleInfixes=titleInfixes,
                    titleSuffix=titleSuffix,
                    summarize=FALSE,
                    xLab=xLab,
                    ...))

#' Create metric rank plots for the method results
#'
#' This function creates metric rank plots for method results.
#'
#' @inheritParams geneSetRankPlots
#'
#' @return A named list of ggplot objects.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' plots <- metricRankPlots(smr)
#'
#' @export
#'
metricRankPlots <- function(smr,
                            titlePrefix = 'Distribution of ranks',
                            titleInfixes = metricNames(),
                            titleSuffix = NULL,
                            rankMethod = 'min',
                            xLab = 'Method',
                            ...)
    return(plotList(args=list(smr, rankMethod),
                    calcFun=allMetricRanks,
                    plotFun=rankPlot,
                    unlistPlotArgs=FALSE,
                    titlePrefix=titlePrefix,
                    titleInfixes=titleInfixes,
                    titleSuffix=titleSuffix,
                    summarize=FALSE,
                    xLab=xLab,
                    ...))

#' Create MDS plots for method results
#'
#' This function creates MDS plots for method results.
#'
#' @inheritParams henna::densityPlot
#' @inheritParams extractCellScores
#' @inheritParams allBenchmarkPlots
#' @inheritParams plotList
#'
#' @param ... Additional arguments passed to \code{henna::densityPlot}.
#'
#' @return A named list of ggplot objects.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' plots <- mdsPlots(scObj, smr)
#'
#' @export
#'
mdsPlots <- function(scObj,
                     smr,
                     titlePrefix = 'MDS plot',
                     titleInfixes = NULL,
                     titleSuffix = NULL,
                     drawScores = TRUE,
                     ...)
    return(plotList(args=list(scObj, smr),
                    calcFun=mdsScoreSummary,
                    plotFun=densityPlot,
                    unlistPlotArgs=TRUE,
                    titlePrefix=titlePrefix,
                    titleInfixes=titleInfixes,
                    titleSuffix=titleSuffix,
                    c(drawScores=drawScores, list(...))))

#' Create correlation plots for method results
#'
#' This function creates correlation plots for method results.
#'
#' @inheritParams mdsPlots
#' @param ... Additional arguments passed to \code{henna::correlationPlot}.
#'
#' @return A named list of ggplot objects.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' plots <- corrPlots(scObj, smr)
#'
#' @export
#'
corrPlots <- function(scObj,
                      smr,
                      titlePrefix = 'Correlation plot',
                      titleInfixes = NULL,
                      titleSuffix = NULL,
                      ...)
    return(plotList(args=list(scObj, smr),
                    calcFun=corrSummary,
                    plotFun=correlationPlot,
                    unlistPlotArgs=FALSE,
                    titlePrefix=titlePrefix,
                    titleInfixes=titleInfixes,
                    titleSuffix=titleSuffix,
                    ...))

#' Create Jaccard tile plots for method binary predictions
#'
#' This function creates Jaccard tile plots for method binary predictions.
#'
#' @inheritParams predJaccards
#' @inheritParams plotList
#' @param ... Additional arguments passed to \code{henna::tilePlot}.
#'
#' @return A named list of ggplot objects.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' plots <- predJaccardPlots(smr$predictions)
#'
#' @export
#'
predJaccardPlots <- function(predictionsSmr,
                             titlePrefix = 'Binary prediction Jaccard plot',
                             titleInfixes = NULL,
                             titleSuffix = NULL,
                             ...)
    return(plotList(args=list(predictionsSmr),
                    calcFun=predJaccards,
                    plotFun=tilePlot,
                    unlistPlotArgs=FALSE,
                    titlePrefix=titlePrefix,
                    titleInfixes=titleInfixes,
                    titleSuffix=titleSuffix,
                    ...))
