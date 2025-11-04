#' @importFrom ggplot2 aes cut_number element_text geom_point ggplot ggtitle labs scale_color_manual theme theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices rainbow
#' @importFrom henna classPlot correlationPlot densityPlot rankSummary rankPlot
#' @importFrom reshape2 melt
#' @importFrom rlang .data
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
#' @param pointSize Point size.
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' scorePlot(smr[[1]][[1]])
#'
#' @export
#'
scorePlot <- function(scoreDF,
                      title = NULL,
                      xLabel = 'Score',
                      legendTitle = 'Gene set',
                      pointSize = 1.5){
    scoreDF <- scoreDF[order(scoreDF$avg), ]
    longDF <- reshape2::melt(as.matrix(scoreDF [, -ncol(scoreDF),
                                                drop=FALSE]))
    pal <- rainbow(length(colnames(scoreDF)))
    p <- ggplot(data=longDF) +
        geom_point(mapping=aes(x=.data[['value']],
                               y=.data[['Var1']],
                               color=.data[['Var2']]),
                   size=pointSize) +
        labs(x=xLabel, y ='Method', color=legendTitle, title=title) +
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
#' @param efBenchmark A list of dataframes generated with efficiencyBenchmark,
#' containing an element labeled 'time'.
#' @param titleSuffix Plot title suffix to be appended to the default title.
#' @param ... Additional parameters to be passed to \code{scorePlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' timePlot(smr[[5]])
#'
#' @export
#'
timePlot <- function(efBenchmark, titleSuffix = NULL, ...){
    title <- suffixedTitle('Running times', titleSuffix)
    p <- scorePlot(efBenchmark$time, title, 'Running time (s)', ...)
    return(p)
}

#' Plot a data frame consisting of gene set analysis method peak memory usage
#
#' This function plots data frame consisting of method peak memory usages with
#' methods as rows, gene sets and the gene set average as columns.
#'
#' @param efBenchmark A list of dataframes generated with efficiencyBenchmark,
#' containing an element labeled 'space'.
#' @inheritParams timePlot
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' memoryPlot(smr[[5]])
#'
#' @export
#'
memoryPlot <- function(efBenchmark, titleSuffix = NULL, ...){
    title <- suffixedTitle('Peak memory usage', titleSuffix)
    p <- scorePlot(efBenchmark$space, title, 'Peak memory usage (MiB)', ...)
    return(p)
}

#' Plot a list of data frame scores
#'
#' This function plots a list of dataframe scores with methods as rows,
#' gene sets and the average of scores across all gene sets as columns.
#'
#' @param smr List of summary data frames, whether boundary, MCC or global.
#' @inheritParams timePlot
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
benchmarkPlots <- function(smr, titleSuffix = NULL, ...){

    if ('sensitivity' %in% names(smr)){
        titleTemplates <- c('Sensitivity', 'Specificity', 'Precision',
                            'Accuracy', 'Size proximity','Score coverage',
                            'Class boundary determination gene set summary',
                            'Class boundary determination metric summary')
        prefix <- 'boundary_'
    } else {
        titleTemplates <- c('AUROC', 'PRAUC', 'Label rank alignment',
                            'Silhouette rank alignment', 'Centrality',
                            'Label Jaccard score', 'Label cosine score',
                            'Global evaluation gene set summary',
                            'Global evaluation metric summary')
        prefix <- 'global_'
    }

    names(titleTemplates) <- names(smr)

    if(!is.null(titleSuffix))
        titleTemplates <- setNames(paste0(titleTemplates, titleSuffix),
                                   names(titleTemplates))

    plots <- lapply(seq_len(length(smr)), function(i)
        scorePlot(smr[[i]], titleTemplates[names(smr)[i]], ...))

    plots[[length(plots)]] <- plots[[length(plots)]] + labs(color='Metric')
    names(plots) <- names(smr)
    slice <- c(length(plots) - 1, length(plots))
    names(plots)[slice] <- paste0(prefix, names(plots)[slice])
    return(plots)
}

#' Plot the complete list of benchmark summaries
#'
#' This function plots the complete list of benchmark summaries.
#'
#' @param smr Complete summary list generated with \code{allBenchmarkResults}.
#' @inheritParams timePlot
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
                                  paste0('MCC with boundary threshold',
                                         titleSuffix), ...),
               scorePlot(smr$MCC$direct,
                         paste0('Comprehensive MCC', titleSuffix),
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
geneSetRankPlots <- function(smr, titleSuffix = NULL, rankMethod = 'min', ...){
    geneSetNames <- colnames(smr$boundary[[1]])
    geneSetNames <- geneSetNames[seq(length(geneSetNames) - 1)]

    message('Computing gene set ranks...')
    gsRankDFs <- allGeneSetRanks(smr, rankMethod)


    plots <- mapply(function(gsRankDF, gsName){
        title <- suffixedTitle(paste0('Distribution of ranks for ',
                                      gsName, ' gene set'),
                               titleSuffix)
        return(rankPlot(gsRankDF, title, summarize=FALSE, xLab='Method',
                 ...))
    }, gsRankDFs, geneSetNames, SIMPLIFY=FALSE)

    names(plots) <- geneSetNames
    return(plots)
}

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
metricRankPlots <- function(smr, titleSuffix = NULL, rankMethod = 'min', ...){
    metricNames <- c('sensitivity', 'specificity', 'precision',
                     'accuracy', 'size proximity','score coverage',
                     'boundary MCC', 'direct MCC',
                     'AUROC', 'PRAUC', 'label rank alignment',
                     'silhouette rank alignment', 'centrality',
                     'label Jaccard score', 'label cosine score')

    message('Computing metric ranks...')
    metricRanksDFs <- allMetricRanks(smr, rankMethod)


    plots <- mapply(function(metricRankDF, metricName){
        title <- suffixedTitle(paste0('Distribution of ',
                                      metricName, ' ranks'),
                               titleSuffix)
        return(rankPlot(metricRankDF, title, summarize=FALSE, xLab='Method',
                        ...))
    }, metricRanksDFs, metricNames, SIMPLIFY=FALSE)

    names(plots) <- metricNames
    return(plots)
}

#' Create aggregate rank plot from a summary object
#'
#' This function creates an aggregate rank plot from a summary object.
#'
#' @inheritParams allBenchmarkPlots
#' @inheritParams geneSetRankPlots
#' @param sigDigits Number of significant digits used when displaying mean
#' ranks. If \code{NULL}, the mean ranks will not be displayed.
#' @inheritParams aggregateRanks
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' p <- aggregateRankPlot(smr)
#'
#' @export
#'
aggregateRankPlot <- function(smr,
                              titleSuffix = NULL,
                              sigDigits = 2,
                              rankMethod = c('min', 'average', 'first',
                                             'last', 'random', 'max'),
                              ...){
    message('Computing aggregate ranks...')
    aggRanks <- aggregateRanks(smr, rankMethod)
    title <- suffixedTitle('Distribution of aggregate ranks', titleSuffix)

    p <- rankPlot(aggRanks,
                  title,
                  summarize=FALSE,
                  xLab='Method',
                  sigDigits=sigDigits,
                  ...)
    return(p)
}

#' Create ratio rank plot for the method results
#'
#' This function creates a ratio rank plot for method results.
#'
#' @inheritParams allBenchmarkPlots
#' @inheritParams allTopRatios
#' @param xLab Label of the x axis.
#' @param yLab Label of the y axis.
#' @param legendLab Legend title.
#' @param ... Additional arguments passed to \code{henna::classRank}.
#'
#' @return A ggplot object.
#'
#' @examples
#' sPath <- system.file('extdata', 'smr.qs', package='GSABenchmark')
#' smr <- qs::qread(sPath)
#' p <- ratioPlot(smr)
#'
#' @export
#'
ratioPlot <- function(smr,
                      titleSuffix = NULL,
                      nItems = 25,
                      xLab = 'Maximum over mean ratio',
                      yLab = 'Metric',
                      legendLab = 'Method',
                      ...){
    message('Computing ratio ranks...')
    ratioDF <- allTopRatios(smr, nItems)
    title <- suffixedTitle('Top maximum over mean ratios', titleSuffix)
    p <- classPlot(ratioDF,
                   paste0(title, titleSuffix),
                   xLab,
                   yLab,
                   legendLab,
                   ...)
    return(p)
}

#' Create MDS plots for method results
#'
#' This function creates MDS plots for method results.
#'
#' @inheritParams mdsScoreSummary
#' @inheritParams timePlot
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
mdsPlots <- function(scObj, smr, titleSuffix = NULL, ...){
    message('Computing MDS for method results...')

    mdsDFs <- mdsScoreSummary(scObj, smr)
    plotNames <- names(mdsDFs)
    plots <- mapply(function(mdsDF, plotName){
        title <- suffixedTitle(paste0('MDS plot - ', plotName), titleSuffix)
        return(densityPlot(mdsDF, title, drawScores=TRUE, ...))
    }, mdsDFs, plotNames, SIMPLIFY=FALSE)
    return(plots)
}

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
corrPlots <- function(scObj, smr, titleSuffix = NULL, ...){
    message('Computing correlations for method results...')

    corrDFs <- corrSummary(scObj, smr)
    plotNames <- names(corrDFs)
    plots <- mapply(function(corrDF, plotName){
        title <- suffixedTitle(paste0('Correlation plot - ', plotName), titleSuffix)
        correlationPlot(corrDF, title, ...)
    }, corrDFs, plotNames, SIMPLIFY=FALSE)
    return(plots)
}
