#' @importFrom methods is
#' @importFrom stats cmdscale cor dist rnorm setNames
#'
NULL

#' Message information on gene set analysis method run
#'
#' This function messages or prints information on the gene set analysis method
#' running and the group where it is running.
#'
#' @param funStr A supported gene set analysis method.
#' @param group A Seurat identity class.
#' @param messageFun An output function.
#'
#' @return Messaged (default) or printed output.
#'
messageMethod <- function(funStr, group, messageFun = message)
    return(messageFun(paste0(str_replace(funStr,
                                         'Run', 'Running '), ' for identity: ',
                             group, '...')))

#' Calculate the alignment between two numeric vectors
#'
#' This function calculates the alignment between two numeric vectors.
#'
#' @param x A vector.
#' @param y A vector
#'
#' @return Alignment score.
#'
#' @export
#'
alignmentScore <- function(x, y)
    return(sum(x * y) / sum(sort(x) * sort(y)))

#' Calculate the rank alignment between two numeric vectors
#'
#' This function calculates the rank alignment between two numeric vectors.
#'
#' @param x A vector.
#' @param y A vector.
#'
#' @return Alignment score.
#'
#' @export
#'
rankAlignmentScore <- function(x, y){
    x <- safeMinmax(rank(x), 1 / length(x))
    y <- safeMinmax(rank(y), 1 / length(y))
    return(alignmentScore(x, y))
}

#' Extract label from column name
#'
#' This function extracts the class identity label from column name.
#'
#' @param colStr Column name.
#'
#' @return Class identity label.
#'
#' @noRd
#'
labelFromColumn <- function(colStr){
    splitRun <- str_split(colStr, '_')[[1]]
    res <- splitRun[length(splitRun)]
    return(res)
}

#' Extract run name from column name
#'
#' This function extracts the run name from column name.
#'
#' @param colStr Column name.
#'
#' @return Run name.
#'
#' @noRd
#'
runFromColumn <- function(colStr){
    splitRun <- str_split(colStr, '_')[[1]]
    res <- paste0(splitRun[seq(length(splitRun) - 1)], collapse='_')
    return(res)
}

#' Convert metric names to titles
#'
#' This function converts metric names to titles.
#'
#' @return A named vector of characters.
#'
#' @noRd
#'
metricTitles <- function(){
    metrics <- c('sensitivity', 'specificity', 'precision',
                 'accuracy', 'sizeProximity','scoreCoverage',
                 'boundaryMCC', 'directMCC',
                 'AUROC', 'PRAUC', 'labRankAlignment',
                 'silRankAlignment', 'centrality')
    titles <- c('Sensitivity', 'Specificity', 'Precision',
                'Accuracy', 'Size proximity','Score coverage',
                'Boundary MCC', 'Direct MCC',
                'AUROC', 'PRAUC', 'Label rank alignment',
                'Silhouette rank alignment', 'Centrality')
    return(setNames(titles, metrics))
}

#' Attach gene set analysis method scores
#'
#' This function attaches gene set analysis method scores to an object.
#'
#' @details Wrapper around \code{CSOA::attachCellScores}
#'
#' @inheritParams runGSAMethods
#' @param scoreDF Data frame of gene set analysis method scores.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @noRd
#'
attachScores <- function(scObj, geneSets, scoreDF){
    rownames(scoreDF) <- colnames(scObj)
    colnames(scoreDF) <- names(geneSets)
    scObj <- attachCellScores(scObj, scoreDF)
    return(scObj)
}
