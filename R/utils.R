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
#' @noRd
#'
messageMethod <- function(funStr, group, messageFun = message)
    return(messageFun(paste0(str_replace(funStr,
                                         'Run', 'Running '), ' for identity: ',
                             group, '...')))

#' Calculate the alignment between two numeric vectors of the same length
#'
#' This function calculates the alignment between two numeric vectors of the
#' same lenght.
#'
#' @param v A vector.
#' @param w A vector of the same length as
#'
#' @return Alignment score.
#'
#' @examples
#' v <- c(2, 3, 6, 7, 8, 4, 12, 9, 10)
#' w <- c(3, 4, 5, 6, 2, 7, 8, 13, 3)
#' alignmentScore(v, w)
#'
#' @export
#'
alignmentScore <- function(v, w){
    if (length(v) != length(w))
        stop('`v` and `w` must have the same length.')
    return(sum(v * w) / sum(sort(v) * sort(w)))
}

#' Calculate the rank alignment between two numeric vectors of the same length
#'
#' This function calculates the rank alignment between two numeric vectors of
#' the same length.
#'
#' @inheritParams alignmentScore
#'
#' @return Rank alignment score.
#'
#' @examples
#' v <- c(2, 3, 6, 7, 8, 4, 12, 9, 10)
#' w <- c(3, 4, 5, 6, 2, 7, 8, 13, 3)
#' rankAlignmentScore(v, w)
#'
#' @export
#'
rankAlignmentScore <- function(v, w){
    if (length(v) != length(w))
        stop('`v` and `w` must have the same length.')
    v <- safeMinmax(rank(v), 1 / length(v))
    w <- safeMinmax(rank(w), 1 / length(w))
    return(alignmentScore(v, w))
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

#' Return metric names
#'
#' This function returns metric names.
#'
#' @return A vector of characters.
#'
#' @noRd
#'
metricNames <- function(){
    metrics <- c('sensitivity', 'specificity', 'precision',
                 'accuracy', 'size proximity','score coverage',
                 'boundary MCC', 'direct MCC',
                 'AUROC', 'PRAUC', 'label rank alignment',
                 'silhouette rank alignment', 'centrality',
                 'label Jaccard score', 'label cosine score')
    return(metrics)
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
                 'AUROC', 'PRAUC', 'labRankAlignment', 'silRankAlignment',
                 'centrality', 'labJaccard', 'labCosine')
    titles <- c('Sensitivity', 'Specificity', 'Precision',
                'Accuracy', 'Size proximity','Score coverage',
                'Boundary MCC', 'Direct MCC',
                'AUROC', 'PRAUC', 'Label rank alignment',
                'Silhouette rank alignment', 'Centrality',
                'Label Jaccard score', 'Label cosine score')
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

#' Remove single-valued columns from a data frame
#'
#' This function removes single-valued columns from a data frame.
#'
#' @param df A data frame.
#'
#' @return A data frame from which single-valued columns have been removed.
#'
#' @noRd
#'
removeSVCols <- function(df){
    hasDistinctValues <- apply(df, 2, function(x) length(unique(x)) > 1)
    removedCol <- names(which(!hasDistinctValues))
    if(length(removedCol))
        message('Single-valued columns will be removed: ', removedCol, '.')
    df <- df[, hasDistinctValues]
    return(df)
}
