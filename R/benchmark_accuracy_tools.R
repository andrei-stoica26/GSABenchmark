#' @importFrom hammers centerOfMass checkGenes computeSilhouette metadataDF metadataNames nearestNeighbors normalizeSilhouette numCosine proximity safeMinmax scCol scPCAMat scExpMat shuffleGenes tabulateVector timeMemoryExpr
#' @importFrom plyr count
#'
NULL

#' Extract gene set analysis method scores and truth labels from single-cell
#' expression object
#'
#' This function extracts gene set analysis method scores and truth labels from
#' a single-cell expression object.
#'
#' @param scObj A \code{Seurat} or \code{SingleCellExperiment} object.
#' @param labelCol The metadata column containing the ground truth annotation.
#' @param scoreCol The metadata column containing the gene set analysis method
#' score.
#' @param label The identity assessed from the labelCol column.
#'
#' @return A data frame with two columns: truth labels (1 or 0) and
#' gene set analysis method scores.
#'
#' @keywords internal
#'
extractCellScores <- function(scObj, labelCol, scoreCol, label){
    df <- metadataDF(scObj)[, c(labelCol, scoreCol)]
    colnames(df)[1] <- 'label'
    df$label <- as.integer(df[, 1] %in% label)
    df <- df[order(df[, scoreCol], decreasing=TRUE), ]
    return(df)
}

#' Condense repeated gene set analysis method scores
#'
#' This function condenses repeated gene set analysis method scores, summing
#' labels and recording frequencies in the process.
#'
#' @param df A data frame with two columns: truth labels (1 or 0) and gene set
#' analysis method scores.
#'
#' @return A condensed scores data frame.
#'
#' @noRd
#'
condenseRepeatedScores <- function(df){
    denseDF <- plyr::count(df[, 2])
    denseDF$label <- vapply(denseDF[, 1], function(score)
        sum(df[, 1][which(df[, 2] == score)]), numeric(1))
    denseDF <- denseDF[, c(3, 1, 2)]
    denseDF <- denseDF[order(denseDF[, 2], decreasing=TRUE), ]
    colnames(denseDF)[2] <- c(colnames(df)[2])
    return(denseDF)
}

#' Add overall scores to a benchmark data frame
#'
#' This function adds overlap scores to a benchmark data frame.
#'
#' @param df Benchmark data frame.
#' @param startCol Column where metric scores start.
#'
#' @return The input data frame sorted decreasingly by the newly added
#' overall column.
#'
#' @keywords internal
#'
computeMetricMeans <- function(df, startCol){
    df$avg <- rowMeans(df[, seq(startCol, ncol(df))])
    df <- df[order(df$avg, decreasing=TRUE),]
    return(df)
}

#' Add method information to a summary data frame
#'
#' This function adds gene set analysis method information to a summary data
#' frame. It sets the names of the gene set analysis method as the rownames of
#' the data frame, and optionally computes the mean score for each method and
#' decreasingly sorts the data frame by these scores.
#'
#' @inheritParams runGSAMethods
#' @param df A data frame where the values represent the scores obtained by a
#' gene set analysis method (row) on a gene set (column) for a metric.
#' @param doAverage Whether to add an average column to each data frame,
#' sorting each data frame decreasingly by the average column in the process.
#'
#' @return A summary data frame with added method information.
#'
#' @keywords internal
#'
addMethodInfo <- function(df, gsaMethods, doAverage=TRUE){
    rownames(df) <- gsaMethods
    if(doAverage){
        df$avg <- rowMeans(df)
        df <- df[order(df$avg, decreasing=TRUE), ]
    }
    return(df)
}

#' Add metric information to the summary list
#'
#' This function adds metric information to the summary list and optionally
#' extends the list by addding overall results for each metric.
#'
#' @param smr List of result data frames for each metrics. Each data frame
#' contains the results for each tested gene set analysis method for each gene
#' set for the corresponding method.
#' @param metrics Metrics.
#' @inheritParams addMethodInfo
#' @param doSummarize Whether to add a metric summary. Must be set to
#' \code{FALSE} when summarizing a MCC benchmark list of lists.
#'
#' @return Extended summary list with an additional data frame showing the
#' average results obtained for each metric.
#'
#' @keywords internal
#'
addMetricSummary <- function(smr, metrics, gsaMethods, doSummarize=TRUE){
    names(smr) <- metrics
    if (doSummarize){
        df <- data.frame(lapply(smr, function(x) x[gsaMethods, ]$avg))
        rownames(df) <- gsaMethods
        df <- df[order(df$avg, decreasing=TRUE), ]
        smr <- c(smr, list(metricSummary=df))
    }
    return(smr)
}

#' Add gene set averages to the global summary list
#'
#' This function adds gene set averages to the global summary list.
#'
#' @inheritParams runGSAMethods
#' @param globalSmr Global summary.
#'
#' @return Extended summary list with an additional data frame showing the
#' average results obtained for each gene set.
#'
#' @keywords internal
#'
addGlobalAverage <- function(globalSmr, gsaMethods){
    smrNames <- names(globalSmr[[1]])
    df <- data.frame(do.call(cbind, setNames(lapply(smrNames, function(smrCol){
        rowMeans(do.call(cbind, lapply(globalSmr,
                                       function(x) x[gsaMethods, smrCol])))
    }), smrNames)))
    rownames(df) <- gsaMethods
    df <- df[order(df$avg, decreasing=TRUE), ]
    globalSmr$avg <- df
    return(globalSmr)
}
