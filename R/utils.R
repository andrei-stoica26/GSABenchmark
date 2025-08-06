#' @importFrom methods is
#' @importFrom stats cmdscale cor dist setNames
#'
NULL

#' Message information on gene set analysis method run
#'
#' This function messages or prints information on the gene set analysis method
#' running and the group where it is running.
#'
#' @param funStr A supported gene set analysis method
#' @param group A Seurat identity class
#' @param messageFun An output function
#'
#' @return Messaged (default) or printed output
#'
messageMethod <- function(funStr, group, messageFun = message)
    return(messageFun(paste0(str_replace(funStr,
                                         'Run', 'Running '), ' for identity: ',
                             group, '...')))

#' Take the square root of a non-negative number or the additive inverse of the
#' square root of the additive inverse of a negative number
#'
#' This function takes the square root of a non-negative number or the additive inverse of the
#' square root of the additive inverse of a negative number
#'
#' @param x A numeric variable
#'
#' @return The square root of x if positive and of -x if negative
#'
#' @export
#'
sqrt2 <- function(x)
    if (x >= 0) return(sqrt(x)) else return(-sqrt(-x))

#' Convert a vector to a data frame based on input row and column names
#'
#' This function converts a vector to a data frame based on input row and column names.
#' Optionally, it also calculates the row means.
#'
#' @param v A vector.
#' @param rowNames A character vector.
#' @param colNames A character vector.
#' @param addRowMeans Whether to add the row means to the data frame
#' @param sortByRowMeans Whether to sort by row means
#'
#' @return A data frame
#'
#' @export
#'
tabulateVector <- function(v,
                           rowNames,
                           colNames,
                           addRowMeans=FALSE,
                           sortByRowMeans=FALSE){
    df <- data.frame(matrix(v, length(rowNames), length(colNames)))
    rownames(df) <- rowNames
    colnames(df) <- colNames
    if(addRowMeans){
        df$avg <- rowMeans(df)
        if(sortByRowMeans)
        df <- df[order(df$avg), ]
    }
    return(df)
}

#' Calculate the alignment between two numeric vectors
#'
#' This function calculates the alignment between two numeric vectors
#'
#' @param x A vector.
#' @param y A vector
#'
#' @return Alignment score
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
