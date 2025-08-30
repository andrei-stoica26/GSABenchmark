#' @importFrom escape escape.matrix
#'
NULL

#' Run a gene set analysis method using escape
#'
#' This function runs one of the gene set analysis methods supported
#' by \code{escape}.
#'
#' @inheritParams runDecoupleRMethod
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
runEscapeMethod <- function(scObj, geneSets, method, ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    mat <- scExpMat(scObj, 'data')

    scoreDF <- escape.matrix(mat, geneSets, method=method, ...)
    scoreDF <- apply(scoreDF, 2, safeMinmax)

    scObj <- attachCellScores(scObj, scoreDF)
    return(scObj)
}

#' Run AUCell
#'
#' This function runs \code{AUCell} using \code{escape}.
#'
#' @inheritParams runEscapeMethod
#' @param ... Additional parameters passed to \code{runEscapeMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runAUCell <- function(scObj, geneSets, ...)
    return(runEscapeMethod(scObj, geneSets, 'AUCell', ...))

#' Run ssGSEA
#'
#' This function runs \code{ssGSEA} using \code{escape}.
#'
#' @inheritParams runAUCell
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runssGSEA <- function(scObj, geneSets, ...)
    return(runEscapeMethod(scObj, geneSets, 'ssGSEA', ...))

#' Run UCell
#'
#' This function runs \code{UCell} using \code{escape}.
#'
#' @inheritParams runAUCell
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runUCell <- function(scObj, geneSets, ...)
    return(runEscapeMethod(scObj, geneSets, 'UCell', ...))
