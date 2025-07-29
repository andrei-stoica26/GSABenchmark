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
runEscapeMethod <- function(scObj, genes, method, colStr = method, ...){
    mat <- scExpMat(scObj, 'data')
    scores <- escape.matrix(mat, list(set1 = genes), method = method, ...)[, 1]
    scObj[[colStr]] <- safeMinmax(scores)
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
runAUCell <- function(scObj, genes, colStr = 'AUCell', ...)
    return(runEscapeMethod(scObj, genes, 'AUCell', colStr, ...))

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
runssGSEA <- function(scObj, genes, colStr = 'ssGSEA', ...)
    return(runEscapeMethod(scObj, genes, 'ssGSEA', colStr, ...))

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
runUCell <- function(scObj, genes, colStr = 'UCell', ...)
    return(runEscapeMethod(scObj, genes, 'UCell', colStr, ...))
