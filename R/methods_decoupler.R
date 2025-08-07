#' @importFrom decoupleR run_mdt run_mlm run_ora run_udt run_ulm
#'
NULL

#' Run a gene set analysis method using decoupleR
#'
#' This function runs one of the gene set analysis methods supported
#' by \code{decoupleR}.
#'
#' @inheritParams extractCellScores
#' @param genes A vector of genes
#' @param method Gene set analysis method
#' @param colStr Name of the results column
#' @param ... Additional arguments passed to gene set analysis methiod
#'
#' @return A single-cell expression object with the results saved as a metadata column.
#'
runDecoupleRMethod <- function(scObj, genes, method,
                               colStr = method, ...){
    checkGenes(scObj, genes)
    mat <- scExpMat(scObj, 'data')
    scores <- do.call(paste0('run_', tolower(method)),
                      list(mat, network=data.frame(source='geneSet',
                                                   target=genes,
                                                   mor=1,
                                                   ...)))$score
    names(scores) <- colnames(mat)
    scObj[[colStr]] <- safeMinmax(scores)
    return(scObj)
}

#' Run MDT using decoupleR
#'
#' This function runs \code{MDT} using \code{decoupleR}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional parameters passed to \code{runDecoupleRMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runMDT <- function(scObj, genes, colStr = 'MDT', ...)
    return(runDecoupleRMethod(scObj, genes, 'MDT', colStr, ...))

#' Run MLM using decoupleR
#'
#' This function runs \code{MLM} using \code{decoupleR}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional parameters passed to \code{runDecoupleRMethod}.
#'
#' @return A single-cell expression object with the results saved as a
#' metadata column.
#'
#' @export
#'
runMLM <- function(scObj, genes, colStr = 'MLM', ...)
    return(runDecoupleRMethod(scObj, genes, 'MLM', colStr, ...))

#' Run ORA using decoupleR
#'
#' This function runs \code{ORA} using \code{decoupleR}.
#'
#' @inheritParams runMDT
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runORA <- function(scObj, genes, colStr = 'ORA', ...)
    return(runDecoupleRMethod(scObj, genes, 'ORA', colStr, ...))

#' Run UDT using decoupleR
#'
#' This function runs \code{UDT} using \code{decoupleR}.
#'
#' @inheritParams runMDT
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runUDT <- function(scObj, genes, colStr = 'UDT', ...)
    return(runDecoupleRMethod(scObj, genes, 'UDT', colStr, ...))

