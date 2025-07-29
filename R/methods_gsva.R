#' @importFrom GSVA gsva gsvaParam plageParam zscoreParam
#'
NULL

#' Run a gene set analysis method using \code{GSVA}
#'
#' This function runs one of the gene set analysis methods supported
#' by \code{GSVA}.
#'
#' @inheritParams runDecoupleRMethod
#' @param invert Whether to transform the scores from x to 1 - x.
#' @param filter Whether to filter the expression matrix as to contain only
#' signature genes.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
runGSVAMethod <- function(scObj, genes, method, colStr = method,
                          invert = FALSE, filter = FALSE, ...){
    if(filter)
        mat <- scExpMat(scObj, 'data', genes) else
            mat <- scExpMat(scObj, 'data')

    geneSet <- setNames(list(genes), 'sigScore')
    gsvaPar <- do.call(paste0(tolower(method), 'Param'), list(mat, geneSet, ...))
    scores <- gsva(gsvaPar)[1, ]
    scores <- safeMinmax(scores)
    if (invert)
        scores <- 1 - scores
    scObj[[colStr]] <- scores
    return(scObj)
}

#' Run GSVA
#'
#' This function runs \code{GSVA} using \code{GSVA}.
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to \code{runGSVAMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runGSVA <- function(scObj, genes, colStr = 'GSVA', ...)
    return(runGSVAMethod(scObj, genes, 'GSVA', colStr, ...))

#' Run PLAGE
#'
#' This function runs \code{PLAGE} using \code{GSVA}.
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to \code{runGSVAMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runPLAGE <- function(scObj, genes, colStr = 'PLAGE', ...)
    return(runGSVAMethod(scObj, genes, 'PLAGE', colStr, TRUE, TRUE, ...))

#' Run Zscore
#'
#' This function runs \code{Zscore} using \code{GSVA}.
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to \code{runGSVAMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runZscore <- function(scObj, genes, colStr = 'Zscore', ...)
    return(runGSVAMethod(scObj, genes, 'Zscore', colStr, FALSE, TRUE, ...))
