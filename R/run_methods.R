#' @importFrom fabR silently_run
#' @importFrom CSOA runCSOA
#'
NULL

#' Run gene set analysis methods
#'
#' This function runs the gene set analysis methods.
#'
#' @inheritParams extractCellScores
#' @param geneSets A list of gene sets.
#' @param gsaMethods Character vector of gene set analysis methods.
#' @param geneSetNames The names of the gene sets.
#' @param infix Infix to add between method name and gene set name in
#' single-cell expression object. The string consisting of the method name and
#' the infix is separated by the gene set name with a '_' character.
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object
#' with the results of the runs stored as metadata columns.
#'
#' @export
#'
runGSAMethods <- function(scObj, labelCol, geneSets, gsaMethods, geneSetNames,
                          infix = NULL){
    checkSetNames(scObj, labelCol, geneSetNames)
    if (length(geneSets) != length(geneSetNames))
        stop('geneSets and geneSetNames must have the same length.')
    for (i in seq_along(geneSets)){
        setName <- geneSetNames[i]
        for (j in seq_along(gsaMethods)){
            method <- gsaMethods[j]
            message(paste0('Running ', method, ' for ', setName, ' genes...'))
            fun <- eval(as.name(paste0('run', method)))
            colStr <- paste0(method,
                             infix, '_',
                             setName)
            scObj <- silently_run(fun(scObj, geneSets[[i]], colStr))
        }
    }
    return(scObj)
}

#' Show supported methods
#'
#' This function shows methods currently supported by GSABenchmark
#'
#' @return A character vector of supported methods
#'
#' @export
#'
supportedMethods <- function()
    return(c('AUCell', 'CSOA', 'GSVA', 'JASMINE', 'MDT', 'MLM', 'ORA',
             'Pagoda2', 'PLAGE', 'Singscore', 'SiPSiC', 'ssGSEA', 'UCell',
             'UDT', 'VAM', 'Zscore'))
?Seurat::SCTransform
