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
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object
#' with the results of the runs stored as metadata columns.
#'
#' @export
#'
runGSAMethods <- function(scObj, labelCol, geneSets, gsaMethods, geneSetNames){
    checkSetNames(scObj, labelCol, geneSetNames)
    for (i in seq_along(geneSets)){
        setName <- geneSetNames[i]
        for (j in seq_along(gsaMethods)){
            method <- gsaMethods[j]
            message(paste0('Running ', method, ' for ', setName, ' genes...'))
            fun <- eval(as.name(paste0('run', method)))
            colStr <- paste0(method, '_', setName)
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
    return(c('AUCell', 'CSOA', 'GSVA', 'MDT', 'MLM', 'ORA', 'Pagoda2',
             'PLAGE', 'Singscore', 'SiPSiC', 'ssGSEA', 'UCell', 'UDT',
             'VAM', 'Zscore'))

