#' @importFrom fabR silently_run
#' @importFrom CSOA attachCellScores runCSOA
#'
NULL

#' Run gene set analysis methods
#'
#' This function runs the gene set analysis methods.
#'
#' @inheritParams extractCellScores
#' @param geneSets A named list of gene sets.
#' @param gsaMethods Character vector of gene set analysis methods.
#' @param infix Infix to add between method name and gene set name in
#' single-cell expression object. The string consisting of the method name and
#' the infix is separated by the gene set name with a '_' character.
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object
#' with the results of the runs stored as metadata columns.
#'
#' @examples
#' scObj <- qs::qread('inst/testdata/scObj.qs')
#' geneSets <- qs::qread('inst/testdata/geneSets.qs')
#' scObj <- runGSAMethods(scObj, 'label', geneSets, c('CSOA', 'Zscore'))
#'
#' @export
#'
runGSAMethods <- function(scObj, labelCol, geneSets, gsaMethods, infix = NULL){
    geneSetNames <- names(geneSets)
    checkSetNames(scObj, labelCol, geneSetNames)
    for (method in gsaMethods){
        message(paste0('Running ', method, '...'))
        fun <- eval(as.name(paste0('run', method)))
        names(geneSets) <- paste0(method, infix, '_', geneSetNames)
        scObj <- silently_run(fun(scObj, geneSets))
    }
    return(scObj)
}

#' Show supported methods
#'
#' This function shows methods currently supported by GSABenchmark.
#'
#' @return A character vector of supported methods.
#'
#' @examples
#' supportedMethods()
#'
#' @export
#'
supportedMethods <- function()
    return(c('AddModuleScore', 'AUCell', 'CSOA', 'GSVA', 'JASMINE', 'MDT',
             'MLM', 'ORA', 'Pagoda2', 'PLAGE', 'Singscore', 'SiPSiC',
             'ssGSEA', 'UCell', 'UDT', 'VAM', 'Zscore'))
