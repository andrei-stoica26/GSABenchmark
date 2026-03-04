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
#' @param outputFun Choose between silently_run (suppress all warnings and
#' messages) or identity (do not suppress them). Default is silently_run.
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object
#' with the results of the runs stored as metadata columns.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runGSAMethods(scObj, 'label', geneSets, c('CSOA', 'Zscore'))
#'
#' @export
#'
runGSAMethods <- function(scObj,
                          labelCol,
                          geneSets,
                          gsaMethods,
                          infix = NULL,
                          outputFun = silently_run){
    geneSetNames <- names(geneSets)
    checkSetNames(scObj, labelCol, geneSetNames)
    slots <- expSlots(scObj)
    for (method in gsaMethods){
        message('Running ', method, '...')
        fun <- eval(as.name(paste0('run', method)))
        names(geneSets) <- paste0(method, infix, '_', geneSetNames)
        if (method == 'CSOA')
            scObj <- outputFun(fun(scObj, geneSets)) else
                scObj <- outputFun(fun(scObj, geneSets, slots[method]))
    }
    return(scObj)
}

#' Show supported methods
#'
#' This function shows the methods currently supported by GSABenchmark.
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

#' Show the default slot for each supported method
#'
#' This function shows the default slot for each supported method.
#'
#' @inheritParams extractCellScores
#'
#' @return A named character vector with the default slot for each supported
#' method.
#'
#' @noRd
#'
expSlots <- function(scObj){
    slot <- 'data'
    if(is(scObj, 'SingleCellExperiment'))
        slot <- 'logcounts'
    v <- setNames(c(rep(slot, 11), 'counts', rep(slot, 5)),
                  supportedMethods())
    return(v)
}
