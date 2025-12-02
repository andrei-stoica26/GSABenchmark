#' @importFrom stringr str_remove str_replace str_split str_sub str_to_title
#'
NULL

#' Extract running times and peak memory usage for gene set analysis methods
#'
#' This function extracts running times and peak memory usage for gene set
#' analysis methods.
#'
#' @inheritParams allBenchmarkResults
#' @param verbose Whether output should be verbose.
#'
#' @return A list of two data frames, the first comprising running times in
#' seconds, the second comprising peak memory usage in mebibytes.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' ef <- efficiencyBenchmark(scObj, 'label', geneSets, c('CSOA', 'Zscore'))
#'
#' @export
#'
efficiencyBenchmark <- function(scObj,
                                labelCol,
                                geneSets,
                                gsaMethods,
                                checkLabels = TRUE,
                                verbose = TRUE){

    geneSetNames <- names(geneSets)

    if (checkLabels)
        checkSetNames(scObj, labelCol, geneSetNames)

    elapsedTime <- c()
    peakMemory <- c()
    for (i in seq_along(geneSets))
        for (method in gsaMethods){
            if (verbose)
                message('Computing running time and peak memory usage for ',
                           method, ' on ', geneSetNames[i], ' genes...')
            res <- timeMemoryExpr(silently_run(do.call(
                paste0('run', method),
                list(scObj, geneSets[i]))),
                verbose)
            elapsedTime <- c(elapsedTime, res[1])
            peakMemory <- c(peakMemory, res[3])
            }
    time <- tabulateVector(elapsedTime, gsaMethods, geneSetNames, TRUE, TRUE)
    space <- tabulateVector(peakMemory, gsaMethods, geneSetNames, TRUE, TRUE)
    return(list(time = time, space = space))
}
