#' @importFrom stringr str_remove str_replace str_split str_sub str_to_title
#' @importFrom peakRAM peakRAM
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
#' scoPath <- system.file('testdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('testdata', 'geneSets.qs', package='GSABenchmark')
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
            methodCall <- paste0('run', method)
            if (verbose)
                message('Computing running time and peak memory usage for ',
                           method, ' on ', geneSetNames[i], ' genes...')
            geneSet <- geneSets[[i]]
            names(geneSet) <- geneSetNames[i]
            df <- peakRAM(silently_run(do.call(methodCall,
                                               list(scObj, geneSet))))
            if(verbose)
                message('Time elapsed: ', round(df$Elapsed_Time_sec, 4),
                        ' seconds')
            elapsedTime <- c(elapsedTime, df$Elapsed_Time_sec)
            peakMemory <- c(peakMemory, df$Peak_RAM_Used_MiB)
            }
    time <- tabulateVector(elapsedTime, gsaMethods, geneSetNames, TRUE, TRUE)
    space <- tabulateVector(peakMemory, gsaMethods, geneSetNames, TRUE, TRUE)
    return(list(time = time, space = space))
}
