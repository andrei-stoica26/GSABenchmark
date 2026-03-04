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
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
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

    nGeneSets <- length(geneSets)
    nMethods <- length(gsaMethods)
    arraySize <- nGeneSets * nMethods

    elapsedTime <- rep(-1, arraySize)
    peakMemory <- rep(-1, arraySize)

    for (i in seq_along(geneSets))
        for (j in seq_along(gsaMethods)){
            method <- gsaMethods[j]
            safeMessage(paste0(
                'Computing running time and peak memory usage for ',
                method, ' on ', geneSetNames[i], ' genes...'), verbose)
            res <- timeMemoryExpr(silently_run(do.call(
                paste0('run', method),
                list(scObj, geneSets[i]))),
                verbose)

            arrayPos <- (i - 1) * nMethods + j
            elapsedTime[arrayPos] <- res[1]
            peakMemory[arrayPos] <- res[3]
            }
    time <- tabulateVector(elapsedTime, gsaMethods, geneSetNames, TRUE, TRUE)
    space <- tabulateVector(peakMemory, gsaMethods, geneSetNames, TRUE, TRUE)
    return(list(time = time, space = space))
}

#' Convert a vector to a data frame based on input row and column names
#'
#' This function converts a vector to a data frame based on input row and
#' column names. Optionally, it also calculates the row means.
#'
#' @param v A vector.
#' @param rowNames A character vector.
#' @param colNames A character vector.
#' @param addRowMeans Whether to add the row means to the data frame.
#' @param sortByRowMeans Whether to sort by row means.
#'
#' @return A data frame.
#'
#' @examples
#' v <- c(2, 3, 4, 19, 15, 25, 32, 8)
#' res <- tabulateVector(v, paste0('r', seq(4)), paste0('c', seq(2)))
#'
#' @export
#'
tabulateVector <- function(v,
                           rowNames,
                           colNames,
                           addRowMeans=FALSE,
                           sortByRowMeans=FALSE){
    df <- data.frame(matrix(v, length(rowNames), length(colNames)))
    rownames(df) <- rowNames
    colnames(df) <- colNames
    if(addRowMeans){
        df$avg <- rowMeans(df)
        if(sortByRowMeans)
            df <- df[order(df$avg), ]
    }
    return(df)
}
