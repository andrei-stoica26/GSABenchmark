#' Run gene set analysis method on shuffled gene sets
#'
#' This function runs a gene set analysis method on shuffled gene sets.
#'
#' @inheritParams runGSAMethods
#' @param gsaMethod Name of the gene set analysis method.
#' @param loss A numeric vector of gene loss values. Must be in [0, 1).
#' @param noise A numeric vector of noise values. Must be in [0, 1).
#' @param doGrid Whether to run the methods for each loss-noise combination.
#' @param seeds A numeric vector of random seeds passed
#' to \code{hammers::shuffleGenes}. Its length determines the number of
#' replicates.
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object
#' with the results of the runs stored as metadata columns.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' scObj <- runMethodShuffle(scObj, 'label', geneSets, 'CSOA', 0.2, 0.2)
#'
#' @export
#'
runMethodShuffle <- function(scObj,
                             labelCol,
                             geneSets,
                             gsaMethod,
                             loss = c(0, 0.2),
                             noise = c(0, 0.2),
                             doGrid = TRUE,
                             seeds = c(1, 2, 3),
                             outputFun = identity){

    checkSetNames(scObj, labelCol, names(geneSets))
    if(doGrid){
        lossNoiseDF <- expand.grid(loss, noise)
        loss <- lossNoiseDF[, 1]
        noise <- lossNoiseDF[, 2]
    }

    nReplicates <- length(seeds)
    nValues <- length(loss)
    for (i in seq(nValues)){
        for (j in seq(nReplicates)){
            lossPerc <- round(loss[i] * 100, 1)
            noisePerc <- round(noise[i] * 100, 1)
            message('Shuffling genes: gene loss = ', lossPerc,
                    '%, noise = ', noisePerc, '%, replicate = ', j, '.')
            shGeneSets <- lapply(geneSets, function(x)
                shuffleGenes(scObj, x, loss[i], noise[i],
                             seed=seeds[j], verbose=FALSE))
            infix <- paste0('_', lossPerc, '_', noisePerc, '_', j)
            scObj <- runGSAMethods(scObj,
                                   labelCol,
                                   shGeneSets,
                                   gsaMethod,
                                   infix,
                                   outputFun)
        }

    }
    return(scObj)
}

#' Generate all benchmark results for shuffled gene sets
#
#' This function generates all benchmark results for shuffled gene sets.
#'
#' @inheritParams runBenchmark
#' @inheritParams runMethodShuffle
#'
#' @return A list of benchmark results.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' scObj <- runMethodShuffle(scObj, 'label', geneSets, 'CSOA', 0.2, 0.2)
#' smr <- runBenchmarkShuffle(scObj, 'label', geneSets, 'CSOA', FALSE)
#'
#' @export
#'
runBenchmarkShuffle <- function(scObj,
                                labelCol,
                                geneSets,
                                gsaMethod,
                                runEFBenchmark = TRUE){
    runs <- grep(gsaMethod, colnames(metadataDF(scObj)), value=TRUE)
    runs <- unique(vapply(runs, runFromColumn, character(1)))
    return(runBenchmark(scObj,
                        labelCol,
                        geneSets,
                        gsaMethods=runs,
                        runEFBenchmark))
}
