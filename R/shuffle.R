#' @importFrom hammers shuffleGenes
#'
NULL


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
#' to \code{hammers::shuffleGenes}. Its length decides the number of
#' replicates.
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object
#' with the results of the runs stored as metadata columns.
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
                             seeds = c(1, 2, 3)){

    checkSetNames(scObj, labelCol, names(geneSets))
    if(doGrid){
        lossNoiseDF <- expand.grid(loss, noise)
        lossVector <- lossNoiseDF[, 1]
        noiseVector <- lossNoiseDF[, 2]
    }

    nReplicates <- length(seeds)

    expression <- scExpMat(scObj, 'counts', densify=FALSE)
    nValues <- length(lossVector)
    for (i in seq(nValues)){
        for (j in seq(nReplicates)){
            lossPerc <- round(lossVector[i] * 100, 1)
            noisePerc <- round(noiseVector[i] * 100, 1)
            message('Shuffling genes: gene loss = ', lossPerc,
                    '%, noise = ', noisePerc, '%, replicate = ', j, '.')
            shGeneSets <- lapply(geneSets, function(x)
                shuffleGenes(scObj, x, lossVector[i], noiseVector[i],
                             seed=seeds[j], verbose=FALSE))
            infix <- paste0('_', lossPerc, '_', noisePerc, '_', j)
            scObj <- runGSAMethods(scObj,
                                   labelCol,
                                   shGeneSets,
                                   gsaMethod,
                                   infix)
        }

    }
    return(scObj)
}

#' Generate all benchmark results for shuffled gene sets
#
#' This function generates all benchmark results for shuffled gene sets.
#'
#' @inheritParams runGSAMethods
#' @inheritParams runMethodShuffle
#'
#' @return A list of benchmark results.
#'
#' @export
#'
runBenchmarkShuffle <- function(scObj,
                                labelCol,
                                geneSets,
                                gsaMethod){
    runs <- grep(gsaMethod, colnames(metadataDF(scObj)), value=TRUE)
    runs <- unique(vapply(runs, runFromColumn, character(1)))
    return(runBenchmark(scObj,
                        labelCol,
                        geneSets,
                        gsaMethods=runs))
}
