#' @importFrom hammers shuffleGenes
#'
NULL


#' Run gene set analysis method on shuffled gene sets
#'
#' This function runs a gene set analysis method on shuffled gene sets.
#'
#' @inheritParams runGSAMethods
#' @param gsaMethod Name of the gene set analysis method.
#' @param lossVector Gene loss values.
#' @param noiseVector Noise values.
#' @param doGrid Whether to run the methods for each loss-noise combination.
#' @param nReplicates Number of replicates.
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
                             geneSetNames,
                             lossVector = c(0, 0.25, 0.5),
                             noiseVector = c(0, 0.25, 0.5),
                             doGrid = TRUE,
                             nReplicates = 3){

    checkSetNames(scObj, labelCol, geneSetNames)
    if(doGrid){
        lossNoiseDF <- expand.grid(lossVector, noiseVector)
        lossVector <- lossNoiseDF[, 1]
        noiseVector <- lossNoiseDF[, 2]
    }

    expression <- scExpMat(scObj, 'counts')
    nValues <- length(lossVector)
    for (i in seq(nValues)){
        for (j in seq(nReplicates)){
            shGeneSets <- lapply(geneSets, function(x)
                shuffleGenes(scObj, x, lossVector[i], noiseVector[i]))
            infix <- paste0('_',
                            lossVector[i] * 100,
                            '_',
                            noiseVector[i] * 100,
                            '_',
                            j)
            scObj <- runGSAMethods(scObj,
                                   labelCol,
                                   shGeneSets,
                                   gsaMethod,
                                   geneSetNames,
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
                                gsaMethod,
                                geneSetNames){
    runs <- grep(gsaMethod, colnames(metadataDF(scObj)), value=TRUE)
    runs <- unique(vapply(runs, runFromColumn, character(1)))
    return(runBenchmark(scObj,
                        labelCol,
                        geneSets=NULL,
                        gsaMethods=runs,
                        geneSetNames))
}
