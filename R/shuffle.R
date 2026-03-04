#' Run gene set analysis method on shuffled gene sets
#'
#' This function runs a gene set analysis method on shuffled gene sets.
#'
#' @inheritParams runGSAMethods
#' @param gsaMethod Name of the gene set analysis method.
#' @param loss A numeric vector of gene loss values. Must be in [0, 1).
#' @param noise A numeric vector of noise values. Must be in [0, 1).
#' @param doGrid Whether to run the methods for each loss-noise combination.
#' @param averageReplicates Whether to average replicates. If \code{TRUE},
#' a single set of scores will be returned for each loss-noise combination
#' for which the method is run. If \code{FALSE}, distinct sets of scores will
#' be returned for each replicate.
#' @param seeds A numeric vector of random seeds passed
#' to \code{hammers::shuffleGenes}. Its length determines the number of
#' replicates.
#'
#' @return A \code{Seurat} or \code{SingleCellExpression} object
#' with the results of the runs stored as metadata columns.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
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
                             averageReplicates = TRUE,
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
        lossPerc <- round(loss[i] * 100, 1)
        noisePerc <- round(noise[i] * 100, 1)
        runTemplate <- paste0('_', lossPerc, '_', noisePerc)
        replTemplate <- paste0(runTemplate, '_')
        for (j in seq(nReplicates)){
            message('Shuffling genes: gene loss = ', lossPerc,
                    '%, noise = ', noisePerc, '%, replicate = ', j, '.')
            shGeneSets <- lapply(geneSets, function(x)
                shuffleGenes(scObj, x, loss[i], noise[i],
                             seed=seeds[j], verbose=FALSE))
            infix <- paste0(replTemplate, j)
            scObj <- runGSAMethods(scObj,
                                   labelCol,
                                   shGeneSets,
                                   gsaMethod,
                                   infix,
                                   outputFun)
        }

        if (averageReplicates){
            message('Averaging replicates...')
            for (gsName in names(geneSets)){
                replCols <- paste0(gsaMethod, replTemplate, seq(nReplicates),
                                   '_', gsName)
                runScore <- rowMeans(metadataDF(scObj)[, replCols])
                for (col in replCols)
                    scCol(scObj, col) <- NULL
                scCol(scObj, paste0(gsaMethod, runTemplate,
                                    '_', gsName)) <- runScore

            }


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
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
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

#' Replaces genes from vector
#'
#' This function removes and adds genes from vector at random.
#'
#' @inheritParams extractCellScores
#' @param genes A character vector of genes.
#' @param lossFrac Fraction of genes than be removed. Must be in \code{[0, 1]}.
#' @param noiseFrac Amount of noise (random genes) in the final gene vector.
#' Must be in \code{[0, 1)}
#' @param geneCountThresh Minimum number of cells in which newly added genes
#' must be expressed.
#' @param seed Random seed.
#' @param verbose Whether the output should be verbose.
#'
#' @return Genes vector after changes.
#'
#' @noRd
#'
shuffleGenes <- function(scObj, genes, lossFrac, noiseFrac,
                         geneCountThresh = 10, seed = 1,
                         verbose = TRUE){
    nGenes <- length(genes)
    nRemovedGenes <- round(lossFrac * nGenes)
    nRetainedGenes <- nGenes - nRemovedGenes

    expression <- scExpMat(scObj, 'counts')
    freq <- rowSums(expression != 0)

    suitableGenes <- setdiff(names(freq[freq >= geneCountThresh]), genes)

    if(lossFrac > 0){
        genes <- c(with_seed(seed, sample(genes, nRetainedGenes)))
        safeMessage(paste0('Removed ', nRemovedGenes, ' genes.'), verbose)
    } else
        safeMessage('No genes were removed.', verbose)

    if(noiseFrac > 0){
        nAddedGenes <- round(noiseFrac * nRetainedGenes / (1 - noiseFrac))
        genes <- c(genes, with_seed(seed, sample(suitableGenes, nAddedGenes)))
        safeMessage(paste0('Added ', nAddedGenes, ' random genes.'), verbose)
    } else
        safeMessage('No genes were added.', verbose)

    return(genes)
}
