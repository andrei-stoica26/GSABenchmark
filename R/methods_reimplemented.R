#' @importFrom Matrix colMeans rowMeans rowSums t
#' @importFrom withr with_seed
#'
NULL

#' Helper used to run AddModuleScore
#'
#' This function is used to run \code{AddModuleScore}.
#'
#' @inheritParams runDecoupleRMethod
#' @param slot Gene expression slot. Default is 'data'.
#' @param pool The set from which features to be compared with signature genes
#' are selected. Defaults to all features.
#' @param nbin Number of bins of aggregate expression levels for pool features.
#' @param ctrl Number of control features chosen from the same bin for each
#' feature.
#'
#' @return A single-cell expression object with the results saved as a
#' metadata column.
#'
#' @keywords internal
#'
addModuleScoreHelper <- function(scObj,
                                 geneSets,
                                 slot = 'data',
                                 pool = rownames(scObj),
                                 nbin = 24,
                                 ctrl = 100){
    mat <- scExpMat(scObj, slot, densify=FALSE)
    matAvg <- Matrix::rowMeans(mat[pool, ])
    matAvg <- matAvg[order(matAvg)]
    matCut <- cut_number(matAvg + rnorm(length(matAvg)) / 1e30,
                         n=nbin,
                         labels=FALSE,
                         right=FALSE)
    names(matCut) <- names(matAvg)

    scoreDF <- do.call(cbind, lapply(geneSets, function(genes){
        ctrlUse <- c()
        for (gene in genes){
            ctrlPop <- matCut[which(matCut == matCut[gene])]
            ctrlSample <- names(sample(ctrlPop,
                                       size=min(ctrl, length(ctrlPop)),
                                       replace=FALSE))
            ctrlUse <- c(ctrlUse, ctrlSample)
        }
        ctrlUse <- unique(ctrlUse)
        featureScores <- Matrix::colMeans(mat[genes, , drop=FALSE])
        ctrlScores <- Matrix::colMeans(mat[ctrlUse, ])
        scores <- featureScores - ctrlScores
        scores <- safeMinmax(scores)
        return(scores)
    }))

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}

#' Run AddModuleScore
#'
#' This function runs \code{AddModuleScore}.
#'
#' @details Reimplemented from \code{https://github.com/satijalab}.
#'
#' @inheritParams addModuleScoreHelper
#' @param seed Random seed.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scObj <- qs::qread('inst/testdata/scObj.qs')
#' geneSets <- qs::qread('inst/testdata/geneSets.qs')
#' scObj <- runAddModuleScore(scObj, geneSets)
#'
#' @export
#'
runAddModuleScore <- function(scObj,
                              geneSets,
                              slot = 'data',
                              pool = rownames(scObj),
                              nbin = 24,
                              ctrl = 100,
                              seed = 1){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    if (is.null(seed))
        stop('A positive integer seed must be set.')
    return(with_seed(seed, addModuleScoreHelper(scObj,
                                                geneSets,
                                                slot,
                                                pool,
                                                nbin,
                                                ctrl)))
}

#' Compute mean signature gene rank
#'
#' This function computes the average gene rank for input genes in a cell,
#' normalized by the number of distinct genes detected for the cell.
#'
#' @param cellVector A numeric vector representing the expression data for a
#' cell.
#' @param genes A character vector.
#' @return A numeric value.
#'
#' @noRd
#'
meanGeneRank <- function(cellVector, genes){
    posCellVector <- cellVector[cellVector != 0]
    geneRanks <- rank(posCellVector)
    sigGeneRanks <- geneRanks[which(names(geneRanks) %in% genes)]
    if (!length(sigGeneRanks))
        return(0)
    return(mean(sigGeneRanks) / length(posCellVector))
}

#' Run JASMINE
#'
#' This function runs \code{JASMINE}.
#'
#' @details Reimplemented from \code{https://github.com/NNoureen/JASMINE}.
#' The method's paper can be found at
#' \code{https://doi.org/10.7554/eLife.71994}.
#'
#' @inheritParams runDecoupleRMethod
#' @param method One of 'oddsratio' and 'likelihood'.
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scObj <- qs::qread('inst/testdata/scObj.qs')
#' geneSets <- qs::qread('inst/testdata/geneSets.qs')
#' scObj <- runJASMINE(scObj, geneSets)
#'
#' @export
#'
runJASMINE <- function(scObj,
                       geneSets,
                       method = c('oddsratio', 'likelihood')){

    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    method <- match.arg(method, c('oddsratio', 'likelihood'))

    mat <- scExpMat(scObj, 'data')
    scoreDF <- do.call(cbind, lapply(geneSets, function(genes){
        ranks <- apply(mat, 2, function(x) meanGeneRank(x, genes))
        ranks <- safeMinmax(ranks)

        sigMat <- mat[genes,]
        complMat <- mat[setdiff(rownames(mat), genes), ]

        sigGeneFreq <- apply(sigMat, 2, function(x) sum(x != 0))
        complGeneFreq <- apply(complMat, 2, function(x) sum(x != 0))

        absentSigGeneFreq <- nrow(sigMat) - sigGeneFreq
        absentComplGeneFreq <- nrow(complMat) - complGeneFreq

        absentSigGeneFreq <- replace(absentSigGeneFreq,
                                     absentSigGeneFreq == 0, 1)
        complGeneFreq <- replace(complGeneFreq, complGeneFreq == 0, 1)

        if(method == 'oddsratio')
            res <- (sigGeneFreq * absentComplGeneFreq) /
            (absentSigGeneFreq * complGeneFreq) else
                res <- sigGeneFreq * (complGeneFreq + absentComplGeneFreq) /
            (complGeneFreq * (sigGeneFreq + absentSigGeneFreq))

        res <- safeMinmax(res)

        scores <- (res + ranks) / 2
        scores <- safeMinmax(scores)
        return(scores)
    }))

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}
