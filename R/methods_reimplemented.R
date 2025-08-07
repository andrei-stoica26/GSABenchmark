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
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @keywords internal
#'
addModuleScoreHelper <- function(scObj,
                                 genes,
                                 colStr = 'AddModuleScore',
                                 slot = 'data',
                                 pool = rownames(scObj),
                                 nbin = 24,
                                 ctrl = 100){
    mat <- scExpMat(scObj, slot, densify=FALSE)
    matAvg <- Matrix::rowMeans(mat[pool, ])
    matAvg <- matAvg[order(matAvg)]
    matCut <- cut_number(matAvg + rnorm(length(matAvg))/1e30,
                         n=nbin,
                         labels=FALSE,
                         right=FALSE)
    names(matCut) <- names(matAvg)
    ctrlUse <- c()
    for (gene in genes) {
        ctrlSample <- names(sample(matCut[which(matCut == matCut[gene])],
                         size=ctrl, replace=FALSE))
        ctrlUse <- c(ctrlUse, ctrlSample)
    }
    ctrlUse <- unique(ctrlUse)
    featureScores <- Matrix::colMeans(mat[genes, , drop=FALSE])
    ctrlScores <- Matrix::colMeans(mat[ctrlUse, ])
    scores <- featureScores - ctrlScores
    scObj[[colStr]] <- safeMinmax(scores)
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
#' @export
#'
runAddModuleScore <- function(scObj, genes, colStr = 'AddModuleScore',
                              slot = 'data', pool = rownames(scObj),
                              nbin = 24, ctrl = 100, seed = 1){
    checkGenes(scObj, genes)
    if (is.null(seed))
        stop('A positive integer seed must be set.')
    return(with_seed(seed, addModuleScoreHelper(scObj,
                                                genes,
                                                colStr,
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
#' @param method One of 'oddsratio' and 'likelihood'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runJASMINE <- function(scObj,
                       genes,
                       colStr = 'JASMINE',
                       method = c('oddsratio', 'likelihood')){
    checkGenes(scObj, genes)
    method <- match.arg(method, c('oddsratio', 'likelihood'))

    mat <- scExpMat(scObj, 'data')

    ranks <- apply(mat, 2, function(x) meanGeneRank(x, genes))
    ranks <- safeMinmax(ranks)

    sigMat <- mat[genes,]
    complMat <- mat[setdiff(rownames(mat), genes), ]

    sigGeneFreq <- apply(sigMat, 2, function(x) sum(x != 0))
    complGeneFreq <- apply(complMat, 2, function(x) sum(x != 0))

    absentSigGeneFreq <- nrow(sigMat) - sigGeneFreq
    absentComplGeneFreq <- nrow(complMat) - complGeneFreq

    absentSigGeneFreq <- replace(absentSigGeneFreq, absentSigGeneFreq == 0, 1)
    complGeneFreq <- replace(complGeneFreq, complGeneFreq == 0, 1)

    if(method == 'oddsratio')
        res <- (sigGeneFreq * absentComplGeneFreq) /
        (absentSigGeneFreq * complGeneFreq) else
            res <- sigGeneFreq * (complGeneFreq + absentComplGeneFreq) /
                        (complGeneFreq * (sigGeneFreq + absentSigGeneFreq))

    res <- safeMinmax(res)

    scores <- (res + ranks) / 2
    scObj[[colStr]] <- safeMinmax(scores)

    return(scObj)
}
