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
                       method = c('oddsratio', 'likelihood'),
                       colStr = 'JASMINE'){
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
