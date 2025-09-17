#' Create a mock gene-cell matrix
#'
#' This function creates a mock gene-cell matrix
#'
#' @details This function is intended solely for running examples and tests
#' to illustrate the functionality of gene set analysis methods supported by
#' GSABenchmark. It is not meant to realistically represent the statistical
#' properties of scRNA-seq datasets
#'
#' @param nGenes Number of genes.
#' @param nCells Number of cells.
#' @param nAssignedElem Number of elements assigned a value selected from an
#' uniform distribution bounded by 0 and `maxElem`.
#' @param maxElem Upper bound for the assigned values.
#' @param nSelectedGenes Number of genes for which a peak of gene expression
#' will be simulated.
#' @param peakStart Cell index corresponding to the start of the
#' gene expression peak.
#' @param peakEnd Cell index corresponding to the end of the
#' gene expression peak.
#' @param peakMin The lowest gene expression value within the gene expression
#' peak.
#' @param peakMax The hightest gene expression value within the gene expression
#' peak.
#' @param seed Random seed.
#' @return A matrix with a gene expression peak.
#'
#' @examples
#' mockGCMat(seed=4)
#'
#' @export
#'
mockGCMat <- function(nGenes = 500,
                      nCells = 300,
                      nAssignedElem = 8000,
                      maxElem = 13,
                      nSelectedGenes = 200,
                      peakStart = 20,
                      peakEnd = 50,
                      peakMin = 14,
                      peakMax = 15,
                      seed = 1){

    if (maxElem < 1)
        stop('`maxElem` must be >= 1.')

    totalElem <- nGenes * nCells
    if (nAssignedElem > totalElem){
        warning('`nAssignedElem` is greater than the total number of ',
                'elements in the matrix (', totalElem,
                '). Setting it to the latter.')
        nAssignedElem <- totalElem
    }

    mat <- matrix(0, nGenes, nCells)
    rownames(mat) <- paste0('G', seq(nGenes))
    colnames(mat) <- paste0('C', seq(nCells))
    mat[seq(nAssignedElem)] <- with_seed(seed, runif(nAssignedElem,
                                                     max=maxElem))

    genes <- paste0('G', seq(nSelectedGenes))
    peakLength <- peakEnd - peakStart + 1
    mat[genes,
        seq(peakStart, peakEnd)] <- with_seed(seed,
                                              matrix(runif(
                                              nSelectedGenes * peakLength,
                                              min=peakMin,
                                              max=peakMax),
                                              nrow=nSelectedGenes,
                                              ncol=peakLength))
    return(mat)
}

#' Create mock gene sets
#'
#' This function create mock gene sets.
#'
#' @details This function is intended to be used in conjunction
#' with \code{mockGCMat}.
#'
#' @param indices Vector containing the start index and end index of each gene
#' set.
#'
#' @return A list of mock gene sets represented as character vectors.
#'
#' @examples
#' mockGeneSets(c(13, 23, 14, 49, 8, 16))
#'
#' @export
#'
mockGeneSets <- function(indices = c(1, 150, 50, 200)){
    if (length(indices) %% 2)
        stop('The length of `indices` must be even.')
    odds <- seq(1, length(indices), 2)
    evens <- seq(2, length(indices), 2)
    geneSets <- setNames(mapply(function(x, y)
        paste0('G', seq(indices[x], indices[y])),
        odds, evens, SIMPLIFY=FALSE),
        paste0('Set', seq(length(odds))))
    return(geneSets)
}
