#' Removes genes from vector
#'
#' This function removes a fraction of genes from vector
#'
#' @inheritParams runDecoupleRMethod
#' @param fraction Fraction of genes to be removed
#'
#' @return Genes vector after the removals
#'
#' @export
#'
removeGenes <- function(genes, lossFraction = 0.5){
  nRemovedGenes <- round(lossFraction * length(genes))
  if (nRemovedGenes < 1)
    stop(paste0('No genes can be removed at an input loss fraction of ', lossFraction, '. Choose a higher value'))
  if (nRemovedGenes > length(genes))
    stop(paste0('No genes can be retained at input loss fraction of ', lossFraction, '. Choose a lower value'))
  newGenes <- c(sample(genes, length(genes) - nRemovedGenes))
  return(newGenes)
}

#' Replaces genes from vector
#'
#' This function removes a fraction of genes from vector and adds other genes from
#' the Seurat object to the vector (not necessarily as many as the removed genes).
#'
#' @inheritParams runDecoupleRMethod
#' @param lossfraction Fraction of genes to be replaced
#' @param finalSizeFactor The size of the final vector relative to the original
#' one. Values greater than 1 indicate that more genes will be added than removed,
#' while values greater lower than 1 indicate otherwise
#' @param geneCountThresh Minimum number of cells in which newly added genes must
#' be expressed
#'
#' @return Genes vector after the replacements
#'
#' @export
#'
noisifyGenes <- function(seuratObj, genes, lossFraction = 0.5, finalSizeFactor = 1, geneCountThresh = 10){
  expression <- LayerData(seuratObj, layer='counts')
  freq <- rowSums(expression != 0)
  suitableGenes <- names(freq[freq >= geneCountThresh])
  genesComplement <- setdiff(suitableGenes, genes)
  nRemovedGenes <- round(lossFraction * length(genes))
  if(finalSizeFactor < 1 - lossFraction)
    stop('finalSizeFactor must be greater than 1 - lossFraction')
  if (nRemovedGenes > length(genes))
    stop(paste0('No genes can be retained at input loss fraction of ', lossFraction, '. Choose a lower value'))
  nRetainedGenes <- length(genes) - nRemovedGenes
  nAddedGenes <- round(finalSizeFactor * length(genes)) - nRetainedGenes
  newGenes <- c(sample(genes, nRetainedGenes), sample(genesComplement, nAddedGenes))
  return(newGenes)
}
