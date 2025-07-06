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
removeGenes <- function(genes, fraction = 0.5){
  nRemovedGenes <- round(fraction * length(genes))
  if (nRemovedGenes < 1)
    stop(paste0('No genes can be removed at an input fraction of ', fraction, '. Choose a higher value'))
  if (nRemovedGenes > length(genes))
    stop(paste0('No genes can be retained at input fraction of ', fraction, '. Choose a lower value'))
  newGenes <- c(sample(genes, length(genes) - nRemovedGenes))
  return(newGenes)
}

#' Replaces genes from vector
#'
#' This function replaces a fraction of genes from vector with other genes from
#' the Seurat object
#'
#' @inheritParams runDecoupleRMethod
#' @param fraction Fraction of genes to be replaced
#'
#' @return Genes vector after the replacements
#'
#' @export
#'
replaceGenes <- function(seuratObj, genes, fraction = 0.5){
  genesComplement <- setdiff(rownames(seuratObj), genes)
  nReplacedGenes <- round(fraction * length(genes))
  if (nReplacedGenes < 1)
    stop(paste0('No genes can be replaced at an input fraction of ', fraction, '. Choose a higher value'))
  if (nReplacedGenes > length(genes))
    stop(paste0('No genes can be retained at input fraction of ', fraction, '. Choose a lower value'))
  newGenes <- c(sample(genes, length(genes) - nReplacedGenes),
                sample(genesComplement, nReplacedGenes))
  return(newGenes)
}


