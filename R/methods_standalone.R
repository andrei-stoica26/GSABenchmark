
#' @importFrom pagoda2 score.cells.puram
#' @importFrom singscore rankGenes simpleScore
#' @importFrom SiPSiC getPathwayScores
#' @importFrom VAM vam
#'
NULL

#' Run pagoda2
#'
#' This function runs \code{pagoda2}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to \code{pagoda2::score.cells.puram}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runPagoda2 <- function(scObj, genes, colStr = 'Pagoda2', ...){
    mat <- Matrix::t(scExpMat(scObj, 'data', densify=FALSE))
    scores <- pagoda2::score.cells.puram(mat, genes, ...)
    scObj[[colStr]] <- safeMinmax(scores)
    return(scObj)
}

#' Run Singscore
#'
#' This function runs \code{Singscore}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to \code{singscore::simpleScore}
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @export
#'
runSingscore <- function(scObj, genes, colStr = 'Singscore', ...){
    mat <- scExpMat(scObj, 'data')
    mat <- singscore::rankGenes(mat)
    scores <- singscore::simpleScore(mat, genes, ...)$TotalScore
    scObj[[colStr]] <- safeMinmax(scores)
    return(scObj)
}

#' Run SiPSiC
#'
#' This function runs \code{SiPSiC}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to \code{SiPSiC::getPathwayScores}
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column
#'
#' @export
#'
runSiPSiC <- function(scObj, genes, colStr = 'SiPSiC', ...){
    mat <- scExpMat(scObj, 'counts', densify=FALSE)
    scores <- SiPSiC::getPathwayScores(mat, genes, ...)[[2]]
    scObj[[colStr]] <- safeMinmax(scores)
    return(scObj)
}

#' Run VAM
#'
#' This function runs \code{VAM}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to \code{VAM::vam}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column
#'
#' @export
#'
runVAM <- function(scObj, genes, colStr = 'VAM', ...){
    mat <- t(scExpMat(scObj, 'data', genes))
    v <- vam(mat, ...)
    scObj[[colStr]] <- v$cdf.value
    return(scObj)
}
