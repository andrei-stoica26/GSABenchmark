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
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runPagoda2(scObj, geneSets)
#'
#' @export
#'
runPagoda2 <- function(scObj, geneSets, slot = 'data', ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)
    mat <- Matrix::t(scExpMat(scObj, slot, densify=FALSE))

    scoreDF <- do.call(cbind, lapply(geneSets, function(genes){
        scores <- pagoda2::score.cells.puram(mat, genes, ...)
        scores <- safeMinmax(scores)
        return(scores)
    }))

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}

#' Run Singscore
#'
#' This function runs \code{Singscore}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to \code{singscore::simpleScore}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runSingscore(scObj, geneSets)
#'
#' @export
#'
runSingscore <- function(scObj, geneSets, slot = 'data', ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    mat <- scExpMat(scObj, slot)
    mat <- singscore::rankGenes(mat)

    scoreDF <- do.call(cbind, lapply(geneSets, function(genes){
        scores <- singscore::simpleScore(mat, genes, ...)$TotalScore
        scores <- safeMinmax(scores)
        return(scores)
    }))

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}

#' Run SiPSiC
#'
#' This function runs \code{SiPSiC}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional arguments passed to \code{SiPSiC::getPathwayScores}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runSiPSiC(scObj, geneSets)
#'
#' @export
#'
runSiPSiC <- function(scObj, geneSets, slot = 'counts', ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    mat <- scExpMat(scObj, slot, densify=FALSE)

    scoreDF <- do.call(cbind, lapply(geneSets, function(genes){
        scores <- getPathwayScores(mat, genes, ...)[[2]]
        scores <- safeMinmax(scores)
        return(scores)
    }))

    scObj <- attachScores(scObj, geneSets, scoreDF)
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
#' column.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runVAM(scObj, geneSets)
#'
#' @export
#'
runVAM <- function(scObj, geneSets, slot = 'data', ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    mat <- Matrix::t(scExpMat(scObj, slot, allGenes))
    scoreDF <- do.call(cbind, lapply(geneSets, function(genes){
        scores <- vam(mat[, genes], ...)$cdf.value
        return(scores)
    }))

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}
