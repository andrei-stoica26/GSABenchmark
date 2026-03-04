#' @importFrom escape escape.matrix
#'
NULL

#' Run a gene set analysis method using escape
#'
#' This function runs one of the gene set analysis methods supported
#' by \code{escape}.
#'
#' @inheritParams runDecoupleRMethod
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @keywords internal
#'
runEscapeMethod <- function(scObj, geneSets, method, slot = 'data', ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    mat <- scExpMat(scObj, slot)

    scoreDF <- escape.matrix(mat, geneSets, method=method, ...)
    scoreDF <- apply(scoreDF, 2, safeMinmax)

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}

#' Run AUCell
#'
#' This function runs \code{AUCell} using \code{escape}.
#'
#' @inheritParams runEscapeMethod
#' @param ... Additional parameters passed to \code{runEscapeMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' if (requireNamespace("AUCell", quietly=TRUE)){
#'     scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#'     scObj <- qs2::qs_read(scoPath)
#'     gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#'     geneSets <- qs2::qs_read(gsPath)
#'     scObj <- runAUCell(scObj, geneSets)
#'     }
#'
#' @export
#'
runAUCell <- function(scObj, geneSets, slot = 'data', ...)
    return(runEscapeMethod(scObj, geneSets, 'AUCell', slot, ...))

#' Run ssGSEA
#'
#' This function runs \code{ssGSEA} using \code{escape}.
#'
#' @inheritParams runAUCell
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runssGSEA(scObj, geneSets)
#'
#' @export
#'
runssGSEA <- function(scObj, geneSets, slot = 'data', ...)
    return(runEscapeMethod(scObj, geneSets, 'ssGSEA', slot, ...))

#' Run UCell
#'
#' This function runs \code{UCell} using \code{escape}.
#'
#' @inheritParams runAUCell
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' if (requireNamespace("UCell", quietly=TRUE)){
#'     scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#'     scObj <- qs2::qs_read(scoPath)
#'     gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#'     geneSets <- qs2::qs_read(gsPath)
#'     scObj <- runUCell(scObj, geneSets)
#'     }
#'
#' @export
#'
runUCell <- function(scObj, geneSets, slot = 'data', ...)
    return(runEscapeMethod(scObj, geneSets, 'UCell', slot, ...))
