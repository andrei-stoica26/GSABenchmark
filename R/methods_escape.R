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
runEscapeMethod <- function(scObj, geneSets, method, ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    mat <- scExpMat(scObj, 'data')

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
#'     scoPath <- system.file('testdata', 'scObj.qs', package='GSABenchmark')
#'     scObj <- qs::qread(scoPath)
#'     gsPath <- system.file('testdata', 'geneSets.qs', package='GSABenchmark')
#'     geneSets <- qs::qread(gsPath)
#'     scObj <- runAUCell(scObj, geneSets)
#'     }
#'
#' @export
#'
runAUCell <- function(scObj, geneSets, ...)
    return(runEscapeMethod(scObj, geneSets, 'AUCell', ...))

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
#' scoPath <- system.file('testdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('testdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' scObj <- runssGSEA(scObj, geneSets)
#'
#' @export
#'
runssGSEA <- function(scObj, geneSets, ...)
    return(runEscapeMethod(scObj, geneSets, 'ssGSEA', ...))

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
#'     scoPath <- system.file('testdata', 'scObj.qs', package='GSABenchmark')
#'     scObj <- qs::qread(scoPath)
#'     gsPath <- system.file('testdata', 'geneSets.qs', package='GSABenchmark')
#'     geneSets <- qs::qread(gsPath)
#'     scObj <- runUCell(scObj, geneSets)
#'     }
#'
#' @export
#'
runUCell <- function(scObj, geneSets, ...)
    return(runEscapeMethod(scObj, geneSets, 'UCell', ...))
