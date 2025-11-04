#' @importFrom GSVA gsva gsvaParam plageParam zscoreParam
#'
NULL

#' Run a gene set analysis method using \code{GSVA}
#'
#' This function runs one of the gene set analysis methods supported
#' by \code{GSVA}.
#'
#' @inheritParams runDecoupleRMethod
#' @param invert Whether to transform the scores from x to 1 - x.
#' @param filter Whether to filter the expression matrix as to contain only
#' signature genes.
#' @param ... Additional parameters passed to other functions.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @keywords internal
#'
runGSVAMethod <- function(scObj, geneSets, method, invert = FALSE,
                          filter = FALSE, ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    if(filter)
        mat <- scExpMat(scObj, 'data', allGenes) else
            mat <- scExpMat(scObj, 'data')

    gsvaPar <- do.call(paste0(tolower(method), 'Param'),
                       list(mat, geneSets, ...))
    scoreDF <- t(gsva(gsvaPar))

    scoreDF <- apply(scoreDF, 2, safeMinmax)
    if (invert)
        scoreDF <- apply(scoreDF, 2, function(x) 1 - x)

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}

#' Run GSVA
#'
#' This function runs \code{GSVA} using \code{GSVA}.
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to \code{runGSVAMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scoPath <- system.file('testdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('testdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' scObj <- runGSVA(scObj, geneSets)
#'
#' @export
#'
runGSVA <- function(scObj, geneSets, ...)
    return(runGSVAMethod(scObj, geneSets, 'GSVA', ...))

#' Run PLAGE
#'
#' This function runs \code{PLAGE} using \code{GSVA}.
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to \code{runGSVAMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scoPath <- system.file('testdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('testdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' scObj <- runPLAGE(scObj, geneSets)
#'
#' @export
#'
runPLAGE <- function(scObj, geneSets, ...)
    return(runGSVAMethod(scObj, geneSets, 'PLAGE', TRUE, TRUE, ...))

#' Run Zscore
#'
#' This function runs \code{Zscore} using \code{GSVA}.
#'
#' @inheritParams runGSVAMethod
#' @param ... Additional parameters passed to \code{runGSVAMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scoPath <- system.file('testdata', 'scObj.qs', package='GSABenchmark')
#' scObj <- qs::qread(scoPath)
#' gsPath <- system.file('testdata', 'geneSets.qs', package='GSABenchmark')
#' geneSets <- qs::qread(gsPath)
#' scObj <- runZscore(scObj, geneSets)
#'
#' @export
#'
runZscore <- function(scObj, geneSets, ...)
    return(runGSVAMethod(scObj, geneSets, 'Zscore', FALSE, TRUE, ...))
