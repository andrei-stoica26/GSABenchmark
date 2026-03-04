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
runGSVAMethod <- function(scObj,
                          geneSets,
                          method,
                          slot = 'data',
                          invert = FALSE,
                          filter = FALSE,
                          ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    if(filter)
        mat <- scExpMat(scObj, slot, allGenes) else
            mat <- scExpMat(scObj, slot)

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
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runGSVA(scObj, geneSets)
#'
#' @export
#'
runGSVA <- function(scObj, geneSets, slot = 'data', ...)
    return(runGSVAMethod(scObj, geneSets, 'GSVA', slot, ...))

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
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runPLAGE(scObj, geneSets)
#'
#' @export
#'
runPLAGE <- function(scObj, geneSets, slot = 'data', ...)
    return(runGSVAMethod(scObj, geneSets, 'PLAGE', slot, TRUE, TRUE, ...))

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
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runZscore(scObj, geneSets)
#'
#' @export
#'
runZscore <- function(scObj, geneSets, slot = 'data', ...)
    return(runGSVAMethod(scObj, geneSets, 'Zscore', slot, FALSE, TRUE, ...))
