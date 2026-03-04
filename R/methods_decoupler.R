#' @importFrom decoupleR run_mdt run_mlm run_ora run_udt
#'
NULL

#' Run a gene set analysis method using decoupleR
#'
#' This function runs one of the gene set analysis methods supported
#' by \code{decoupleR}.
#'
#' @inheritParams runGSAMethods
#' @param method Gene set analysis method.
#' @param slot Gene expression slot.
#' @param ... Additional arguments passed to gene set analysis method.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @keywords internal
#'
runDecoupleRMethod <- function(scObj,
                               geneSets,
                               method,
                               slot = 'data',
                               ...){
    allGenes <- Reduce(union, geneSets)
    checkGenes(scObj, allGenes)

    mat <- scExpMat(scObj, slot)
    scoreDF <- do.call(cbind, lapply(geneSets, function(genes){
        funStr <- paste0('run_', tolower(method))
        inputDF <- data.frame(source='geneSet',
                                      target=genes,
                                      mor=1)
        scores <- replace_na(do.call(funStr, list(mat, inputDF, ...))$score)
        scores <- safeMinmax(scores)
        return(scores)
    }))

    scObj <- attachScores(scObj, geneSets, scoreDF)
    return(scObj)
}

#' Run MDT using decoupleR
#'
#' This function runs \code{MDT} using \code{decoupleR}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional parameters passed to \code{runDecoupleRMethod}.
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' if (requireNamespace("ranger", quietly=TRUE)){
#'     scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#'     scObj <- qs2::qs_read(scoPath)
#'     gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#'     geneSets <- qs2::qs_read(gsPath)
#'     scObj <- runMDT(scObj, geneSets)
#' }
#'
#' @export
#'
runMDT <- function(scObj, geneSets, slot = 'data', ...)
    return(runDecoupleRMethod(scObj, geneSets, 'MDT', slot, ...))

#' Run MLM using decoupleR
#'
#' This function runs \code{MLM} using \code{decoupleR}.
#'
#' @inheritParams runDecoupleRMethod
#' @param ... Additional parameters passed to \code{runDecoupleRMethod}.
#'
#' @return A single-cell expression object with the results saved as a
#' metadata column.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runMLM(scObj, geneSets)
#'
#' @export
#'
runMLM <- function(scObj, geneSets, slot = 'data', ...)
    return(runDecoupleRMethod(scObj, geneSets, 'MLM', slot, ...))

#' Run ORA using decoupleR
#'
#' This function runs \code{ORA} using \code{decoupleR}.
#'
#' @inheritParams runMDT
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#' scObj <- qs2::qs_read(scoPath)
#' gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#' geneSets <- qs2::qs_read(gsPath)
#' scObj <- runORA(scObj, geneSets)
#'
#' @export
#'
runORA <- function(scObj, geneSets, slot = 'data', ...)
    return(runDecoupleRMethod(scObj, geneSets, 'ORA', slot, ...))

#' Run UDT using decoupleR
#'
#' This function runs \code{UDT} using \code{decoupleR}.
#'
#' @inheritParams runMDT
#'
#' @return A single-cell expression object with the results saved as a metadata
#' column.
#'
#' @examples
#' if (requireNamespace("rpart", quietly=TRUE)){
#'     scoPath <- system.file('extdata', 'scObj.qs2', package='GSABenchmark')
#'     scObj <- qs2::qs_read(scoPath)
#'     gsPath <- system.file('extdata', 'geneSets.qs2', package='GSABenchmark')
#'     geneSets <- qs2::qs_read(gsPath)
#'     scObj <- runUDT(scObj, geneSets)
#'     }
#'
#' @export
#'
runUDT <- function(scObj, geneSets, slot = 'data', ...)
    return(runDecoupleRMethod(scObj, geneSets, 'UDT', slot, ...))

