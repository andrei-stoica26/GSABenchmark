#' Time the evaluation of an expression
#'
#' This function computes the time required to evaluate an expression.
#'
#' @param expr An expression to be evaluated.
#' @param verbose Whether to display the elapsed time as output.
#'
#' @return A numeric value representing elapsed time.
#'
#' @keywords internal
#'
timeExpr <- function(expr, verbose=FALSE){
    elapsedTime <- as.numeric(system.time(expr)[3])
    if(verbose)
        message('Time elapsed: ', round(elapsedTime, 4), ' seconds.')
    return(elapsedTime)
}

#' Compute the elapsed time, memory and peak memory required to evaluate an
#' expression
#'
#' This function computes the elapsed time, memory and peak memory required to
#' evaluate an expression.
#'
#' @inheritParams timeExpr
#'
#' @return A numeric vector containing the elapsed time, the memory usage and
#' the peak memory usage.
#'
#' @noRd
#'
timeMemoryExpr <- function(expr, verbose=FALSE){
    startMem <- gc(verbose=FALSE, reset=TRUE)
    elapsedTime <- timeExpr(expr, verbose)
    endMem <- gc(verbose=FALSE, reset=FALSE)
    totalMbRAM <- endMem[2, 2] - startMem[2, 2]
    peakMbRAM <- endMem[2, 6] - startMem[2, 6]
    res <- setNames(c(elapsedTime, totalMbRAM, peakMbRAM),
                    c('elapsedTime', 'totalMbRAM', 'peakMbRAM'))
    return(res)
}
