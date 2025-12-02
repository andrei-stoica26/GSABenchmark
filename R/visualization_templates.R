#' Create a title for a plot selected from a list if a suffix is provided
#'
#' This function creates a title to be assigned to a plot selected from a list
#' if a suffix is provided. Otherwise, it returns \code{NULL}.
#'
#' @param prefix Title prefix.
#' @param infix Title infix.
#' @param suffix Title suffix.
#'
#' @return A subplot title or \code{NULL}..
#'
#' @noRd
#'
subplotTitle <- function(prefix, infix, suffix=NULL){
    if(is.null(suffix))
        return(NULL)
    if(is.null(prefix))
        return(paste0(infix, ' - ', suffix))
    return(paste0(prefix, ' (', infix, ') - ', suffix))
}

#' Create a list of plots from a list of arguments
#'
#' This function creates a list of plots from a list of arguments by first
#' converting the arguments into results that can be plotted
#' using \code{calcFun}, then plotting these results using \code{plotFun}.
#'
#' @param args A list of arguments passed to \code{calcFun}.
#' @param calcFun The function used for generating the results that will be
#' plotted
#' @param plotFun The function used for plotting.
#' @param unlistPlotArgs Whether the additional arguments passed to
#' \code{plotFun} have to be unlisted (that is, they are provided as a list).
#' @param titlePrefix Prefix used to create a title with \code{subplotTitle}.
#' Ignored if \code{titleSuffix} is \code{NULL}.
#' @param titleInfixes A character vector of title infixes used to create a
#' title with \code{subplotTitle}. Ignored if \code{titleSuffix}
#' is \code{NULL}.
#' @param titleSuffix Suffix used to create plot titles If \code{NULL}, the
#' plots will have no titles.
#' @param ... Additional arguments passed to \code{plotFun}.
#'
#' @return A named list of ggplot objects.
#'
#' @keywords internal
#'
plotList <- function(args,
                     calcFun,
                     plotFun,
                     unlistPlotArgs=FALSE,
                     titlePrefix=NULL,
                     titleInfixes=NULL,
                     titleSuffix=NULL,
                     ...){

    plotInput <- do.call(calcFun, args)
    plotNames <- names(plotInput)
    plotArgs <- list(...)
    if(unlistPlotArgs)
        plotArgs <- unlist(plotArgs, recursive=FALSE)

    if(is.null(titleInfixes))
        titleInfixes <- plotNames

    plots <- setNames(mapply(function(df, titleInfix) {
        title <- subplotTitle(titlePrefix, titleInfix, titleSuffix)
        plotArgs <- c(list(df, title), plotArgs)
        return(do.call(plotFun, plotArgs))
    }, plotInput, titleInfixes, SIMPLIFY=FALSE), plotNames)
    return(plots)
}
