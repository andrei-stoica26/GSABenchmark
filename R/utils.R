#' Message information on gene set analysis method run
#'
#' This function messages or prints information on the gene set analysis method
#' running and the group where it is running.
#'
#' @param funStr A supported gene set analysis method
#' @param group A Seurat identity class
#' @param messageFun An output function
#'
#' @return Messaged (default) or printed output
#'
#'
messageMethod <- function(funStr, group, messageFun = message)
  return(messageFun(paste0(str_replace(funStr, 'Run', 'Running '), ' for identity: ', group, '...')))

#' End the last output diversion
#'
#' This function ends the output diversion to a file name, if one exists.
#' @param fileName File where the output is currently diverted
#'
safeSink <- function(fileName)
  if (!is.null(fileName))
    sink()
