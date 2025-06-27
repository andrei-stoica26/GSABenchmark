#' @importFrom stringr str_remove str_split str_sub str_replace
#' @importFrom pracma Reshape
#'
NULL

#' Extract running times for gene set analysis methods
#'
#' This function extracts running times for gene set analysis methods from a
#' text file generated using clusterRun
#'
#' @param fileName A text file generated with clusterRun
#'
#' @return An expression matrix of matrix class
#'
#' @export
#'
clusterRunTimes <- function(fileName){
  logLines <- readLines(fileName)
  logLines <- logLines[grep('Using', logLines, invert = T)]
  messages <- logLines[seq(1, length(logLines), 2)]
  methods <- unique(sapply(messages, function(x) str_split(x, ' ')[[1]][3], USE.NAMES = F))
  clusters <- paste0('Cluster', unique(sapply(messages, function(x) str_sub(str_split(x, ' ')[[1]][6], end = -5), USE.NAMES = F)))
  runTimes <- sapply(logLines[seq(2, length(logLines), 2)], function(x) str_remove(x, 'Time difference of '), USE.NAMES = F)
  df <- data.frame(t(pracma::Reshape(runTimes, length(methods), length(clusters))))
  rownames(df) <- clusters
  colnames(df) <- methods
  return(df)
}

#' Time code to run a function and optionally save the results in a log file
#'
#' This function times the code to run a function and either prints the output
#' into the console or saves in a logFile
#'
#' @param fun The function that will be tamed
#' @param logPrint Whether to print to a log file
#' @param ... Parameters passed to the time function
#'
#' @export
#'
timeCode <- function(fun, logPrint = FALSE, ...){
  x <- Sys.time()
  res <- fun(...)
  y <- Sys.time()
  if(!logPrint)
    print(y - x) else log_print(y - x, console=F)
  return(res)
}
