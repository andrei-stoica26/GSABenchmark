#' @importFrom stringr str_remove str_split str_sub str_replace
#' @importFrom pracma Reshape
#' @importFrom nnspat euc.dist
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

#' Print a function's running time
#'
#' This function prints a function's running time, also offering an option to
#' message the time. The latter allows users to redirect the output to a file
#' while still visualizing it as a message in the console.
#'
#' @param fun Function
#' @param doMessage Whether the time should be messaged in addition to being
#' printed
#' @param ... Function parameters
#'
#' @return What the input function returns
#'
#' @export
#'
timeCode <- function(fun, doMessage = FALSE, ...){
  x <- Sys.time()
  res <- fun(...)
  y <- Sys.time()
  print(y - x)
  if(doMessage)
    message(y - x)
  return(res)
}
