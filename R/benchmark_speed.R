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
#' @return A data frame listing run times
#'
#' @export
#'
extractRunTimes <- function(fileName){
  logLines <- readLines(fileName)
  messages <- logLines[seq(1, length(logLines), 2)]
  methods <- unique(sapply(messages, function(x) str_split(x, ' ')[[1]][2]))
  geneSets <- unique(sapply(messages, function(x) str_split(x, ' ')[[1]][4]))
  runTimes <- sapply(logLines[seq(2, length(logLines), 2)], function(x) str_remove(x, 'Time difference of '), USE.NAMES=FALSE)
  runTimes <- sapply(runTimes, timeToSeconds, USE.NAMES=FALSE)
  df <- data.frame(pracma::Reshape(runTimes, length(methods), length(geneSets)))
  rownames(df) <- methods
  colnames(df) <- geneSets
  df$avg <- rowMeans(df)
  return(df)
}

timeToSeconds <- function(timeStr){
  splitTime <- str_split(timeStr, ' ')[[1]]
  timeDict <- setNames(c(1, 60, 3600), c('secs', 'mins', 'hours'))
  return(as.numeric(splitTime[1]) * timeDict[splitTime[2]])
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
