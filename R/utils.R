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

#' Add two objects as if using the + operator
#'
#' This function adds two objects as if using the + operator. Useful to use in
#' functions like Reduce
#'
#' @param x First object
#' @param y Second object
#'
#' @return Sum of the two objects
#'
plus <- function(x, y)
  return(x + y)


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

#' Take the square root of a non-negative number or the additive inverse of the
#' square root of the additive inverse of a negative number
#'
#' This function takes the square root of a non-negative number or the additive inverse of the
#' square root of the additive inverse of a negative number
#'
#' @param x A numeric variable
#'
#' @export
#'
sqrt2 <- function(x)
  if (x >= 0) return(sqrt(x)) else return(-sqrt(-x))

