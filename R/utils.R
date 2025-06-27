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


#' Take the square root of a non-negative number or the additive inverse of the
#' square root of the additive inverse of a negative number
#'
#' This function takes the square root of a non-negative number or the additive inverse of the
#' square root of the additive inverse of a negative number
#'
#' @param x A numeric variable
#'
#' @return The square root of x if positive and of -x if negative
#'
#' @export
#'
sqrt2 <- function(x)
  if (x >= 0) return(sqrt(x)) else return(-sqrt(-x))

