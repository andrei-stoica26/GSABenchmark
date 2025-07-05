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


#' Find the last appearances of elements in vector
#'
#' This function finds the indices at which each element in the vector makes its
#' last appearance.
#'
#' @param v A vector.
#'
#' @return A vector of indices.
#'
#' @export
#'
findLastApps <- function(v)
  return(sapply(lapply(unique(v), function(x) which(v %in% x)), function(y) y[length(y)]))

#' Convert a vector to a data frame based on input row and column names
#'
#' This function converts a vector to a data frame based on input row and column names.
#' Optionally, it also calculates the row means.
#'
#' @param v A vector.
#' @param rowNames A character vector.
#' @param colNames A character vector.
#' @param addRowMeans Whether to add the row means to the data frame
#'
#' @return A data frame
#'
#' @export
#'
tabulateVector <- function(v, rowNames, colNames, addRowMeans=FALSE){
  df <- data.frame(matrix(v, length(rowNames), length(colNames)))
  rownames(df) <- rowNames
  colnames(df) <- colNames
  if(addRowMeans)
    df$avg <- rowMeans(df)
  return(df)
}
