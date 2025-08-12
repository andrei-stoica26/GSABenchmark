#' Check if any gene set names are not found in the identity class column
#'
#' This function check if any gene set names are not found in the identity
#' class column.
#'
#' @inheritParams runBenchmark
#'
#' @return None. This function is called for its side effect.
#'
#' @noRd
#'
checkSetNames <- function(scObj, labelCol, geneSetNames){
    if (is(labelCol)[1] != 'character')
        stop('geneSetNames must be a character vector.')
    if(max(grepl(' ', geneSetNames)))
        stop('No spaces are allowed in any gene set name.')
    if(max(grepl('_', geneSetNames)))
        stop('No underscores are allowed in any gene set name.')
    extraNames <- setdiff(geneSetNames,
                          as.character(unique(scCol(scObj, labelCol))))
    if (length(extraNames))
        stop('All gene set names must exist in the ', labelCol,
             ' column of scObj.')
}

#' Remove single-valued columns from a data frame
#'
#' This function remove single-valued columns from a data frame.
#'
#' @param df A data frame.
#'
#' @return A data frame with any single-valued columns removed.
#'
#' @noRd
removeSVCols <- function(df){
    hasDistinctValues <- apply(df, 2, function(x) length(unique(x)) > 1)
    removedCol <- names(which(!hasDistinctValues))
    if(length(removedCol))
        message('Single-valued columns will be removed: ', removedCol, '.')
    df <- df[, hasDistinctValues]
    return(df)
}
