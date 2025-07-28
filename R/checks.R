#' Check if any gene set names are not found in the identity class column
#'
#' This function check if any gene set names are not found in the identity
#' class column.
#'
#' @inheritParams runBenchmark
#'
#' @return None. This function is called for its side effect
#' (showing an error message).
#'
#' @noRd
#'
checkSetNames <- function(scObj, labelCol, geneSetNames){
    if (is(labelCol)[1] != 'character')
        stop('geneSetNames must be a character vector.')
    extraNames <- setdiff(geneSetNames,
                          as.character(unique(scCol(scObj, labelCol))))
    if (length(extraNames))
        stop('All gene set names must exist in the ', labelCol,
             ' column of scObj.')
}
