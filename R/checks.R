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
        stop('No spaces are allowed in any gene set name')
    extraNames <- setdiff(geneSetNames,
                          as.character(unique(scCol(scObj, labelCol))))
    if (length(extraNames))
        stop('All gene set names must exist in the ', labelCol,
             ' column of scObj.')
}
