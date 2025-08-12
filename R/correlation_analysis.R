#' Calculate correlation matrix for method results
#'
#' This function calculates the correlation matrix for all the methods for an
#' identity class.
#'
#' @inheritParams mdsScoreSummary
#' @param corMethod Correlation method.
#'
#' @return A list of correlation matrices.
#'
#' @export
#'
corrSummary <- function(scObj, smr, corMethod = 'spearman'){

    gsaMethods <- sort(rownames(smr$boundary[[1]]))
    geneSetNames <- colnames(smr$boundary[[1]])
    geneSetNames <- geneSetNames[seq(length(geneSetNames) - 1)]

    doAgg <- TRUE

    corrDFs <- lapply(geneSetNames, function(gsName){
        setDF <- as.matrix(metadataDF(scObj)[, paste0(gsaMethods,
                                                      '_', gsName)])
        setDF <- removeSVCols(setDF)
        colnames(setDF) <- str_remove(colnames(setDF), paste0('_', gsName))

        if(length(setdiff(gsaMethods, colnames(setDF)))){
            message('At least one single-valued columns has been removed.',
                    ' Aggregate correlation will not be computed.')
            doAgg <- FALSE
        }

        return(round(cor(setDF, method=corMethod), 2))
    })

    if(doAgg){
        aggCorr <- round(Reduce(`+`, corrDFs) / length(corrDFs), 2)
        corrDFs <- c(corrDFs, list(aggCorr))
        names(corrDFs) <- c(geneSetNames, 'aggregate')
    }

    return(corrDFs)
}
