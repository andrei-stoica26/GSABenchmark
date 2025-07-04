#' @importFrom abdiv cosine_distance euclidean
#' @importFrom spatstat.utils revcumsum
#' @importFrom MLmetrics AUC Gini KS_Stat PRAUC
#'
NULL

#' Compute the class boundary metrics
#'
#' This function computes the class boundary metrics
#'
#' @param df A data frame consisting of cells with a binary label column (1) and
#' a column containing the scores obtained by a gene set analysis method (2)
#' @param normSilDF Data frame of normalized silhouettes
#' @param dimMat UMAP dimensionality reduction matrix of the Seurat object
#' @param metrics Metrics used to evaluate the accuracy of the class boundary.
#' Must be a character vector containing at least two methods among 'sensitivity',
#' 'specificity', 'precision', 'accuracy', 'sizeProximity' and 'scoreSpecificity'.
#' Additionally, silhouette coverage will be automatically computed if normSilDF
#' is provided, and centrality will be automatically computed if dimMat is also
#' provided.
#' @return A benchmark data frame with cells as row names, labels in the first column,
#' gene set analysis method scores in the second column, and computed metric scores
#' in the remaining columns.
#'
#' @export
#'
computeBoundaryMetrics <- function(df,
                                   normSilDF = NULL,
                                   dimMat = NULL,
                                   metrics = c('sensitivity', 'specificity',
                                               'precision','accuracy',
                                               'sizeProximity','scoreSpecificity')){
  supportedMetrics <- c('sensitivity', 'specificity', 'precision',
                      'accuracy', 'sizeProximity', 'scoreSpecificity')
  metrics <- checkMetrics(metrics, supportedMetrics)

  truePos <- cumsum(df[, 1])
  totalTrue <- truePos[nrow(df)]
  if ('sensitivity' %in% metrics)
    df$sensitivity <- truePos / totalTrue

  trueNeg <- revcumsum(1 - df[, 1])
  totalFalse <- nrow(df) - totalTrue
  if ('specificity' %in% metrics)
    df$specificity <- trueNeg / totalFalse

  pos <- seq_len(nrow(df))
  if ('precision' %in% metrics)
    df$precision <- truePos / pos

  if ('accuracy' %in% metrics)
    df$accuracy <- (truePos + trueNeg) / nrow(df)

  if ('sizeProximity' %in% metrics){
    maxDiff <- max(totalTrue, nrow(df) - totalTrue)
    df$sizeProximity <- 1 - abs(pos - totalTrue) / maxDiff
  }

  if ('scoreSpecificity' %in% metrics){
    score <- cumsum(df[, 2])
    totalScore <- score[length(score)]
    df$scoreSpecificity <- score / totalScore
  }

  if(!is.null(normSilDF)){
    label <- str_split(colnames(df)[2], '_')[[1]][2]
    sil <- normSilDF[label]
    sil <- sil[order(sil[, 1], decreasing=TRUE), drop=FALSE, ]

    silPos <- subset(sil, sil[, 1] >= 0)
    maxSilWS <- sqrt(sum(df[seq_len(nrow(silPos)), 2], silPos[, 1]))

    silNeg <- subset(sil, sil[, 1] < 0)
    silNeg <- silNeg[order(silNeg[, 1]), drop=FALSE, ]
    minSilWS <- sqrt2(sum(df[seq_len(nrow(silNeg)), 2], silNeg[, 1]))

    sil <- sil[rownames(df), drop = F, ]
    silScore <- sapply(cumsum(df[, 2] * sil[, 1]), sqrt2)

    df$silhouetteCoverage <- sapply(silScore, function(x) liver::minmax(c(minSilWS, x, maxSilWS))[2])

    if (!is.null(dimMat)){
      maxDist <- max(dist(dimMat))
      silCM <- centerOfMass(dimMat[rownames(silPos), ], silPos[, 1])

      centers <- centerOfMassV(dimMat[rownames(df), ], df[, 2])
      distances <- apply(centers, 1, function(x) euclidean(x, silCM))
      df$centrality <- 1 - distances / maxDist
    }
  }

  df <- computeMetricOveralls(df, 3)
  return(df)
}

#' Compute the distribution metrics
#'
#' This function computes the distribution metrics
#'
#' @inheritParams computeBoundaryMetrics
#' @param metrics Metrics used to evaluate the accuracy of the class boundary.
#' Must be a character vector containing at least two methods among 'AUC',
#' 'Gini', 'KS_Stat' and 'PRAUC'.
#'
#' @return A benchmark data frame with one row and computed metric scores in the
#' columns.
#'
#' @export
#'
computeDistributionMetrics <- function(df, metrics = c('AUC', 'Gini', 'KS_Stat', 'PRAUC')){
  supportedMetrics <- c('AUC', 'Gini', 'KS_Stat', 'PRAUC')
  metrics <- checkMetrics(metrics, supportedMetrics)
  df <- data.frame(lapply(metrics, function(x) do.call(x, list(df[, 2], df[, 1]))))
  colnames(df) <- metrics
  df <- computeMetricOveralls(df, 1)
  return(df)
}
