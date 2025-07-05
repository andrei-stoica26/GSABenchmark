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
#' @param dimMat UMAP dimensionality reduction matrix of the Seurat object. Unused
#' if normSilDF is NULL.
#' @param maxDist Maximum UMAP distance in the Seurat object. Unused in normSilDF
#' or dimMat is NULL
#' @param metrics Metrics used to evaluate the accuracy of the class boundary.
#' Must be a character vector containing at least two methods among "sensitivity",
#' "specificity", "precision", "accuracy", "sizeProximity" and "scoreSpecificity".
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
                                   maxDist = NULL,
                                   metrics = c('sensitivity', 'specificity',
                                               'precision','accuracy',
                                               'sizeProximity','scoreSpecificity')){
  supportedMetrics <- c('sensitivity', 'specificity', 'precision',
                      'accuracy', 'sizeProximity', 'scoreSpecificity')
  metrics <- checkMetrics(metrics, supportedMetrics)

  denseDF <- condenseRepeatedScores(df)

  truePos <- cumsum(denseDF[, 1])
  totalTrue <- truePos[nrow(denseDF)]
  if ('sensitivity' %in% metrics)
    df$sensitivity <- rep(truePos / totalTrue, denseDF[, 3])

  trueNeg <- revcumsum(denseDF[, 3] - denseDF[, 1])
  totalFalse <- nrow(df) - totalTrue
  if ('specificity' %in% metrics)
    df$specificity <- rep(trueNeg / totalFalse, denseDF[, 3])

  pos <- cumsum(denseDF[, 3])
  if ('precision' %in% metrics)
    df$precision <- rep(truePos / pos, denseDF[, 3])

  if ('accuracy' %in% metrics)
    df$accuracy <- rep((truePos + trueNeg) / nrow(df), denseDF[, 3])

  if ('sizeProximity' %in% metrics){
    maxDiff <- max(totalTrue, nrow(df) - totalTrue)
    df$sizeProximity <- rep(1 - abs(pos - totalTrue) / maxDiff, denseDF[, 3])
  }

  if ('scoreSpecificity' %in% metrics){
    score <- cumsum(denseDF[, 2] * denseDF[, 3])
    totalScore <- score[length(score)]
    df$scoreSpecificity <- rep(score / totalScore, denseDF[, 3])
  }

  if(!is.null(normSilDF)){
    label <- str_split(colnames(df)[2], '_')[[1]][2]
    sil <- normSilDF[label]
    sil <- sil[order(sil[, 1], decreasing=TRUE), drop=FALSE, ]

    silPos <- subset(sil, sil[, 1] > 0)
    maxSilWS <- sum(df[seq_len(nrow(silPos)), 2] * silPos[, 1])

    silNeg <- subset(sil, sil[, 1] < 0)
    silNeg <- silNeg[order(silNeg[, 1]), drop=FALSE, ]
    minSilWS <- sum(df[seq_len(nrow(silNeg)), 2] * silNeg[, 1])

    silDF <- data.frame(sil[rownames(df), 1], df[[2]])
    colnames(silDF) <- colnames(df)[c(1, 2)]
    denseSilDF <- condenseRepeatedScores(silDF)
    silScore <- cumsum(denseSilDF[, 1] * denseSilDF[, 2])

    df$silhouetteCoverage <- rep(sapply(silScore,
                                        function(x) liver::minmax(c(minSilWS, x, maxSilWS))[2]), denseSilDF[, 3])
    if (!is.null(dimMat) & !is.null(maxDist)){
      silCM <- centerOfMass(dimMat[rownames(silPos), ], silPos[, 1])
      centers <- centerOfMassV(dimMat[rownames(df), ], df[, 2])
      centers <- centers[findLastApps(df[, 2]), ]
      distances <- apply(centers, 1, function(x) euclidean(x, silCM))
      df$centrality <- rep(1 - distances / maxDist, denseDF[, 3])
    }
  }

  df <- computeMetricMeans(df, 3)
  return(df)
}

#' Compute the distribution metrics
#'
#' This function computes the distribution metrics
#'
#' @inheritParams computeBoundaryMetrics
#' @param metrics Metrics used to evaluate the accuracy of the class boundary.
#' Must be a character vector containing at least two methods among "AUC",
#' "Gini", "KS_Stat" and "PRAUC".
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
  df <- computeMetricMeans(df, 1)
  return(df)
}
