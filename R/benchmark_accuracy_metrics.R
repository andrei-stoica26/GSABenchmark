#' @importFrom abdiv cosine_distance euclidean
#' @importFrom spatstat.utils revcumsum
#' @importFrom MLmetrics AUC LogLoss PRAUC
#' @importFrom mltools mcc
#'
NULL

#' Compute the class boundary metrics
#'
#' This function computes the class boundary metrics
#'
#' @param df A data frame consisting of cells with a binary label column (1) and
#' a column containing the scores obtained by a gene set analysis method (2)
#' @return A benchmark data frame with cells as row names, labels in the first column,
#' gene set analysis method scores in the second column, and computed metric scores
#' in the remaining columns.
#'
#' @export
#'
computeBoundaryMetrics <- function(df){

  denseDF <- condenseRepeatedScores(df)

  truePos <- cumsum(denseDF[, 1])
  totalTrue <- truePos[nrow(denseDF)]
  df$sensitivity <- rep(truePos / totalTrue, denseDF[, 3])

  trueNeg <- revcumsum(denseDF[, 3] - denseDF[, 1])
  totalFalse <- nrow(df) - totalTrue
  df$specificity <- rep(trueNeg / totalFalse, denseDF[, 3])

  pos <- cumsum(denseDF[, 3])
  df$precision <- rep(truePos / pos, denseDF[, 3])

  df$accuracy <- rep((truePos + trueNeg) / nrow(df), denseDF[, 3])

  maxDiff <- max(totalTrue, nrow(df) - totalTrue)
  df$sizeProximity <- rep(1 - abs(pos - totalTrue) / maxDiff, denseDF[, 3])

  score <- cumsum(denseDF[, 2] * denseDF[, 3])
  totalScore <- score[length(score)]
  df$scoreCoverage <- rep(score / totalScore, denseDF[, 3])

  df <- computeMetricMeans(df, 3)
  return(df)
}

#' Compute the global evaluation metrics
#'
#' This function computes the global evaluation metrics
#'
#' @inheritParams computeBoundaryMetrics
#' @param normSilDF Data frame of normalized silhouettes
#' @param dimMat UMAP dimensionality reduction matrix of the Seurat object. Unused
#' if normSilDF is NULL.
#' @param maxDist Maximum UMAP distance in the Seurat object. Unused if normSilDF
#' or dimMat is NULL
#'
#' @return A benchmark data frame with one row and computed metric scores in the
#' columns.
#'
#' @export
#'
computeGlobalMetrics <- function(df, normSilDF = NULL, dimMat = NULL, maxDist = NULL){
  resDF <- data.frame(AUROC = AUC(df[, 2], df[, 1]),
                      PRAUC = PRAUC(df[, 2], df[, 1]),
                      rankLogScore = 1 - LogLoss(liver::minmax(rank(df[, 2])), df[, 1]),
                      labRankAlignment = rankAlignmentScore(df[, 1], df[, 2]))
  if(!is.null(normSilDF)){
    label <- str_split(colnames(df)[2], '_')[[1]][2]
    sil <- normSilDF[label]
    resDF$silRankAlignment <- rankAlignmentScore(sil[rownames(df), 1], df[, 2])
    if (!is.null(dimMat) & !is.null(maxDist)){
      silCM <- centerOfMass(dimMat, sil[rownames(dimMat), 1])
      scoreCM <- centerOfMass(dimMat, df[rownames(dimMat), 2])
      resDF$centrality <- 1 - euclidean(scoreCM, silCM) / maxDist
    }
  }
  resDF <- computeMetricMeans(resDF, 1)
  return(resDF)
}

#' Compute the Matthews correlation coefficient
#'
#' This function computes the Matthews correlation coefficient for an input
#' data frame generated through class boundary determination benchmarking. The
#' first column contains labels while the second is binarized based on the previously
#' determined cutoff (the first element in the second column).
#'
#' @param boundaryResDF
#'
#' @return The Matthews correlation coefficient
#'
#' @export
#'
computeMCC <- function(boundaryResDF){
  boundaryResDF <- boundaryResDF[, c(1, 2)]
  boundaryResDF$prediction <- 0
  boundaryResDF[which(boundaryResDF[, 2] >= boundaryResDF[1, 2]), 3] <- 1
  return(mltools::mcc(boundaryResDF[, 1], boundaryResDF[, 3]))
}
