#' @importFrom abdiv euclidean
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

  labelSumsPerScore <- denseDF[, 1]
  scoreThresholds <- denseDF[, 2]
  frequencies <- denseDF[, 3]
  nThresholds <- nrow(denseDF)

  Total <- nrow(df)
  TP <- cumsum(labelSumsPerScore)
  FP <- cumsum(frequencies) - TP

  True <- TP[nThresholds]
  False <- Total - True
  TN <- False - FP

  Positive <- TP + FP
  Largest <- max(True, False)

  df$sensitivity <- rep(TP / True, frequencies)
  df$specificity <- rep(TN / False, frequencies)
  df$precision <- rep(TP / Positive, frequencies)
  df$accuracy <- rep((TP + TN) / Total, frequencies)
  df$sizeProximity <- rep(1 - abs(Positive - True) / Largest, frequencies)

  score <- cumsum(scoreThresholds * frequencies)
  totalScore <- score[length(score)]
  if (totalScore > 0)
    df$scoreCoverage <- rep(score / totalScore, frequencies) else
      df$scoreCoverage <- rep(0, frequencies)

  df <- computeMetricMeans(df, 3)
  return(df)
}

#' Compute the MCC at each threshold
#'
#' This function computes the MCC at each threshold
#'
#' @inheritParams computeBoundaryMetrics
#'
#' @return A benchmark data frame with cells as row names, labels in the first column,
#' gene set analysis method scores in the second column, and MCC scores in the third
#' column
#'
#' @export
#'
computeMCCMetric <- function(df){
  denseDF <- condenseRepeatedScores(df)

  labelSumsPerScore <- denseDF[, 1]
  scoreThresholds <- denseDF[, 2]
  frequencies <- denseDF[, 3]
  nThresholds <- nrow(denseDF)

  Total <- nrow(df)
  Positives <- cumsum(frequencies)
  Negatives <- Total - Positives
  TP <- cumsum(labelSumsPerScore)
  FP <- Positives - TP

  True <- TP[nThresholds]
  False <- Total - True
  classGeomMean <- sqrt(True * False)

  denseDF$x <- TP * False - FP * True
  denseDF$y <- classGeomMean * sqrt(Positives * Negatives)
  denseDF$y[nrow(denseDF)] <- -1

  thresholdMCC <- denseDF$x / denseDF$y
  df$MCC <- rep(thresholdMCC, frequencies)
  df <- df[order(df$MCC, decreasing=TRUE),]

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
                      labRankAlignment = rankAlignmentScore(df[, 1], df[, 2]))
  if(!is.null(normSilDF)){
    label <- str_split(colnames(df)[2], '_')[[1]][2]
    sil <- normSilDF[label]
    resDF$silRankAlignment <- rankAlignmentScore(sil[rownames(df), 1], df[, 2])
    if (!is.null(dimMat) & !is.null(maxDist)){
      scoreCM <- centerOfMass(dimMat, df[rownames(dimMat), 2])
      if (length(intersect(scoreCM, NaN)))
        resDF$centrality <- 0 else{
          silCM <- centerOfMass(dimMat, sil[rownames(dimMat), 1])
          resDF$centrality <- 1 - euclidean(scoreCM, silCM) / maxDist
      }
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
#' @param boundaryResDF A data frame generated with the boundary benchmark
#' @param threshold Threshold that determines class boundary
#'
#' @return The Matthews correlation coefficient
#'
#' @export
#'
computeMCC <- function(boundaryResDF, threshold = boundaryResDF[1, 2]){
  boundaryResDF <- boundaryResDF[, c(1, 2)]
  boundaryResDF$prediction <- 0
  boundaryResDF$prediction[which(boundaryResDF[, 2] >= threshold)] <- 1
  return(mltools::mcc(boundaryResDF[, 1], boundaryResDF[, 3]))
}
