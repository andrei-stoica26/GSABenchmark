#' @importFrom abdiv cosine_distance euclidean
#' @importFrom spatstat.utils revcumsum
#'
NULL

#' Add the cell count based-metrics to the benchmark
#'
#' This function adds the cell count based-based metrics to the benchmark
#'
#' @param df A benchmark data frame
#'
#' @return A benchmark data frame with added cell count-based metrics
#'
addCellCountMetrics <- function(df){
  truePos <- cumsum(df$label)
  totalTrue <- truePos[nrow(df)]
  df$sensitivity <- truePos / totalTrue

  trueNeg <- revcumsum(1 - df$label)
  totalFalse <- nrow(df) - totalTrue
  df$specificity <- trueNeg / totalFalse

  pos <- seq_len(nrow(df))
  df$precision <- truePos / pos

  df$accuracy <- (truePos + trueNeg) / nrow(df)

  maxDiff <- max(totalTrue, nrow(df) - totalTrue)
  df$sizeProximity <- 1 - abs(pos - totalTrue) / maxDiff
  return(df)
}

#' Add the score-based metrics to the benchmark
#'
#' This function adds the score-based metrics to the benchmark
#'
#' @param df A benchmark data frame
#'
#' @return A benchmark data frame with added score-based metrics
#'
addScoreMetric <- function(df){
  score <- cumsum(df[, 2])
  totalScore <- score[length(score)]
  df$scoreSpecificity <- score / totalScore
  return(df)
}

#' Add the silhouette-based metrics to the benchmark
#'
#' This function adds the silhouette-based metrics to the benchmark
#'
#' @param df A benchmark data frame
#' @param normSilDF Data frame of normalized silhouettes
#' @param dimMat UMAP dimensionality reduction matrix of the Seurat object
#'
#' @return A benchmark data frame with added silhouette-based metrics
#'
addCentralityMetrics <- function(df, normSilDF, dimMat){
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

  maxDist <- max(dist(dimMat))
  silCM <- centerOfMass(dimMat[rownames(silPos), ], silPos[, 1])

  centers <- centerOfMassV(dimMat[rownames(df), ], df[, 2])
  distances <- apply(centers, 1, function(x) euclidean(x, silCM))
  df$centrality <- 1 - distances / maxDist

  return(df)
}
