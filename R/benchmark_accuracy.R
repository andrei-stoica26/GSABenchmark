#' @importFrom abdiv cosine_distance euclidean
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
  nLabel <- cumsum(df$label)
  totalCells <- nrow(df)
  nCells <- seq_len(totalCells)
  df$specificity <- nLabel / nCells * 100

  totalLabel <- nLabel[length(nLabel)]
  df$coverage <- nLabel / totalLabel * 100

  maxDiff <- max(totalLabel, totalCells - totalLabel)
  df$sizeProximity <- (1 - abs(nCells - totalLabel) / maxDiff) * 100
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
  df$scoreCoverage <- score / totalScore * 100
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

  df$silhouetteCoverage <- sapply(silScore, function(x) liver::minmax(c(minSilWS, x, maxSilWS))[2] * 100)

  maxDist <- max(dist(dimMat))
  silCM <- centerOfMass(dimMat[rownames(silPos), ], silPos[, 1])

  centers <- centerOfMassV(dimMat[rownames(df), ], df[, 2])
  distances <- apply(centers, 1, function(x) euclidean(x, silCM))
  df$centrality <- (1 - distances / maxDist) * 100

  return(df)
}
