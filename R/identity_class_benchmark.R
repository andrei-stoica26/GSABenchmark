#' @importFrom fabR silently_run
#' @importFrom logr log_close log_open log_print
NULL


#' Run gene set analysis methods
#'
#' This function runs the gene set analysis methods.
#'
#' @param seuratObj A Seurat object
#' @param markerList A list of gene sets
#' @param markerNames The names of the gene sets
#' @param gsaMethods List of gene set analysis methods
#' @param logFile Name of the log file saving running times
#'
#' @return A Seurat object with the results of the runs stored as metadata columns
#'
#' @export
#'
runGSAMethods <- function(seuratObj, markerList, markerNames, gsaMethods, logFile = 'gsa'){
  log_open(logFile, logdir=FALSE, show_notes=FALSE, compact=TRUE, header=FALSE)
  for (i in seq_along(markerList)){
    setName <- markerNames[i]
    for (j in seq_along(gsaMethods)){
      method <- gsaMethods[j]
      runningMsg <- paste0('Running ', method, ' for ', setName, ' genes...')
      message(runningMsg)
      log_print(runningMsg, console=F)
      fun <- eval(as.name(paste0('run', method)))
      colStr <- paste0(method, '_', setName)
      seuratObj <- silently_run(timeCode(fun, logPrint=TRUE, seuratObj, markerList[[i]], colStr))
    }
  }
  log_close(footer=F)
  return(seuratObj)
}

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

#' Add the silhouette-based metric to the benchmark
#'
#' This function adds the silhouette-based metrics to the benchmark
#'
#' @param df A benchmark data frame
#' @param normSilDF Data frame of normalized silhouettes
#'
#' @return A benchmark data frame with added silhouette-based metrics
#'
addCentralityMetric <- function(df, normSilDF){
  label <- str_split(colnames(df)[2], '_')[[1]][2]
  sil <- normSilDF[label]
  sil <- sil[order(sil[, 1], decreasing=TRUE), drop = FALSE, ]

  silPos <- sil[sil >= 0]
  maxSilWS <- sqrt(sum(df[seq_along(silPos), 2], silPos))

  silNeg <- rev(sil[sil < 0])
  minSilWS <- sqrt2(sum(df[seq_along(silNeg), 2], silNeg))

  sil <- sil[rownames(df), drop = F, ]
  silScore <- sapply(cumsum(df[, 2] * sil[, 1]), sqrt2)

  df$centrality <- sapply(silScore, function(x) liver::minmax(c(minSilWS, x, maxSilWS))[2] * 100)
  return(df)
}

#' Add the score-based metrics to the benchmark
#'
#' This function adds the score-based metrics to the benchmark
#'
#' @param seuratObj A Seurat object
#' @param labelCol The Seurat metadata column containing the ground truth annotation
#' @param scoreCol The Seurat metadata column containing the gene set analysis method score
#' @param label The identity assessed from the labelCol column
#' @param nMetrics Number of metrics
#' @inheritParams addCentralityMetric
#'
#' @return A benchmark data frame
#'
#' @export
#'
identityClassMatch <- function(seuratObj, labelCol, scoreCol, normSilDF, label, nMetrics = 5){
  df <- seuratObj@meta.data[, c(labelCol, scoreCol)]
  df$label <- as.integer(df[, 1] %in% label)
  df <- df[order(df[, scoreCol], decreasing=TRUE),]

  df <- addCellCountMetrics(df)
  df <- addScoreMetric(df)
  df <- addCentralityMetric(df, normSilDF)

  df$accuracy <- rowMeans(df[, 2 + seq_len(nMetrics)])
  df <- df[order(df$accuracy, decreasing=TRUE),]
  return(df)
}

#' Performs a comparison of different gene set analysis methods using annotations
#'
#' This function compares the ability of different gene set analysis methods to
#' correctly identity cells with a certain label.
#'
#' @param seuratObj A Seurat object
#' @param labelCol The column of labels
#' @param markerNames The names of marker sets
#' @param gsaMethods The assessed gene set analysis methods
#' @param labels The assessed labels. Must be a subset of the values in labelCol
#' @inheritParams identityClassMatch
#'
#' @return A benchmark data frame
#'
#' @export
#'
identityClassBenchmark <- function(seuratObj, labelCol, markerNames, gsaMethods, normSilDF, labels = markerNames){
  res <- lapply(seq_along(markerNames), function(i) {
    setName <- markerNames[i]
    setRes <- lapply(seq_along(gsaMethods), function(j){
      method <- gsaMethods[j]
      scoreCol <- paste0(method, '_', setName)
      df <- identityClassMatch(seuratObj, labelCol, scoreCol, normSilDF, labels[i])
      return(df)
    })
    names(setRes) <- gsaMethods
    return(setRes)
  })
  names(res) <- labels
  return(res)
}

#' Compares the identification of different labels through different gene set
#' analysis methods
#'
#' This function compares the identification of different labels through different gene set
#' analysis methods.
#'
#' @param icBenchmark A list of benchmark data frames
#' @inheritParams identityClassMatch
#'
#' @return Summary data frames
#'
#' @export
#'
identityClassBenchmarkSummary <- function(icBenchmark, nMetrics = 5){
  markerNames <- names(icBenchmark)
  gsaMethods <- names(icBenchmark[[1]])
  metrics <- colnames(icBenchmark[[1]][[1]])[2 + seq_len(nMetrics + 1)]
  smr <- lapply(metrics, function(metric){
    df <- data.frame(Reduce(rbind, lapply(gsaMethods, function(method)
      sapply(markerNames, function(setName)
        icBenchmark[[setName]][[method]][[metric]][1]))))
    rownames(df) <- gsaMethods
    df$avg <- rowMeans(df)
    df <- df[order(df$avg, decreasing=TRUE), ]
    return(df)
  })
  names(smr) <- metrics
  df <- data.frame(lapply(smr, function(x) x[gsaMethods, ]$avg))
  rownames(df) <- gsaMethods
  df <- df[order(df$accuracy, decreasing=TRUE), ]
  smr <- c(smr, list(total = df))
  return(smr)
}
