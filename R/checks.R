#' Check the input metrics
#'
#' This function checks the input metrics agaisnt supported metrics
#'
#' @param metrics Benchmark metrics
#' @param supportedMetrics Supported metrics
#'
#' @return Metrics after filtering out any errors in the input
#'
checkMetrics <- function(metrics, supportedMetrics){
  uniqueMetrics <- unique(metrics)
  if (length(metrics) > length(uniqueMetrics)){
    warning('The metrics vector contains non-unique values. The duplicates will be ignored')
    metrics <- uniqueMetrics
  }
  unknownMetrics <- setdiff(metrics, supportedMetrics)
  if(length(unknownMetrics) > 0){
    warning('The metrics vector contains unknown metric names. These will be ignored')
    metrics <- intersect(metrics, supportedMetrics)
  }
  if(length(metrics) < 2)
    stop('Not enough valid metrics provided; including at least 2 is required')
  return(metrics)
}
