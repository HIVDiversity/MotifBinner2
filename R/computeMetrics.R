#' Computes metrics useful for reporting
#' @param result The result produced from a specific operation
#' @inheritParams applyOperation
#' @export

computeMetrics <- function(result, config, seq_dat)
{
  UseMethod('computeMetrics')
}

