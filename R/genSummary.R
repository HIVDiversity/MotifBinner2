#' Generates a summary for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genSummary <- function(result, config)
{
  UseMethod('genSummary')
}


