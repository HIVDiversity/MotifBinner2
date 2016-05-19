#' Save the relevant results of an operation to disk
#' @param result The result produced from a specific operation
#' @inheritParams applyOperation
#' @export

saveToDisk <- function(result, config)
{
  UseMethod('saveToDisk')
}

#' Generates a summary for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genSummary <- function(result, config)
{
  UseMethod('genSummary')
}


