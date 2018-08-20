#' Prepares config list for inclusion in final report
#'
#' Converts the configuration used to call MotifBinner to a data.frame that can
#' be included in the final report.
#'
#' The configuration used to specify which opertation should be applied and
#' what settings to use is contained in a nested list of lists. These lists are
#' typically built using the buildConfig function. By including these settings
#' and the version of MotifBinner in the final report an easy way to document
#' the exact procedures that were used to produce a dataset.
#'
#' The main action performed by this prepConfig operation is just to convert
#' the nested list of lists that holds the config to a data.frame and to
#' provide a convenient way, using the standard approach of MotifBinner, to add
#' this information into the final report.
#'
#' @inheritParams applyOperation
#' @export

prepConfig <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)


#  input_file <- config$operation_list[[op_number]]$data_source
#  if (is.null(input_file)){
#    stop('Input file must be specified')
#  }
#  if (!file.exists(input_file)){
#    stop(paste('input file does not exists: ', input_file, sep = ''))
#  }

  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = as.data.frame(config))
  class(result) <- 'prepConfig'
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

computeMetrics.prepConfig <- function(result, config, seq_dat)
{
  return(result)
}

saveToDisk.prepConfig <- function(result, config, seq_dat)
{
  return(result)
}

print.prepConfig <- function(result, config)
{
  cat('\n-------------------\n')
  cat(  'Operation: prepConfig\n')
  cat(  '-------------------\n')
#  cat('\nLoaded Sequences: Bases:\n')
#  print(result$seq_dat@sread)
#  cat('\nLoaded Sequences: Qualities:\n')
#  print(result$seq_dat@quality@quality)
#  cat('\nSummary:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  invisible(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
