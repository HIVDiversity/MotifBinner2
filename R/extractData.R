#' Extracts required dataset by navigating down into a list
#' @inheritParams applyOperation
#' @export

extractData <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)

  data_set <- all_results[[data_source_indx]][[op_args$extract_levels[1]]]
  for (i in 2:length(op_args$extract_levels)){
    data_set <- data_set[[op_args$extract_levels[2]]]
  }
  seq_dat <- data_set

  per_read_metrics <- data.frame('read_exists' = rep(1, length(seq_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'extractData'
  if (op_args$cache){
    result$seq_dat <- seq_dat
  }
  result$input_dat <- seq_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

computeMetrics.extractData <- function(result, config, seq_dat)
{
  return(result)
}

saveToDisk.extractData <- function(result, config, seq_dat)
{
  return(result)
}

print.extractData <- function(result, config)
{
  cat('\n-------------------\n')
  cat(  'Operation: extractData\n')
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
