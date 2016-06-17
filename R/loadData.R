#' Loads the fwd and rev fastq reads
#' @inheritParams applyOperation
#' @export

loadData <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  input_file <- config$operation_list[[op_number]]$input_file
  if (is.null(input_file)){
    stop('Input file must be specified')
  }
  if (!file.exists(input_file)){
    stop(paste('input file does not exists: ', input_file, sep = ''))
  }

  seq_dat <- readFastq(input_file)
  per_read_metrics <- data.frame('read_exists' = rep(1, length(seq_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'loadData'
  if (op_args$cache){
    result$seq_dat <- seq_dat
  }
  return(result)
}

computeMetrics.loadData <- function(result, config)
{
  return(result)
}

saveToDisk.loadData <- function(result, config)
{
  return(result)
}

print.loadData <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: loadData')
  cat('\n-------------------')
  cat('\nLoaded Sequences:\n')
  print(result$summary[,c('parameters', 'seqs_kept', 'mean_length_kept', 'mean_qual_kept')])
  return(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
