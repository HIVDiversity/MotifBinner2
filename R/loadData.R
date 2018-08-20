#' Loads a fastq file
#'
#' Given a file name, it calls ShortRead::readFastq and populates the result list
#'
#' The action function for the class loadData
#'
#' @inheritParams applyOperation
#' @export

loadData <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)


#  input_file <- config$operation_list[[op_number]]$data_source
  input_file <- op_args$data_source
  if (is.null(input_file)){
    stop('Input file must be specified')
  }
  if (!file.exists(input_file)){
    stop(paste('input file does not exists: ', input_file, sep = ''))
  }

  max_seq <- op_args$max_seq
  seq_dat <- readFastq(input_file)
  total_seq <- length(seq_dat)
  if (!is.null(max_seq)){
    if (max_seq > 0 & length(seq_dat) > max_seq){
      seq_dat <- seq_dat[1:max_seq]
    }
  }
  per_read_metrics <- data.frame('read_exists' = rep(1, length(seq_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics, total_seq = total_seq))
  class(result) <- 'loadData'
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

computeMetrics.loadData <- function(result, config, seq_dat)
{
  return(result)
}

saveToDisk.loadData <- function(result, config, seq_dat)
{
  return(result)
}

print.loadData <- function(result, config)
{
  cat('\n-------------------\n')
  cat(  'Operation: loadData\n')
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
