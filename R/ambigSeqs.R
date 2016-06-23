#' Get the data kept after a trimming step
#' @export

getKept <- function(result, seq_dat=NULL)
{
  if (is.null(seq_dat) & (!('input_dat' %in% names(result))))
  {
    stop('data to be trimmed must be supplied either via the seq_dat argument or as an element named input_dat in the result list')
  }

  trim_criteria <- matrix(TRUE, nrow = length(seq_dat), ncol = length(result$trim_steps))
  i <- 0
  for (trim_step in result$trim_steps)
  {
    i <- i+1
    trim_dat <- result$metrics$per_read_metrics[,trim_step$name,drop=T]
    stopifnot(length(trim_dat) == length(seq_dat))
    if (length(trim_step$breaks) == 1)
    {
      stopifnot(all(trim_dat == trim_step$threshold))
      kept <- seq_dat
      trimmed <- seq_dat[0]
    } else {
      if ('comparator' %in% names(trim_step))
      {
        comparator = trim_step$comparator
      } else {
        comparator = `<=`
      }
      trim_criteria[,i] <- comparator(trim_dat, trim_step$threshold)
    }
  }
  seq_dat <- seq_dat[apply(trim_criteria, 1, all)]
  return(seq_dat)
}

#' Get the data trimmed at a trimming step
#' @export

getTrimmed <- function(seq_dat, kept_dat)
{
  seq_dat[!(as.character(seq_dat@id) %in% as.character(kept_dat@id))]
}

#' Removes sequences with too many ambig bases
#' @inheritParams applyOperation
#' @export

ambigSeqs <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat
  
  per_read_metrics <- ambigSeqs_internal(seq_dat)

  threshold <- op_args$threshold
  if (is.null(threshold))
  {
    threshold <- 0.02
  }
  trim_steps <- list(step1 = list(name = 'perc_ambig',
                                  threshold = threshold,
                                  breaks = c(-Inf, 0, 0.01, 0.02, 0.1, 0.5, 1)
                                  )
                    )

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- getKept(result, seq_dat=seq_dat)
  class(result) <- 'ambigSeqs'
  if (op_args$cache){
    result$seq_dat <- kept_dat
  }
  result$input_dat <- seq_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

ambigSeqs_internal <- function(seq_dat)
{
  counts <- alphabetFrequency(seq_dat@sread)
  tmp <- data.frame(counts)
  names(tmp) <- paste('c',attr(counts, 'dimnames')[[2]],sep='')
  counts <- tmp
  rm(tmp)

  ambigCols <- !(gsub('^c','', names(counts)) %in% c('A','C','G','T','-'))
  counts$seq_len <- width(seq_dat@sread)
  counts$ambig <- apply(counts[,ambigCols], 1, sum)
  counts$perc_ambig <- counts$ambig/counts$seq_len
  per_read_metrics <- counts[,c('seq_len', 'ambig', 'perc_ambig')]
  return(per_read_metrics)
}

saveToDisk.ambigSeqs <- function(result, config, seq_dat)
{
  kept <- getKept(result, seq_dat)
  trimmed <- getTrimmed(seq_dat = seq_dat, kept_dat = kept)

  if (length(kept) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_kept_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(kept, tmp_name, compress=F)
  }
  if (length(trimmed) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_trimmed_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(trimmed, tmp_name, compress=F)
  }
  return(result)
}

computeMetrics.ambigSeqs <- function(result, config, seq_dat)
{
  return(result)
}

print.ambigSeqs <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: ambigSeqs')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
