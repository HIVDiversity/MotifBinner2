getKept <- function(result, seq_dat=NULL)
{
  if (is.null(seq_dat) & (!('input_dat' %in% names(result))))
  {
    stop('data to be trimmed must be supplied either via the seq_dat argument or as an element named input_dat in the result list')
  }

  if (length(result$trim_steps) > 1){stop('implement multiple trimming columns now')}
  for (trim_step in result$trim_steps)
  {
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
      seq_dat <- seq_dat[comparator(trim_dat, trim_step$threshold)]
    }
  }
  return(seq_dat)
}

getTrimmed <- function(seq_dat, kept_dat)
{
  seq_dat[!(seq_dat@id %in% kept_dat@id)]
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

saveToDisk.ambigSeqs <- function(result, config)
{
  for (data_set_name in names(result$kept)){
    seq_dat <- result$kept[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      writeFastq(seq_dat, file.path(result$op_dir, 
        paste(config$prefix_for_names, '_kept_', data_set_name, '.fastq', sep = '')), compress=F)
    }
  }
  for (data_set_name in names(result$trimmed)){
    seq_dat <- result$trimmed[[data_set_name]]
    if (length(seq_dat) > 0)
    {
      writeFastq(seq_dat, file.path(result$op_dir, 
        paste(config$prefix_for_names, '_trimmed_', data_set_name, '.fastq', sep = '')), compress=F)
    }
  }
  return(result)
}

genSummary.ambigSeqs <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'ambigSeqs',
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$kept$fwd_reads,
                        trimmed_seq_dat = result$trimmed$fwd_reads),
    genSummary_internal(operation = 'ambigSeqs',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$kept$rev_reads,
                        trimmed_seq_dat = result$trimmed$rev_reads))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$op_dir, 'ambigSeqs_summary.csv'), row.names=FALSE)
  return(result)
}

computeMetrics.ambigSeqs <- function(result, config)
{
  return(result)
}

print.ambigSeqs <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: ambigSeqs')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameters', 'seqs_kept', 'mean_length_kept', 'mean_qual_kept')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameters', 'seqs_trimmed', 'mean_length_trimmed', 'mean_qual_trimmed')])
  return(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
