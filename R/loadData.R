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

genSummary_comb <- function(kept, trimmed, op, parameter)
{
  cbind(
    data.frame(op = op,
               parameter = parameter),
    genSummary_one(kept, 'k'),
    genSummary_one(trimmed, 't')
    )
}

genSummary_one <- function(seq_dat, names_pref)
{
  if (length(seq_dat) == 0){
    ret_val <-
    data.frame(seqs = 0,
               q25_length = NA,
               mean_length = NA,
               q75_length = NA,
               q25_qual = NA,
               mean_qual = NA,
               q75_qual = NA,
               prop_gaps = NA,
               prop_non_ACGT = NA,
               prop_AT = NA)
  } else {
    read_widths <- width(seq_dat)
    seq_qual <- seq_dat@quality
    per_read_quality <- letterFrequency(seq_qual@quality , names(encoding(seq_qual)))
    per_read_quality <- (per_read_quality %*% matrix(encoding(seq_qual), ncol=1))[,1]
    per_read_quality <- per_read_quality / read_widths
    conMat <- consensusMatrix(seq_dat@sread)
    tot_lets <- sum(conMat)
    prop_gaps <- round(sum(conMat['-',])/tot_lets,5)
    prop_non_ACGT <- round(sum(conMat[5:18,])/tot_lets,5)
    prop_AT <- round(sum(conMat[c('A','T'),])/tot_lets,5)

    ret_val <-
    data.frame(seqs = length(seq_dat),
               q25_length = round(quantile(read_widths, 0.25, names=FALSE),2),
               mean_length = round(mean(read_widths, na.rm=TRUE),2),
               q75_length = round(quantile(read_widths, 0.75, names=FALSE),2),
               q25_qual = round(quantile(per_read_quality, 0.25, names = FALSE),2),
               mean_qual = round(mean(per_read_quality, na.rm=T),2),
               q75_qual = round(quantile(per_read_quality, 0.75, names = FALSE),2),
               prop_gaps = round(prop_gaps,4),
               prop_non_ACGT = round(prop_non_ACGT,4),
               prop_AT = round(prop_AT,4))
  }
  names(ret_val) <- paste(names_pref, names(ret_val), sep = '_')
  return(ret_val)
}

genSummary.loadData <- function(result, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  # load data
  if (!op_args$cache)
  {
    stop('implement getKept now')
  } else 
  {
    seq_dat <- result$seq_dat
  }

  for (trim_step in result$trim_steps)
  {
    trim_dat <- result$metrics$per_read_metrics[,trim_step$name,drop=T]
    stopifnot(length(trim_dat) == length(seq_dat))
    if (length(trim_step$breaks) == 1)
    {
      stopifnot(all(trim_dat == trim_step$threshold))
      kept <- seq_dat
      trimmed <- seq_dat[0]
      parameter <- paste(trim_step$name, ' = ', trim_step$threshold, sep = '')
      summary_tab <- genSummary_comb(kept = kept,
                                     trimmed = trimmed,
                                     op = class(result),
                                     parameter = parameter)
    } else {
      stop('implement trimming with multiple categories now')
    }
  }
  result$summary <- summary_tab
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
