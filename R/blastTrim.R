#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

blastTrim <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  stopifnot(length(op_args$data_source) == 1)
  stopifnot(op_args$cache)

  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  
  seq_dat <- all_results[[data_source_indx]]$seq_dat
  seq_dat_fasta <- seq_dat@sread
  names(seq_dat_fasta) <- seq_dat@id

  seq_file_name <- file.path(op_dir, 'seq.fasta')
  writeXStringSet(seq_dat_fasta,
             seq_file_name,
             compress=FALSE)

  blast_hits_file_name <- file.path(op_dir, 'blast_hits.tab')
  blast_command <- 
         paste('blastn', 
               ' -db ', op_args$database, 
               ' -query ', seq_file_name,
               ' -out ', blast_hits_file_name,
               ' -num_threads ', config$ncpu,
               ' -outfmt 6',
               ' -max_hsps 1',
               ' -max_target_seqs 1', 
               sep = '')
  print(blast_command)
  system(blast_command)
  
  per_read_metrics <- try(read.delim(blast_hits_file_name, header = F,
                                     stringsAsFactors = F)[,c(1,2,3),],
                          silent = TRUE)
  if (class(per_read_metrics) == 'try-error'){
    if (grepl("no lines av", attr(per_read_metrics, "condition"))){
      per_read_metrics <- data.frame(query = character(0),
                                     target = character(0),
                                     pident = numeric(0))
      
      no_hits <- as.character(seq_dat@id)
    } else {
      print(per_read_metrics)
      stop('something failed while blasting')
    }
  } else {
    names(per_read_metrics) <- c('query', 'target', 'pident')
    no_hits <- as.character(seq_dat@id)[!(as.character(seq_dat@id) %in% per_read_metrics$query)]
  }
  file.remove(seq_file_name)
  if (length(no_hits) > 0){
    per_read_metrics <- rbind(as.data.frame(per_read_metrics),
      data.frame(query = no_hits,
                 target = NA,
                 pident = 0,
                 stringsAsFactors = F))
  }
  per_read_metrics <- per_read_metrics[order(per_read_metrics$query),]
  seq_dat <- seq_dat[order(as.character(seq_dat@id))]

  if (op_args$action == 'trim_hits') {
    trim_steps <- list(step1 = list(name = 'pident',
                                    threshold = op_args$threshold,
                                    comparator = `<=`,
                                    breaks = c(-Inf, op_args$threshold - 5,
                                               op_args$threshold,
                                               op_args$threshold + 5, Inf)))
  } else {
    trim_steps <- list(step1 = list(name = 'pident',
                                    threshold = op_args$threshold,
                                    comparator = `>=`,
                                    breaks = c(Inf, op_args$threshold + 5,
                                               op_args$threshold,
                                               op_args$threshold - 5, -Inf)))
  }

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'blastTrim'
  result$seq_dat <- getKept(result, seq_dat)
  result$input_dat <- seq_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.blastTrim <- function(result, config, seq_dat = NULL)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')
  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)

  seq_file_name <- paste(config$base_for_names, '_kept_', result$config$op_args$name, '.fastq', sep = '')
  writeFastq(result$seq_dat, file.path(op_dir, seq_file_name), compress = F)

  trim_file_name <- paste(config$base_for_names, '_trimmed_', result$config$op_args$name, '.fastq', sep = '')
  trim_dat <- getTrimmed(seq_dat, result$seq_dat)
  writeFastq(trim_dat, file.path(op_dir, trim_file_name), compress = F)

  clustering_file <- paste(config$base_for_names, '_blast_hits_', result$config$op_args$name, '.csv', sep = '')
  write.csv(result$metrics$per_read_metrics, file.path(op_dir, clustering_file), row.names=F)
  return(result)
}

computeMetrics.blastTrim <- function(result, config, seq_dat)
{
  return(result)
}

print.blastTrim <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: blastTrim')
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
