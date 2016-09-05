#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

regionSplit <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat

  aln_profile <- readDNAStringSet(op_args$profile)
  char_aln_profile <- as.character(aln_profile)
  attr(char_aln_profile, "names") <- NULL

  region_map <- readLines(op_args$region_map)

  dir.create(file.path(op_dir, 'mapped_reads'), showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp('/home/phillipl/projects/MotifBinner2/code/MotifBinner2/src/trimFront.cpp')

  i <- 1
  registerDoMC(cores = config$ncpu)
  split_regions <- foreach (i = 1:length(seq_dat), .combine = bind_rows) %dopar% {
    cur_seq <- seq_dat@sread[i]
    names(cur_seq) <- seq_dat@id[i]
    
    seq_file_name <- file.path(op_dir, 'mapped_reads', paste('read_', sprintf("%06d",i), '.fasta', sep = ''))
    writeXStringSet(cur_seq,
                    seq_file_name,
                    width=20000)
    aligned_file_name <- file.path(op_dir, 'mapped_reads', paste('read_', sprintf("%06d",i), '_mapped.fasta', sep = ''))
    system(paste('mafft --quiet --addfragments ', seq_file_name, ' ', op_args$profile,
                 ' > ', aligned_file_name, sep = ''))
    stopifnot(file.exists(aligned_file_name))
    aligned_seqs <- readDNAStringSet(aligned_file_name)
    char_aligned_seqs <- as.character(aligned_seqs)
    attr(char_aligned_seqs, "names") <- NULL

    regionSplit_cpp

    regions <- regionSplit_cpp(char_aligned_seqs, char_aln_profile, region_map)$regions

    regions_df <- data.frame(to_delete = 1)
    for (j in names(regions)){
      regions_df <-
        cbind(regions_df, 
              data.frame(tmp = regions[j],
                         stringsAsFactors = F))
      names(regions_df)[ncol(regions_df)] <- intToUtf8(j)
    }
    regions_df$to_delete <- NULL
    regions_df
  }









  
  per_read_metrics <- data.frame(seq_len = width(seq_dat@sread))

  threshold <- op_args$threshold
  if (is.null(threshold))
  {
    threshold <- 295
  }
  trim_steps <- list(step1 = list(name = 'seq_len',
                                  threshold = threshold,
                                  comparator = `>=`,
                                  breaks = c(Inf, 300, 297, 295, 290, 280, -Inf)
                                  )
                    )

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- getKept(result, seq_dat=seq_dat)
  class(result) <- 'regionSplit'
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

saveToDisk.regionSplit <- function(result, config, seq_dat)
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

computeMetrics.regionSplit <- function(result, config, seq_dat)
{
  return(result)
}

print.regionSplit <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: regionSplit')
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
