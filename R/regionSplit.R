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

  registerDoMC(cores = config$ncpu)
  split_regions <- foreach (i = 1:length(seq_dat), .combine = bind_rows) %dopar%
#  timing <- list()
#  split_regions <- list()
#  for (i in 1:100)#length(seq_dat))
  {

#    ptm <- proc.time()
    
    cur_seq <- seq_dat@sread[i]
    names(cur_seq) <- seq_dat@id[i]
    
    seq_file_name <- file.path(op_dir, 'mapped_reads', paste('read_', sprintf("%06d",i), '.fasta', sep = ''))
    writeXStringSet(cur_seq,
                    seq_file_name,
                    width=20000)
    aligned_file_name <- file.path(op_dir, 'mapped_reads', paste('read_', sprintf("%06d",i), '_mapped.fasta', sep = ''))
#    timing$alignment_prep <- ifelse(is.null(timing$alignment_prep), proc.time() - ptm,
#                                    timing$alignment_prep + proc.time() - ptm)
#    ptm <- proc.time()
    system(paste('mafft --quiet --addfragments ', seq_file_name, ' ', op_args$profile,
                 ' > ', aligned_file_name, sep = ''))
#    timing$alignment <- ifelse(is.null(timing$alignment), proc.time() - ptm,
#                                    timing$alignment + proc.time() - ptm)
#    ptm <- proc.time()
    stopifnot(file.exists(aligned_file_name))
    aligned_seqs <- readDNAStringSet(aligned_file_name)
    char_aligned_seqs <- as.character(aligned_seqs)
    attr(char_aligned_seqs, "names") <- NULL

    gapped_qual <- transfer_gaps_cpp(char_aligned_seqs[length(char_aligned_seqs)],
                                     as.character(seq_dat@quality@quality[i]),
                                     -1)
    gapped_qual <- gapped_qual$quals

#    timing$alignment_proc <- ifelse(is.null(timing$alignment_proc), proc.time() - ptm,
#                                    timing$alignment_proc + proc.time() - ptm)
#    ptm <- proc.time()
    tmp_split <- regionSplit_cpp(char_aligned_seqs, char_aln_profile, region_map, gapped_qual)
#    timing$splitting <- ifelse(is.null(timing$splitting), proc.time() - ptm,
#                                    timing$splitting + proc.time() - ptm)
#    ptm <- proc.time()

    regions_df <- data.frame(to_delete = 1)
    regions_qual_df <- data.frame(to_delete = 1)
    for (j in names(tmp_split$regions)){
      regions_df <-
        cbind(regions_df, 
              data.frame(tmp = tmp_split$regions[j],
                         stringsAsFactors = F))
      names(regions_df)[ncol(regions_df)] <- intToUtf8(j)
      regions_qual_df <-
        cbind(regions_qual_df, 
              data.frame(tmp = tmp_split$regions_qual[j],
                         stringsAsFactors = F))
      names(regions_qual_df)[ncol(regions_qual_df)] <- paste("qual_", intToUtf8(j), sep = '')
    }
    regions_df$to_delete <- NULL
    regions_qual_df$to_delete <- NULL
#    timing$proc_splits <- ifelse(is.null(timing$proc_splits), proc.time() - ptm,
#                                    timing$proc_splits + proc.time() - ptm)
    #split_regions[[i]] <- cbind(regions_df, regions_qual_df)
    cbind(regions_df, regions_qual_df)
  }

  region_names <- names(split_regions)[!grepl("^qual_", names(split_regions))]
  res_seq_dat <- list()
  for (region_name in region_names){
    res_seq_dat[[region_name]] <- ShortReadQ(sread = DNAStringSet(split_regions[,region_name]),
                                             id = BStringSet(paste(as.character(seq_dat@id), 
                                                                   region_name, sep = '_')),
                                             quality = BStringSet(split_regions[,paste('qual_', 
                                                                                       region_name, sep = '')])
                                             )
  }

  per_read_metrics <- data.frame('read_exists' = rep(1, length(seq_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'regionSplit'
  if (op_args$cache){
    result$seq_dat <- res_seq_dat
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
  for (i in names(result$seq_dat))
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_kept_', result$config$op_args$name, '_region_', i, '.fastq', sep = ''))
    writeFastq(result$seq_dat[[i]], tmp_name, compress=F)
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
