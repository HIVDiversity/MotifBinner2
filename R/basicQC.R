# Adapted from FastQC:
# Andrews S. (2010). FastQC: a quality control tool for high throughput
# sequence data. Available online at:
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc


#' Performs basic quality assessment
#' @inheritParams applyOperation
#' @export

basicQC <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  seq_dat <- all_results[[data_source_indx]]$seq_dat
  
  per_read_metrics <- data.frame('read_exists' = rep(1, length(seq_dat)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'basicQC'
  if (op_args$cache){
    stop('Do not cache data in steps that do not alter the datasets')
  } else {
    result$tmp <- seq_dat
  }
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.basicQC <- function(result, config)
{
  return(result)
}

#genSummary.basicQC <- function(result, config)
#{
#  summary_tab <- rbind(
#    genSummary_internal(operation = 'basicQC',
#                        parameters = 'fwd_reads',
#                        kept_seq_dat = result$kept$fwd_reads,
#                        trimmed_seq_dat = DNAStringSet(NULL)),
#    genSummary_internal(operation = 'basicQC',
#                        parameters = 'rev_reads',
#                        kept_seq_dat = result$kept$rev_reads,
#                        trimmed_seq_dat = DNAStringSet(NULL)))
#  result$summary <- summary_tab
#  write.csv(summary_tab, file.path(result$op_dir, 'basicQC_summary.csv'), row.names=FALSE)
#  return(result)
#}

computeMetrics.basicQC <- function(result, config, seq_dat)
{
  fastq_name_headers <- c("instrument", "run_num", "flowcell", "lane", "tile",
                          "xpos", "ypos", "read", "is_filtered", "control_num",
                          "index_seq")
  fastq_numeric <- c("run_num", "lane", "tile", "xpos", "ypos", "read",
                     "control_num", "index_seq")

  qual_mat <- as(FastqQuality(quality(quality(seq_dat))), 'matrix')
  per_read_quality <- apply(qual_mat, 1, sum, na.rm=T)
  pos_qual <- melt(qual_mat)
  names(pos_qual) <- c('read_num', 'cycle', 'qual')
  rm(qual_mat)
  
  seq_df <- data.frame(read_num = 1:length(seq_dat),
                       seq_name = as.character(seq_dat@id),
                       read_widths = width(seq_dat),
                       qual = per_read_quality / width(seq_dat),              
                       stringsAsFactors = F)
  rm(per_read_quality)
  
  seq_df <- seq_df %>%
    mutate(seq_name = gsub(' ', ':', seq_name)) %>%
    separate(seq_name, fastq_name_headers, ":") %>%
    mutate_each_(funs(as.numeric), fastq_numeric) %>%
    mutate(xcat = round(xpos/max(xpos), 2)) %>%
    mutate(ycat = round(ypos/max(ypos), 2))
  
  zone_qual <- seq_df %>%
    group_by(xcat, ycat) %>%
    summarize(zone_quality = mean(qual, na.rm=T),
              n_seq = sum(!is.na(qual)))

  tile_qual <- pos_qual %>%
    inner_join(seq_df %>% select(read_num, tile), by = 'read_num') %>%
    mutate(tile_indx = dense_rank(tile)) %>%
    select(tile_indx, cycle, qual) %>%
    group_by(tile_indx, cycle) %>%
    summarize(avg_qual = mean(qual, na.rm = T))

  pos_qual <- pos_qual %>%
    mutate(cycle_cat = trunc((cycle+1)/2)*2) %>%
    select(cycle_cat, qual) %>%
    group_by(cycle_cat, qual) %>%
    summarize(count = n())

  result$metrics <- list(
    seq_df = seq_df,
    zone_qual = zone_qual,
    tile_qual = tile_qual,
    pos_qual = pos_qual
  )
  
  return(result)
}

print.basicQC <- function(result, config)
{
  cat('\n------------------')
  cat('\nOperation: basicQC')
  cat('\n------------------')
  cat('\nLoaded Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs')])
  invisible(result)
}

