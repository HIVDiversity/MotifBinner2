#' Tallies to true sizes of clusters of clusters
#' @inheritParams applyOperation
#' @export

doubleClusterTally <- function(all_results, config)
{

  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)

  s1_strip_pattern <- op_args$s1_strip_pattern # "_[0-9]*$"
  s2_strip_pattern <- op_args$s2_strip_pattern # "_[0-9]*_*[xyz]*_*[0-9]*$"
  
  s1_data_source <- op_args$data_source['s1']
  s2_data_source <- op_args$data_source['s2']

  data_source_indx <- grep(s1_data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)

  step1 <- all_results[[data_source_indx]]$metrics$per_read_metrics
  names(step1) <- paste("s1", names(step1), sep = "_")
  step1_sizes <- step1[step1$s1_assigned_to == "Centroid", c("s1_clus_num", "s1_size")]
  names(step1_sizes)[2] <- "s1_assigned_to_size"
  step1 <- merge(step1, step1_sizes)
  
  step1$s1_cent <- ifelse(step1$s1_assigned_to == "Centroid", 
                            gsub(s1_strip_pattern, "", step1$s1_seq_name), 
                            gsub(s1_strip_pattern, "", step1$s1_assigned_to))

  data_source_indx <- grep(s2_data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)

  step2 <- all_results[[data_source_indx]]$metrics$per_read_metrics
  s2_seq_dat <- all_results[[data_source_indx]]$seq_dat
  names(step2) <- paste("s2", names(step2), sep = "_")
  step2$s1_like_name <- gsub(s2_strip_pattern, "", step2$s2_seq_name)
  step2_sizes <- step2[step2$s2_assigned_to == "Centroid", c("s2_clus_num", "s2_size")]
  names(step2_sizes)[2] <- "s2_assigned_to_size"
  step2 <- merge(step2, step2_sizes)
  step2$s2_cent <- ifelse(step2$s2_assigned_to == "Centroid", 
                            gsub(s2_strip_pattern, "", step2$s2_seq_name), 
                            gsub(s2_strip_pattern, "", step2$s2_assigned_to))

  all_clusts <- merge(step1, step2,
                      by.x = 's1_cent', by.y='s1_like_name',
                      all.x = T, all.y = T)
  stopifnot(nrow(all_clusts) == nrow(step1))
  
  #table(table(all_clusts$s1_seq_name, useNA='always'), useNA='always')
  #table(step2$s2_is_centroid)
  
  true_cluster_sizes <-
  all_clusts %>%
    group_by(s2_clus_num) %>%
    summarize(n_s1_sequences = n(),
              cent_name = unique(s2_cent),
#              s2_seq_name_orig = unique(s2_seq_name),
              n_s2_clusts = unique(s2_assigned_to_size))

  true_cluster_sizes <- subset(true_cluster_sizes, !is.na(cent_name))

  true_cluster_sizes$new_name <-
    paste(true_cluster_sizes$cent_name, 
          '_', true_cluster_sizes$n_s2_clusts,
          '_', true_cluster_sizes$n_s1_sequences,
          sep = '')

  stopifnot(all(gsub(s2_strip_pattern, '', as.character(s2_seq_dat@id)) %in% 
                  true_cluster_sizes$cent_name))

  new_names_in_order <-
    true_cluster_sizes$new_name[match(gsub(s2_strip_pattern, '', as.character(s2_seq_dat@id)), 
                                      true_cluster_sizes$cent_name)]
  stopifnot(all(gsub('_[0-9]*_[0-9]*$', '', new_names_in_order) == 
                gsub(s2_strip_pattern, '', as.character(s2_seq_dat@id))))

  s2_seq_dat@id <- BStringSet(new_names_in_order)

  seq_dat <- s2_seq_dat
  
  per_read_metrics <- true_cluster_sizes
  per_read_metrics$read_exists <- 1

  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  kept_dat <- getKept(result, seq_dat=seq_dat)
  class(result) <- 'doubleClusterTally'
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

saveToDisk.doubleClusterTally <- function(result, config, seq_dat)
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

computeMetrics.doubleClusterTally <- function(result, config, seq_dat)
{
  return(result)
}

print.doubleClusterTally <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: doubleClusterTally')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}
