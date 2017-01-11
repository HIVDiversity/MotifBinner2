#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

vsearchCluster <- function(all_results, config)
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

  id <- op_args$id
  min_clus_size <- op_args$min_clus_size
  if (is.null(min_clus_size)){
    min_clus_size <- 1
  }
  
  seq_dat <- all_results[[data_source_indx]]$seq_dat

  seq_file_name <- file.path(op_dir, 'seq.fastq')
  writeFastq(seq_dat,
             seq_file_name,
             compress=FALSE)

  lookup_file_name <- file.path(op_dir, 'lookup.tab')
  system(paste('vsearch', 
               ' --cluster_size ', seq_file_name, 
               ' --id ', id,
               ' --sizeout ',
               ' --uc ', lookup_file_name,
               sep = ''))
  
  assignments <- read.delim(lookup_file_name, header = F,
                            stringsAsFactors = F)
  cluster_sizes <- subset(assignments, 
                          V1 == 'C',
                          select = c("V1", "V2", "V3", "V9"))
  names(cluster_sizes) <- c('type', 'clus_id', 'size', 'centroid')
  centroid_names <- cluster_sizes[,c('centroid', 'size')]
  centroid_names$new_name <- paste(centroid_names$centroid, centroid_names$size, sep = '_')
  centroid_names$size <- NULL
  assignments <- assignments[,c(1,2,9,10),]
  
  names(assignments) <- c('type', 'clus_num', 'seq_name', 'assigned_to')
  assignments <- merge(assignments, centroid_names, 
                       by.x = 'seq_name', by.y = 'centroid',
                       all.x = TRUE)
  print(str(assignments))
  assignments$new_name <- ifelse(is.na(assignments$new_name), assignments$seq_name, assignments$new_name)
  per_read_metrics <- subset(assignments, type != 'C')

  seq_dat <- seq_dat[order(as.character(seq_dat@id))]
  per_read_metrics <- per_read_metrics[order(per_read_metrics$seq_name),]
  stopifnot(all(as.character(seq_dat@id) == per_read_metrics$seq_name))
  seq_dat@id <- BStringSet(per_read_metrics$new_name)
  per_read_metrics$seq_name <- per_read_metrics$new_name
  per_read_metrics$new_name <- NULL
  per_read_metrics$size <- as.numeric(gsub('^.*_', '', per_read_metrics$seq_name))

  per_read_metrics$assigned_to[per_read_metrics$type == 'S'] <- 'Centroid'
  per_read_metrics$is_centroid <- ifelse(per_read_metrics$type == 'S', 1, 0)
#  write.csv(assignments, gsub('.tab','.csv', lookup_file_name),
#            row.names=F)

  trim_steps <- list(step1 = list(name = 'is_centroid',
                                  threshold = 1,
                                  comparator = `>=`,
                                  breaks = c(Inf, 1, -Inf)),
                     step2 = list(name = 'size',
                                  threshold = min_clus_size,
                                  comparator = `>=`,
                                  breaks = rev(sort(unique(c(0, 1, 100, 1000, min_clus_size, min_clus_size-1, min_clus_size+1))))))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'vsearchCluster'
  result$seq_dat <- getKept(result, seq_dat)
  result$input_dat <- seq_dat
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  file.remove(lookup_file_name)
  file.remove(seq_file_name)
  return(result)
}

saveToDisk.vsearchCluster <- function(result, config, seq_dat = NULL)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')
  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)

  seq_file_name <- paste(config$base_for_names, '_kept_', result$config$op_args$name, '.fastq', sep = '')
  writeFastq(result$seq_dat, file.path(op_dir, seq_file_name), compress = F)
  clustering_file <- paste(config$base_for_names, '_clustering_', result$config$op_args$name, '.csv', sep = '')
  write.csv(result$metrics$per_read_metrics, file.path(op_dir, clustering_file), row.names=F)
  return(result)
}

computeMetrics.vsearchCluster <- function(result, config, seq_dat)
{
  return(result)
}

print.vsearchCluster <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: vsearchCluster')
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
