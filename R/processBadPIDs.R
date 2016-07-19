#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

processBadPIDs <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_source_indx <- grep(op_args$data_source, names(all_results))
  stopifnot(length(data_source_indx) == 1)
  stopifnot(all(sort(names(all_results[[data_source_indx]]$seq_dat)) == c("fwd", "rev")))
  seq_dat_fwd <- all_results[[data_source_indx]]$seq_dat$fwd
  seq_dat_rev <- all_results[[data_source_indx]]$seq_dat$rev
  stopifnot(all(seq_dat_fwd@id == seq_dat_rev@id))

  pid_len <- nchar(bin_metrics$pid[1])

  per_read_metrics <- data.frame(read_name = as.character(seq_dat_fwd@id),
                                 stringsAsFactors = F)
  per_read_metrics$pid <- gsub("^.*_PID:" , "", per_read_metrics$read_name)
  per_read_metrics$clean_pid <- gsub("_" , "", per_read_metrics$pid)

### group on pids and compare to consensus cutoff
  bin_metrics <- data.frame(pid = names(table(per_read_metrics$clean_pid)),
                            raw_size = as.numeric(table(per_read_metrics$clean_pid)),
                            stringsAsFactors = F)
  cc <- get_consensus_cutoff(pid_len, max(bin_metrics$raw_size),
                             get_cutoff_models(), phi = 1/100)
  bin_metrics$big_enough <- bin_metrics$raw_size > 2*cc

### find potential parents
  bin_metrics$dist_to_big <- -1
  big_enough_pids <- bin_metrics$pid[bin_metrics$big_enough]
  bin_metrics$parent_pid <- ""
  bin_metrics$parent_size <- -1
  bin_metrics$parentage_conflict <- -1
  for (i in 1:nrow(bin_metrics)){
    if (!bin_metrics$big_enough[i]){
      dist_to_bigs <- stringdist::stringdist(bin_metrics$pid[i], big_enough_pids, 'hamming')
      bin_metrics$dist_to_big[i] <- min(dist_to_bigs)
      if (min(dist_to_bigs) <= 2){
        can_pids <- big_enough_pids[which(dist_to_bigs == min(dist_to_bigs))]
        can_pids_size <- bin_metrics$raw_size[match(can_pids, bin_metrics$pid)]
        biggest_can_id <- which(can_pids_size == max(can_pids_size))
        bin_metrics[i,'parent_pid'] <- can_pids[biggest_can_id[1]]
        bin_metrics[i,'parent_size'] <- can_pids_size[biggest_can_id[1]]
        if (length(biggest_can_id)>1){
          bin_metrics$parentage_conflict[i] <- length(biggest_can_id)
        }
      }
    }
  }

  seq_err <- 1/100
  no_seq_err <- 99/100
  min_bin_size <- 3
  max_chimeric_content <- 0.3
  offspring_prob_cutoff <- sum(bin_metrics$big_enough)/(4^pid_len)


  bin_metrics$off_spring_prob <- NA_real_
### compute prob that content from parent
  for (i in 1:nrow(bin_metrics)){
    if (!bin_metrics$big_enough[i]){
      if (bin_metrics$parent_pid[i] == ""){
        off_spring_prob <- 0
      } else {
        success_prob <- ((seq_err/3)^bin_metrics$dist_to_big[i])*(no_seq_err^(pid_len - bin_metrics$dist_to_big[i]))
        off_spring_prob <- 1-pbinom(trunc(bin_metrics$raw_size[i]*max_chimeric_content), round(bin_metrics$parent_size[i]/no_seq_err,0), success_prob)
      }
      bin_metrics$off_spring_prob[i] <- off_spring_prob 
    }
  }
  bin_metrics$is_offspring <- bin_metrics$off_spring_prob > offspring_prob_cutoff
  bin_metrics$true_bin <- bin_metrics$big_enough
  for (i in 1:nrow(bin_metrics)){
    if (!bin_metrics$big_enough[i]){
      if (!bin_metrics$is_offspring[i] & bin_metrics$raw_size[i] > min_bin_size){
        bin_metrics$true_bin[i] <- TRUE
      }
    }
  }
  bin_metrics[is.na(bin_metrics$off_spring_prob),'off_spring_prob'] <- 0

  per_read_metrics <- merge(per_read_metrics, bin_metrics[,c('pid', 'off_spring_prob', 'true_bin', 'raw_size')], 
                            by.x = 'clean_pid', by.y = 'pid', all.x=T)
  stopifnot(sum(is.na(per_read_metrics$off_spring_prob)) == 0)
  per_read_metrics <- per_read_metrics[match(as.character(seq_dat_fwd@id), per_read_metrics$read_name),]
  stopifnot(all(per_read_metrics$read_name == as.character(seq_dat_fwd@id)))

  trim_steps <- list(step1 = list(name = 'off_spring_prob',
                                  threshold = offspring_prob_cutoff,
                                  comparator = `<=`,
                                  breaks = c(-Inf, offspring_prob_cutoff*c(0.01, 0.1, 1, 10, 100), Inf)
                                  ),
                     step2 = list(name = 'raw_size',
                                  threshold = min_bin_size,
                                  comparator = `>`,
                                  breaks = c(Inf, min_bin_size + 3:-3, -Inf)
                                  )
                    )

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  
  kept_dat <- list('fwd' = seq_dat_fwd[per_read_metrics$true_bin],
                   'rev' = seq_dat_rev[per_read_metrics$true_bin])
  trim_dat <- list('fwd' = seq_dat_fwd[!per_read_metrics$true_bin],
                   'rev' = seq_dat_rev[!per_read_metrics$true_bin])
  
  class(result) <- 'processBadPIDs'
  if (!op_args$cache){stop('must cache on processBadPIDs')}

  result$seq_dat <- kept_dat
  result$trim_dat <- trim_dat
  result$input_dat <- list(fwd = seq_dat_fwd,
                           rev = seq_dat_rev)

  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

saveToDisk.processBadPIDs <- function(result, config, seq_dat)
{
  kept <- result$seq_dat
  trimmed <- result$trim_dat

  if (length(kept[['fwd']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_fwd_kept_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(kept[['fwd']], tmp_name, compress=F)
  }
  if (length(kept[['rev']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_rev_kept_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(kept[['rev']], tmp_name, compress=F)
  }
  if (length(trimmed[['fwd']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_fwd_trimmed_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(trimmed[['fwd']], tmp_name, compress=F)
  }
  if (length(trimmed[['rev']]) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_rev_trimmed_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(trimmed[['rev']], tmp_name, compress=F)
  }
  return(result)
}

genSummary_matchPairs <- function(result)
{

}


computeMetrics.processBadPIDs <- function(result, config, seq_dat)
{
  return(result)
}

print.processBadPIDs <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: processBadPIDs')
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
