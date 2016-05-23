# Adapted from FastQC:
# Andrews S. (2010). FastQC: a quality control tool for high throughput
# sequence data. Available online at:
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc


#' Performs basic quality assessment
#' @inheritParams applyOperation
#' @export

basicQC <- function(all_results, config)
{
  op_dir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), '_basicQC', sep = ''))
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  
  result <- list(final = all_results[[length(all_results)]]$final,
                 step_num = length(all_results)+1,
                 op_dir = op_dir)
  
  class(result) <- 'basicQC'
  return(result)
}

saveToDisk.basicQC <- function(result, config)
{
  return(result)
}

genSummary.basicQC <- function(result, config)
{
  summary_tab <- rbind(
    genSummary_internal(operation = 'basicQC',
                        parameters = 'fwd_reads',
                        kept_seq_dat = result$final$fwd_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)),
    genSummary_internal(operation = 'basicQC',
                        parameters = 'rev_reads',
                        kept_seq_dat = result$final$rev_reads,
                        trimmed_seq_dat = DNAStringSet(NULL)))
  result$summary <- summary_tab
  write.csv(summary_tab, file.path(result$op_dir, 'basicQC_summary.csv'), row.names=FALSE)
  return(result)
}

computeMetrics.basicQC <- function(result, config)
{
  result$metrics <- list()
  for (data_set_name in names(result$final)){
    seq_dat <- result$final[[data_set_name]]
    
    fastq_name_headers <- c("instrument", "run_num", "flowcell", "lane", "tile",
                            "xpos", "ypos", "read", "is_filtered", "control_num",
                            "index_seq")
    fastq_numeric <- c("run_num", "lane", "tile", "xpos", "ypos", "read",
                       "control_num", "index_seq")
    
    per_read_quality <- letterFrequency(seq_dat@quality@quality, 
                                        names(encoding(seq_dat@quality)))
    per_read_quality <- (per_read_quality %*% matrix(encoding(seq_dat@quality), ncol=1))[,1]
    
    seq_df <- data.frame(seq_name = as.character(seq_dat@id),
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

    result$metrics[[data_set_name]] <- list(
      seq_df = seq_df,
      zone_qual = zone_qual
    )
  }
  return(result)
}

print.basicQC <- function(result, config)
{
  cat('\n------------------')
  cat('\nOperation: basicQC')
  cat('\n------------------')
  cat('\nLoaded Sequences:\n')
  print(result$summary[,c('parameters', 'seqs_kept')])
  return(result)
}

