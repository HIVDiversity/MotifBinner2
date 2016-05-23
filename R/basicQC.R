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
  return(result)
}

print.basicQC <- function(result, config)
{
  cat('\nbasicQC\n')
  print(result$summary)
  return(result)
}

