# Adapted from FastQC:
# Andrews S. (2010). FastQC: a quality control tool for high throughput
# sequence data. Available online at:
# http://www.bioinformatics.babraham.ac.uk/projects/fastqc


#' Performs basic quality assessment
#' @inheritParams applyOperation
#' @export

basicQC <- function(all_results, config)
{
  opdir <- file.path(config$output_dir, config$prefix_for_names,
                      paste('n', sprintf("%03d", length(all_results)+1), 'basicQC', sep = ''))
  dir.create(opdir, showWarnings = FALSE, recursive = TRUE)
  
  result <- list(final = all_results[[length(all_results)]]$final,
                 step_num = length(all_results)+1,
                 opdir = opdir)
  
  class(result) <- 'basicQC'
  return(result)
}

