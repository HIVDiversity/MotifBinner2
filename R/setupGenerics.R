#' Save the relevant results of an operation to disk
#' @param result The result produced from a specific operation
#' @inheritParams applyOperation
#' @export

saveToDisk <- function(result, config)
{
  UseMethod('saveToDisk')
}

#' Generates a report for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genReport <- function(result, config)
{
  new_report_name_md <- paste(config$prefix_for_names, 'bin_report.md', sep = '_')                         
  new_report_name_html <- paste(config$prefix_for_names, 'bin_report.html', sep = '_')                     
  new_report_name_pdf <- paste(config$prefix_for_names, 'bin_report.pdf', sep = '_')                       

  report_file_name <- file.path(find.package('MotifBinner2'),
                                'reports',
                                paste(class(result), '.Rmd', sep = ''))
  if (!file.exists(report_file_name))
  {
    report_file_name <- file.path(find.package('MotifBinner2'),
                                  'inst', 'reports',
                                  paste(class(result), '.Rmd', sep = ''))
  }

  cwd <- getwd()                                                                                    
  setwd(output)                                                                                     
  knit(knitr_file_location, new_report_name_md)                                                     

  if ('html' %in% config$report_type){
    render(report_file_name, output_format = 'html_document',
           output_file = new_report_name_html, clean=FALSE)
  }
  if ('pdf' %in% config$report_type){
    render(report_file_name, output_format = 'pdf_document',
           output_file = new_report_name_pdf, clean=FALSE)
  }
  setwd(cwd)
  return(result)
}

#' Generates a summary for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genSummary <- function(result, config)
{
  UseMethod('genSummary')
}


