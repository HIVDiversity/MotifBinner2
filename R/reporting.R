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
  if (class(result) == 'all_results'){
    setwd(file.path(config$output_dir))
  } else {
    setwd(result$op_dir)
  }
  knit(report_file_name, new_report_name_md)

  if ('html' %in% config$report_type){
    render(new_report_name_md, output_format = 'html_document',
           output_file = new_report_name_html, clean=FALSE)
  }
  if ('pdf' %in% config$report_type){
    render(new_report_name_md, output_format = 'pdf_document',
           output_file = new_report_name_pdf, clean=FALSE)
  }
  setwd(cwd)
  return(result)
}

