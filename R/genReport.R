#' Generates a report for a operation or the full process
#' @inheritParams saveToDisk
#' @export

genReport <- function(result, config)
{
  if (class(result) == 'allResults')
  { 
    result <- genReport_internal_allResults(result, config)
  } else {
    result <- genReport_internal_generic(result, config)
  }
  
  return(result)
}

genReport_internal_allResults <- function(all_results, config)
{
  new_report_name_md   <- paste(config$base_for_names, 'bin_report.md', sep = '_')
  new_report_name_html <- paste(config$base_for_names, 'bin_report.html', sep = '_')
  new_report_name_pdf  <- paste(config$base_for_names, 'bin_report.pdf', sep = '_')

  report_file_name <- file.path(find.package('MotifBinner2'),
                                'reports',
                                paste(class(all_results), '.Rmd', sep = ''))
  if (!file.exists(report_file_name))
  {
    report_file_name <- file.path(find.package('MotifBinner2'),
                                  'inst', 'reports',
                                  paste(class(all_results), '.Rmd', sep = ''))
  }

  cwd <- getwd()
  setwd(file.path(config$output_dir, config$base_for_names))
  knit(report_file_name, new_report_name_md, quiet = TRUE)

  if ('html' %in% config$report_type){
    render(new_report_name_md, output_format = 'html_document',
           output_file = new_report_name_html, clean=TRUE, quiet=TRUE)
  }
  if ('pdf' %in% config$report_type){
    warning('pdf reports not implemented yet')
    # pandoc md -> tex; manual hacking in .tex; pandoc tex -> pdf
#    render(report_file_name, output_format = 'pdf_document',
#           output_file = new_report_name_pdf, clean=TRUE, quiet=TRUE)
  }
  setwd(cwd)
  return(all_results)
}

genReport_internal_generic <- function(result, config)
{
  new_report_name_md   <- paste(config$base_for_names, class(result), 'bin_report.md', sep = '_')
  new_report_name_html <- paste(config$base_for_names, class(result), 'bin_report.html', sep = '_')
  new_report_name_pdf  <- paste(config$base_for_names, class(result), 'bin_report.pdf', sep = '_')

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
  setwd(result$config$op_dir)
  knit(report_file_name, new_report_name_md, quiet = TRUE)

  if ('html' %in% config$report_type){
    render(new_report_name_md, output_format = 'html_document',
           output_file = new_report_name_html, clean=TRUE, quiet=TRUE)
  }
  if ('pdf' %in% config$report_type){
    warning('pdf reports not implemented yet')
    # pandoc md -> tex; manual hacking in .tex; pandoc tex -> pdf
    #render(new_report_name_md, output_format = 'pdf_document',
    #       output_file = new_report_name_pdf, clean=TRUE, quiet=TRUE)
  }
  setwd(cwd)
  return(result)
}
