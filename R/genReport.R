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
  new_report_name_md   <- paste(config$base_for_names, 'report.md', sep = '_')
  new_report_name_html <- paste(config$base_for_names, 'report.html', sep = '_')
  new_report_name_pdf  <- paste(config$base_for_names, 'report.pdf', sep = '_')


  dir.create(file.path(config$output_dir, config$base_for_names, 'full_report'),
             showWarnings = FALSE,
             recursive = TRUE)
  op_name <- names(all_results)[2]
  for (op_name in names(all_results))
  {
    report_file_name <- file.path(all_results[[op_name]]$config$op_dir,
                                  paste(config$base_for_names, all_results[[op_name]]$config$op_full_name, 'report.md', sep = '_'))
    stopifnot(file.exists(report_file_name))
    file.copy(report_file_name,
              file.path(config$output_dir, 
                        config$base_for_names, 
                        'full_report',
                        basename(report_file_name)),
              overwrite = TRUE)
    op_folders <- list.dirs(all_results[[op_name]]$config$op_dir)
    fig_folders <- op_folders[grep('/figure[^/]*$', op_folders)]
    if (length(fig_folders) > 0){
      command <- paste('cp -r', 
                       fig_folders, 
                       file.path(config$output_dir, 
                                 config$base_for_names, 
                                 'full_report'),
                       sep = ' ')
      system(command)
    }
  }
  all_report_files <-
  list.files(file.path(config$output_dir, 
                       config$base_for_names, 
                       'full_report'),
             '.md',
             full.names = FALSE)

  rmarkdown_templates_folder <- file.path(find.package('rmarkdown'),
                                          'rmd',
                                          'h')

  cwd <- getwd()
  setwd(file.path(config$output_dir, config$base_for_names, 'full_report'))

  html_header_extra <- c(
paste('<script src="', rmarkdown_templates_folder, '/jquery-1.11.3/jquery.min.js"></script>', sep = ''),
paste('<meta name="viewport" content="width=device-width, initial-scale=1" />', sep = ''),
paste('<link href="', rmarkdown_templates_folder,  '/bootstrap-3.3.5/css/bootstrap.min.css" rel="stylesheet" />', sep = ''),
paste('<script src="', rmarkdown_templates_folder, '/bootstrap-3.3.5/js/bootstrap.min.js"></script>', sep = ''),
paste('<script src="', rmarkdown_templates_folder, '/bootstrap-3.3.5/shim/html5shiv.min.js"></script>', sep = ''),
paste('<script src="', rmarkdown_templates_folder, '/bootstrap-3.3.5/shim/respond.min.js"></script>', sep = '')
)
  writeLines(html_header_extra, 'header_extra_stuff.html')

  pandoc_command <- paste(
    "/usr/bin/pandoc +RTS -K512m -RTS ",
    paste(all_report_files, collapse = ' '), " ",
    "--to html --mathjax --variable 'mathjax-url:https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML' ",
    "--from markdown+autolink_bare_uris+ascii_identifiers+tex_math_single_backslash ",
    "--output ", paste(config$base_for_names, '_full_report.html', sep = ''), " ",
#    "--output ", file.path(config$output_dir, config$base_for_names, paste(config$base_for_names, '_full_report.html', sep = '')),
    "--smart --email-obfuscation none --self-contained --standalone --section-divs ",
    "--toc ",
    "--variable 'theme:bootstrap' ",
    "--include-in-header header_extra_stuff.html ", 
    "--template ", rmarkdown_templates_folder, '/default.html ', 
    "--no-highlight --variable highlightjs=", rmarkdown_templates_folder, "/highlight ",
    "--variable navigationjs=", rmarkdown_templates_folder, "/navigation-1.1",
    sep = "")

  system(pandoc_command)
  file.copy(paste(config$base_for_names, '_full_report.html', sep = ''),
            paste('../', config$base_for_names, '_full_report.html', sep = ''))
  setwd(cwd)
  return(all_results)
}

genReport_internal_generic <- function(result, config)
{
  new_report_name_md   <- paste(config$base_for_names, result$config$op_full_name, 'report.md', sep = '_')
  new_report_name_html <- paste(config$base_for_names, result$config$op_full_name, 'report.html', sep = '_')
  new_report_name_pdf  <- paste(config$base_for_names, result$config$op_full_name, 'report.pdf', sep = '_')

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
           output_file = new_report_name_html, clean=FALSE, quiet=TRUE)
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
