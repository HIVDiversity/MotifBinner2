store_configs <- function()
{
#buildConfig <- function(fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
#                        rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
#                        fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
#                        profile_file, output_dir, base_for_names,
#                        intermediate_reports = TRUE,
#                        verbosity = 3,
#                        report_type = c('html')
#                        )
#  buildConfig(fwd_file = 
#              fwd_primer_seq = 
#              fwd_primer_lens = 
#              fwd_min_score = 
#              rev_file = 
#              rev_primer_seq = 
#              rev_primer_lens = 
#              rev_min_score = 
#              fwd_pid_in_which_fragment = 
#              rev_pid_in_which_fragment = 
#              profile_file = 
#              output_dir = 
#              base_for_names = 
#              )
  config <-
  buildConfig(overlap = FALSE,
              fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R1.fastq",
              fwd_primer_seq = 'CTCTTTTGACCCAATTCCTATACATTATTG',
              fwd_primer_lens = 30,
              fwd_min_score = 20,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R2.fastq",
              rev_primer_seq = 'NNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT',
              rev_primer_lens = c(14, 28),
              rev_min_score = 28,
              fwd_pid_in_which_fragment = "NULL",
              rev_pid_in_which_fragment = 1,
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP129_2040_009wpi_C2C3",
              erase_history = FALSE
              )

  config <-
  buildConfig(fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R1.fastq",
              fwd_primer_seq = 'CTCTTTTGACCCAATTCCTATACATTATTG',
              fwd_primer_lens = 30,
              fwd_min_score = 20,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R2.fastq",
              rev_primer_seq = 'NNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT',
              rev_primer_lens = c(14, 28),
              rev_min_score = 28,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "/fridge/data/MotifBinner2_test/c2c3_mol_clock_profile_1.fasta",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP129_2040_009wpi_C2C3",
              erase_history = FALSE
              )
  
  config <-
  buildConfig(fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq",
              fwd_primer_seq = 'TATGGGAYSAAAGYCTMAARCCATGTG',
              fwd_primer_lens = 27,
              fwd_min_score = 22,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R2.fastq",
              rev_primer_seq = 'CACACGCTCAGNNNNNNNNNATTCCATGTGTACATTGTACTGTRCTG',
              rev_primer_lens = c(11, 9, 27),
              rev_min_score = 40,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 2,
              profile_file = "/fridge/data/MotifBinner2_test/v1v2_profile1.fasta",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP256_3100_030wpi_v1v2"
              )
  config <-  # mol clock
  buildConfig(fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP129_3100_028wpi_C2C3_R1.fastq",
              fwd_primer_seq = 'NNNNNNNCTCTTTTGACCCAATTCCTATACATTATTG',
              fwd_primer_lens = 37,
              fwd_min_score = 32,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP129_3100_028wpi_C2C3_R2.fastq",
              rev_primer_seq = 'NNNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT',
              rev_primer_lens = c(15, 28),
              rev_min_score = 37,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "/fridge/data/MotifBinner2_test/v1v2_profile1.fasta",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP129_3100_028wpi_c2c3",
              erase_history = FALSE
              )
  config <-  # mol clock
  buildConfig(fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP008_2050_010wpi_C2C3_R1.fastq",
              fwd_primer_seq = 'NNNNNNNCTCTTTTGACCCAATTCCTATACATTATTG',
              fwd_primer_lens = 37,
              fwd_min_score = 32,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP008_2050_010wpi_C2C3_R2.fastq",
              rev_primer_seq = 'NNNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT',
              rev_primer_lens = c(15, 28),
              rev_min_score = 37,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "/fridge/data/MotifBinner2_test/v1v2_profile1.fasta",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP008_2050_010wpi_c2c3",
              erase_history = FALSE
              )
  config <-  # mol clock
  buildConfig(fwd_file = "/fridge/data/MotifBinner2_test/raw/CAP239_3090_022wpi_C2C3_R1.fastq",
              fwd_primer_seq = 'NNNNNNNCTCTTTTGACCCAATTCCTATACATTATTG',
              fwd_primer_lens = 37,
              fwd_min_score = 32,
              rev_file = "/fridge/data/MotifBinner2_test/raw/CAP239_3090_022wpi_C2C3_R2.fastq",
              rev_primer_seq = 'NNNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT',
              rev_primer_lens = c(15, 28),
              rev_min_score = 37,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "/fridge/data/MotifBinner2_test/v1v2_profile1.fasta",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "CAP239_3090_022wpi_c2c3",
              erase_history = FALSE
              )
  config <-    ## Zhou et al 2015 env region
  buildConfig(fwd_file = "/fridge/data/zhou2015/data/SRR1761912/SRR1761912_1.fastq",
              fwd_primer_seq = "NNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA",
              fwd_primer_lens = 34,
              fwd_min_score = 23,
              rev_file = "/fridge/data/zhou2015/data/SRR1761912/SRR1761912_2.fastq",
              rev_primer_seq = "NNNNNNNNNCAGTCCATTTTGCTCTACTAATGTTACAATGTGC",
              rev_primer_lens = c(9, 34),
              rev_min_score = 29,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "/fridge/data/zhou2015/NA_aligned.fasta",
              fwd_profile_file = "/fridge/data/zhou2015/NA_aligned_fwd.fasta",
              rev_profile_file = "/fridge/data/zhou2015/NA_aligned_rev.fasta",
              pattern_to_chop_from_names = "\\.[12] [0-9]+ length=[0-9].*",
              output_dir = "/fridge/data/zhou2015/binned",
              base_for_names = "zhou_2015_v1v3_1761912",
              erase_history = FALSE,
              ncpu = 2,
              bins_to_process = Inf
              )
  config <-    ## Zhou et al 2016 env region
  buildConfig(fwd_file = "/fridge/data/zhou2016/raw/SRR3221805_1.fastq",
              fwd_primer_seq = "NNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA",
              fwd_primer_lens = 34,
              fwd_min_score = 23,
              rev_file = "/fridge/data/zhou2016/raw/SRR3221805_2.fastq",
              rev_primer_seq = "NNNNNNNNNCAGTCCATTTTGCTCTACTAATGTTACAATGTGC",
              rev_primer_lens = c(9, 34),
              rev_min_score = 29,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              min_read_length = 245,
              profile_file = "/fridge/data/zhou2015/NA_aligned.fasta",
              fwd_profile_file = "/fridge/data/zhou2015/NA_aligned_fwd.fasta",
              rev_profile_file = "/fridge/data/zhou2015/NA_aligned_rev.fasta",
              pattern_to_chop_from_names = "\\.[12] [0-9]+ length=[0-9].*",
              output_dir = "/fridge/data/zhou2016/binned",
              base_for_names = "zhou_2016_v1v3_3221805",
              erase_history = FALSE,
              ncpu = 4,
              bins_to_process = Inf
              )
  config <-    ## Zhou et al 2016 env region
  buildConfig(fwd_file = "/fridge/data/zhou2016/raw/SRR3221806_1.fastq",
              fwd_primer_seq = "NNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA",
              fwd_primer_lens = 34,
              fwd_min_score = 23,
              rev_file = "/fridge/data/zhou2016/raw/SRR3221806_2.fastq",
              rev_primer_seq = "NNNNNNNNNCAGTCCATTTTGCTCTACTAATGTTACAATGTGC",
              rev_primer_lens = c(9, 34),
              rev_min_score = 29,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              min_read_length = 245,
              profile_file = "/fridge/data/zhou2015/NA_aligned.fasta",
              fwd_profile_file = "/fridge/data/zhou2016/profile/SRR3221806_1_profile.fasta",
              rev_profile_file = "/fridge/data/zhou2016/profile/SRR3221806_2_profile.fasta",
              pattern_to_chop_from_names = "\\.[12] [0-9]+ length=[0-9].*",
              output_dir = "/fridge/data/zhou2016/binned",
              base_for_names = "zhou_2016_v1v3_3221806",
              erase_history = FALSE,
              ncpu = 4,
              bins_to_process = Inf
              )
  config <-    ## Zhou et al 2016 env region
  buildConfig(fwd_file = "/fridge/data/zhou2016/raw/SRR3221818_1.fastq",
              fwd_primer_seq = "NNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA",
              fwd_primer_lens = 34,
              fwd_min_score = 23,
              rev_file = "/fridge/data/zhou2016/raw/SRR3221818_2.fastq",
              rev_primer_seq = "NNNNNNNNNCAGTCCATTTTGCTCTACTAATGTTACAATGTGC",
              rev_primer_lens = c(9, 34),
              rev_min_score = 29,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              min_read_length = 245,
              profile_file = "/fridge/data/zhou2015/NA_aligned.fasta",
              fwd_profile_file = "/fridge/data/zhou2016/profile/SRR3221818_1_profile.fasta",
              rev_profile_file = "/fridge/data/zhou2016/profile/SRR3221818_2_profile.fasta",
              pattern_to_chop_from_names = "\\.[12] [0-9]+ length=[0-9].*",
              output_dir = "/fridge/data/zhou2016/binned",
              base_for_names = "zhou_2016_v1v3_3221818",
              erase_history = FALSE,
              ncpu = 4,
              bins_to_process = Inf
              )

  config <-    ## Anthony 2016 env V3-C5 region
  buildConfig(fwd_file = '/fridge/data/MotifBinner2_test/raw/CAP008_1070_002wpi_V3B_C5A_R1.fastq',
              fwd_primer_seq = "NNNNNNNNGGACCAGGACAAACATTCTATGC",
              fwd_primer_lens = 31,
              fwd_min_score = 27,
              rev_file = '/fridge/data/MotifBinner2_test/raw/CAP008_1070_002wpi_V3B_C5A_R2.fastq',
              rev_primer_seq = "NNNNNNNNNNNNNNNGTCCYTCATATYTCCTCCTYCAGG",
              rev_primer_lens = c(15, 24),
              rev_min_score = 34,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              min_read_length = 295,
              profile_file = "/fridge/data/MotifBinner2_test/profile/CAP008_1070_002wpi_V3B_C5A_R1_profile.fasta",
              fwd_profile_file = "",
              rev_profile_file = "",
              output_dir = "/fridge/data/MotifBinner2_test/binned",
              base_for_names = "CAP008_1070_002wpi_V3B_C5A",
              erase_history = FALSE,
              ncpu = 6,
              bins_to_process = Inf
              )


  config <-    ## Anthony 2016 PCR recomb experiments
  buildConfig(fwd_file = "/fridge/data/colins_clone_pcr_recomb_data/raw/ddPCR_Pool_4_long_R1.fastq",
              fwd_primer_seq = "NNNNNNNGCTGGTTATGCGATTCTAAAGTG",
              fwd_primer_lens = 30,
              fwd_min_score = 25,
              fwd_pid_in_which_fragment = NULL,
              rev_file = "/fridge/data/colins_clone_pcr_recomb_data/raw/ddPCR_Pool_4_long_R2.fastq",
              rev_primer_seq = "NNNNNNNNNNNNNNNTGTGTTGTAAYTTCTAGRTC",
              rev_primer_lens = c(15,20),
              rev_min_score = 30,
              rev_pid_in_which_fragment = 1,
              min_read_length = 295,
              output_dir = "/fridge/data/colins_clone_pcr_recomb_data",
              base_for_names ="ddPCR_pool_4_long",
              ncpu = 6,
              fwd_profile_file = "",
              rev_profile_file = "",
              erase_history = FALSE,
              bins_to_process = Inf
             )

  config <-  # mol clock 2
  buildConfig(fwd_file = "/fridge/data/molecular_clock2/raw/HVTN503_159451817_1012_C2C3_R1.fastq",
              fwd_primer_seq = 'NNNNNNNCTCTTTTGACCCAATTCCTATACATTATTG',
              fwd_primer_lens = 37,
              fwd_min_score = 32,
              rev_file = "/fridge/data/molecular_clock2/raw/HVTN503_159451817_1012_C2C3_R2.fastq",
              rev_primer_seq = 'NNNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT',
              rev_primer_lens = c(15, 28),
              rev_min_score = 37,
              fwd_pid_in_which_fragment = NULL,
              rev_pid_in_which_fragment = 1,
              profile_file = "",
              output_dir = "/fridge/data/MotifBinner2_test",
              base_for_names = "HVTN503_159451817_1012_C2C3",
              erase_history = FALSE
              )
}

dummy_test_debug <- function()
{
  options(error=NULL)
  getwd()
  library(devtools)
  setwd('~/projects/MotifBinner2/code/MotifBinner2')
  load_all(quiet=TRUE)
  unlink(file.path(config$output_dir, config$base_for_names), recursive=T)
#  load_all()
  all_results <- list()
  class(all_results) <- 'allResults'

  all_results <- applyOperation(all_results, config, op_number = 'n001') # fwd_loadData
  all_results <- applyOperation(all_results, config, op_number = 'n002') # fwd_basicQC
  all_results <- applyOperation(all_results, config, op_number = 'n003') # fwd_ambigSeqs
  all_results <- applyOperation(all_results, config, op_number = 'n004') # fwd_primerDimer
  all_results <- applyOperation(all_results, config, op_number = 'n005') # fwd_seqLength
  all_results <- applyOperation(all_results, config, op_number = 'n006') # fwd_qualTrim
  all_results <- applyOperation(all_results, config, op_number = 'n007') # fwd_trimAffixes
  all_results <- applyOperation(all_results, config, op_number = 'n008') # rev_loadData
  all_results <- applyOperation(all_results, config, op_number = 'n009') # rev_basicQC
  all_results <- applyOperation(all_results, config, op_number = 'n010') # rev_ambigSeqs
  all_results <- applyOperation(all_results, config, op_number = 'n011') # rev_primerDimer
  all_results <- applyOperation(all_results, config, op_number = 'n012') # rev_seqLength
  all_results <- applyOperation(all_results, config, op_number = 'n013') # rev_qualTrim
  all_results <- applyOperation(all_results, config, op_number = 'n014') # rev_trimAffixes
  all_results <- applyOperation(all_results, config, op_number = 'n015') # fwdExtractPIDs
  all_results <- applyOperation(all_results, config, op_number = 'n016') # revExtractPIDs
  all_results <- applyOperation(all_results, config, op_number = 'n017') # matchPairs
  all_results <- applyOperation(all_results, config, op_number = 'n018') # processBadPIDs
  all_results <- applyOperation(all_results, config, op_number = 'n019') # 
  all_results <- applyOperation(all_results, config, op_number = 'n020') # 
  all_results <- applyOperation(all_results, config, op_number = 'n021') # 
  all_results <- applyOperation(all_results, config, op_number = 'n022') # 
  all_results <- applyOperation(all_results, config, op_number = 'n023') #
  all_results <- applyOperation(all_results, config, op_number = 'n024') #
  all_results <- applyOperation(all_results, config, op_number = 'n025') #
  all_results <- applyOperation(all_results, config, op_number = 'n026') #
  
  all_results <- applyOperation(all_results, config, op_number = 'n030') # primerSeqErr
  all_results <- applyOperation(all_results, config, op_number = 'n031') # binSeqErr
  
  all_results <- applyOperation(all_results, config, op_number = 'n100') # dataTracing

  genReport(all_results, config)
  op_number <- 'n031'
  config$current_op_number <- op_number

  result <- all_results

  all_results <- applyOperation(all_results, config, op_number = 'n019') # alignBins
  all_results <- applyOperation(all_results, config, op_number = 'n020') # buildConsensus
  all_results <- applyOperation(all_results, config, op_number = 'n021') # primerSeqErr
  all_results <- applyOperation(all_results, config, op_number = 'n022') # binSeqErr
  all_results <- applyOperation(all_results, config, op_number = 'n023') # fwd_alignBinsSP
  all_results <- applyOperation(all_results, config, op_number = 'n024') # rev_alignBinsSP
  all_results <- applyOperation(all_results, config, op_number = 'n025') # fwd_buildConsensus
  all_results <- applyOperation(all_results, config, op_number = 'n026') # rev_buildConsensus
  all_results <- applyOperation(all_results, config, op_number = 'n027') # binSeqErr_fwd_rev

  load('/fridge/data/MotifBinner2_test/binned/CAP008_1070_002wpi_V3B_C5A/CAP008_1070_002wpi_V3B_C5A.Rdata')
  
  
  op_number <- 'n021'
  config$current_op_number <- op_number

  timing <- list()
  ptm <- proc.time()
  op_number <- 'n015'
  config$current_op_number <- op_number
  op <- get(config$operation_list[[op_number]]$op)
  result <- op(all_results, config)
  seq_dat <- result$input_dat
  result$input_dat <- NULL
  timing$main <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- genSummary(result, config, seq_dat = seq_dat)
  timing$summary <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- computeMetrics(result, config, seq_dat = seq_dat)
  timing$metrics <- proc.time() - ptm
  
  ptm <- proc.time()
  result$timing <- timing
  result <- genReport(result, config)
  timing$report <- proc.time() - ptm
  
  ptm <- proc.time()
  result <- saveToDisk(result, config, seq_dat)
  timing$saveToDisk <- proc.time() - ptm

  ptm <- proc.time()

  print('bye')

#  config <- list(fwd_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R1.fastq",
#                 rev_reads_file = "/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R2.fastq",
#                 output_dir = "/fridge/data/MotifBinner2_test",
#                 prefix_for_names = "",
#                 fwd_primer = 'TATGGGAYSAAAGYCTMAARCCATGTG',
#                 rev_primer = 'CACACGCTCAGNNNNNNNNNATTCCATGTGTACATTGTACTGTRCTG',
#                 operation_list = c('loadData', 'basicQC', 'ambigSeqs', 'primerDimer', 'seqLength'),
#                 ambigSeqs = list(max_ambig = 5),
#                 primerDimer = list(primer_dimer_len = 80),
#                 seqLength = list(short_seq_len = 295,
#                                  long_seq_len = 305),
#                 intermediate_reports = TRUE,
#                 verbosity = 3,
#                 report_type = c('html'))
#
#  x <- do.call(processPrimers, config)
#  all_results <- x
#  all_results$summary <- NULL
#  
#  operation_function <- trimAffixes
#  config$operation_number <- length(all_results)
#  result <- operation_function(all_results, config)
#
#  str(result$fwd_result$kept$trim_stats)
#  print(result$fwd_result$kept$seq_dat@sread)
#
#  all_results <- list()
#  class(all_results) <- 'allResults'
#
#  dir.create(file.path(config$output_dir, config$prefix_for_names),
#             showWarnings = FALSE, recursive = TRUE)
#
#  all_results <- applyOperation('loadData', all_results, config)
#  all_results <- applyOperation('basicQC', all_results, config)
#  all_results <- applyOperation('ambigSeqs', all_results, config)
#  all_results <- applyOperation('primerDimer', all_results, config)
#  all_results <- applyOperation('seqLength', all_results, config)
#  all_results <- applyOperation('trimEnds', all_results, config)
#
#  setwd('~/projects/MotifBinner2/code/MotifBinner2')
#  unlink('/fridge/data/MotifBinner2_test/CAP256_3100_030wpi_v1v2_20k/n005_seqLength', recursive=T)
#  load_all(quiet=TRUE)
#
#  ptm <- proc.time()
#  timing <- list()
#  operation_function <- trimAffixes
#  config$operation_number <- length(all_results)
#
#  result <- operation_function(all_results, config)
#  timing$main <- proc.time() - ptm
#  ptm <- proc.time()
#
#  result <- saveToDisk(result, config)
#  timing$saveToDisk <- proc.time() - ptm
#
#  ptm <- proc.time()
#
#  result <- genSummary(result, config)
#  timing$genSummary <- proc.time() - ptm
#  ptm <- proc.time()
#
#  result <- computeMetrics(result, config)
#  timing$computeMetrics <- proc.time() - ptm
#  result$timing <- timing
#  ptm <- proc.time()
#
#  result <- genReport(result, config)
#  timing$genReport <- proc.time() - ptm
#  ptm <- proc.time()
#
#  result <- print(result, config)
#  timing$print <- proc.time() - ptm
#  result$timing <- timing
#
#  all_results[[basename(result$op_dir)]] <- result
#
#  return(all_results)
#
}

alignment_tests <- function()
{
  library(devtools)
  setwd('~/projects/MotifBinner2/code/MotifBinner2')
  Rcpp::compileAttributes()
#  load_all(quiet=TRUE)
  cat('\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n'); load_all()
  seq_dat <- readFastq('~/projects/MotifBinner2/code/MotifBinner2/inst/test_dat.fastq')
  seq_dat@sread

  prefix <- 'TATGGGAYSAAAGYCTMAARCCATGTG'

  trimmed <- trimFront_cpp(as.character(seq_dat@sread),
                          as.character(seq_dat@quality@quality),
                          prefix, c(13, 7, 7))



  trimmed$sread
}

memory_optimization <- function()
{
  for (i in 1:length(all_results)){
    res <- all_results[[i]]
    cat(paste(names(all_results)[i], ':\n', sep = ''))
    for (j in names(res)){
      size_rep <- round(pryr::object_size(res[j])/1000000,0)
      if (size_rep > 1){
        cat(paste(' - ', j, ' : ', size_rep, '\n', sep = ''))
      }
    }
  }
}

build_dependency_map <- function(operation_list){
  dependencies <- list()
  for (i in 1:length(operation_list)){
    dependencies[[names(operation_list)[i]]] <- NULL
    for (j in operation_list[[i]]$data_source){
      if (grepl("^n[0-9]*$", j)){
        dependencies[[j]] <- c(dependencies[[j]], names(operation_list)[i])
      }
    }
  }
  ops_performed <- gsub("_.*$", "", names(all_results))
  for (purge_step in names(dependencies)){
    if (all(dependencies[[purge_step]] %in% ops_performed)){
      all_results$seq_dat <- NULL
      all_results$trim_dat <- NULL
    } else {
    }
  }
}

purge_all_results <- function(all_results){
}
