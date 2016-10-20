#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option("--fwd_file", help = "The file containing the FWD reads"),
make_option("--fwd_primer_seq", 
            help = "The full sequence of the FWD primer - include degeneracies used for phasing"),
make_option("--fwd_primer_lens",
            help = paste("The lengths of the different fragments of the fwd primer. ",
                         "It is only important to specify the primer in more than one fragment if ",
                         "the PID is included in this primer. Provide a comma seperated list of ",
                         "integers - no spaces!! example: 12,9,15",
                         sep = "")),
make_option("--fwd_primer_min_score",
            help = paste("The minimum alignment score between the primer sequence and the read. ",
                         "A match counts 1, a gap -1 and a mismatch 0. Note that any matches to ",
                         "ambiguous letters also count 1, so each N in the primer guarentees a ",
                         "score of 1 for that base.",
                         sep = "")),
make_option("--fwd_pid_in_which_fragment",
            help = paste("The index of the fragment of the FWD primer containing the PID. ",
                         "NULL if this primer does not contain a PID. The first index is ",
                         "1 - not 0. ",
                         sep = "")),

make_option("--rev_file", help = "The file containing the REV reads"),
make_option("--rev_primer_seq", 
            help = "The full sequence of the REV primer - include degeneracies used for phasing"),
make_option("--rev_primer_lens",
            help = paste("The lengths of the different fragments of the rev primer. ",
                         "It is only important to specify the primer in more than one fragment if ",
                         "the PID is included in this primer. Provide a comma seperated list of ",
                         "integers - no spaces!! example: 12,9,15",
                         sep = "")),
make_option("--rev_primer_min_score",
            help = paste("The minimum alignment score between the primer sequence and the read. ",
                         "A match counts 1, a gap -1 and a mismatch 0. Note that any matches to ",
                         "ambiguous letters also count 1, so each N in the primer guarentees a ",
                         "score of 1 for that base.",
                         sep = "")),
make_option("--rev_pid_in_which_fragment",
            help = paste("The index of the fragment of the REV primer containing the PID. ",
                         "NULL if this primer does not contain a PID. The first index is ",
                         "1 - not 0. ",
                         sep = "")),

make_option("--min_read_length",
            default = 295,
            help = paste("Require reads to be longer than this in the seqLength trimming step. ",
                         "It is recommended to set a strict value for Illumina machines - ",
                         "consider expected read length - 5.",
                         sep = "")),
#make_option("--profile_file",
#            help = paste("The file containing the profile to use for mapping the reads in the ", 
#                         "bins. Specifiy this option if you want both FWD and REV reads to be ",
#                         "mapped to the same profile. Always manually inspect some of the ",
#                         "aligned bins if you have no or small overlaps or if you are sequencing ",
#                         "a highly variable sample - like variable loops in ENV",
#                         sep = "")),
make_option("--output_dir",
            help = paste("Path to the folder where the results are to be stored. ",
                         "Note that a subfolder will be made in this directory",
                         sep = "")),
make_option("--base_for_names",
            help = paste("This will be prepended to the names all output files.",
                         sep = "")),
make_option("--ncpu",
            help = paste("Number of cores to use. NOTE: Each additional core specified ",
                         "significantly increases memory consumption.",
                         sep = "")))

opt <- parse_args(OptionParser(option_list = option_list,
  description = "Bin Illumina reads produced with Primer ID approach",
  epilogue = "Example Call:
./MotifBinner2.R --fwd_file=/fridge/data/molecular_clock2/raw/HVTN503_159451817_1012_C2C3_R1.fastq --fwd_primer_seq=NNNNNNNCTCTTTTGACCCAATTCCTATACATTATTG --fwd_primer_lens=37 --fwd_primer_min_score=32 --rev_file=/fridge/data/molecular_clock2/raw/HVTN503_159451817_1012_C2C3_R2.fastq --rev_primer_seq=NNNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT --rev_primer_lens=15,28 --rev_primer_min_score=37 --fwd_pid_in_which_fragment=NULL --rev_pid_in_which_fragment=1 --output_dir=/fridge/data/molecular_clock2/pipeline1 --base_for_names=HVTN503_159451817_1012_C2C3 --ncpu=6 "))

#######################
# old commands examples
#######################
# ?unknown
#./MotifBinner2.R --fwd_file=/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R1.fastq --fwd_primer_seq=CTCTTTTGACCCAATTCCTATACATTATTG --fwd_primer_lens=30 --fwd_primer_min_score=26 --rev_file=/fridge/data/MotifBinner2_test/raw/CAP129_2040_009wpi_C2C3_R2.fastq --rev_primer_seq=NNNNNNNNNNNNNNTGCAATAGAAAAATTCTCCTCTACAATT --rev_primer_lens=14,28 --rev_primer_min_score=36 --fwd_pid_in_which_fragment=NULL --rev_pid_in_which_fragment=1 --output_dir=/fridge/data/MotifBinner2_test --base_for_names=CAP129_2040_009wpi_C2C3 --ncpu=6 "))

# CAP256
#./MotifBinner2.R --fwd_file=/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq --fwd_primer_seq=TATGGGAYSAAAGYCTMAARCCATGTG --fwd_primer_lens=27 --fwd_primer_min_score=22 --rev_file=/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R2.fastq --rev_primer_seq=CACACGCTCAGNNNNNNNNNATTCCATGTGTACATTGTACTGTRCTG --rev_primer_lens=11,9,27 --rev_primer_min_score=42 --fwd_pid_in_which_fragment=NULL --rev_pid_in_which_fragment=2 --output_dir=/fridge/data/MotifBinner2_test --base_for_names=CAP256_3100_030wpi_v1v2 --ncpu=6


# still used profiles
# ./MotifBinner2.R --fwd_file=/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq --fwd_primer_seq=TATGGGAYSAAAGYCTMAARCCATGTG --fwd_primer_lens=27 --fwd_primer_min_score=22 --rev_file=/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R2.fastq --rev_primer_seq=CACACGCTCAGNNNNNNNNNATTCCATGTGTACATTGTACTGTRCTG --rev_primer_lens=11,9,27 --rev_primer_min_score=42 --fwd_pid_in_which_fragment=NULL --rev_pid_in_which_fragment=2 --profile_file=/fridge/data/MotifBinner2_test/v1v2_profile1.fasta --output_dir=/fridge/data/MotifBinner2_test --base_for_names=CAP256_3100_030wpi_v1v2 --ncpu=6

suppressPackageStartupMessages(library(MotifBinner2))

buildConfig <- function(fwd_file, fwd_primer_seq, fwd_primer_lens, fwd_min_score,
                        rev_file, rev_primer_seq, rev_primer_lens, rev_min_score,
                        fwd_pid_in_which_fragment, rev_pid_in_which_fragment,
                        profile_file, fwd_profile_file = NULL, rev_profile_file = NULL,
                        min_read_length = 295,
                        pattern_to_chop_from_names = ' [0-9]:N:[0-9]*:[0-9]*$',
                        output_dir = "/fridge/data/MotifBinner2_test",
                        base_for_names = "CAP129_2040_009wpi_C2C3",
                        intermediate_reports = TRUE,
                        erase_history = FALSE,
                        verbosity = 3,
                        report_type = c('html'),
                        ncpu = 4,
                        bins_to_process = Inf
                        )
{
  if (fwd_pid_in_which_fragment == "NULL"){fwd_pid_in_which_fragment <- NULL}
  if (rev_pid_in_which_fragment == "NULL"){rev_pid_in_which_fragment <- NULL}
  operation_list = list(
    'n001' = 
      list(name = 'fwd_loadData',
        op = 'loadData',
        data_source = fwd_file,
        cache_data = TRUE),
    'n002' =
      list(name = 'fwd_basicQC',
        op = 'basicQC',
        data_source = "n001",
        cache_data = FALSE),
    'n003' =
      list(name = 'fwd_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n001",
        threshold = 0.02,
        cache_data = TRUE),
    'n004' =
      list(name = 'fwd_primerDimer',
        op = 'primerDimer',
        data_source = "n003",
        threshold = 80,
        cache_data = TRUE),
    'n005' =
      list(name = 'fwd_seqLength',
        op = 'seqLength',
        data_source = "n004",
        threshold = min_read_length,
        cache_data = TRUE),
    'n006' =
      list(name = 'fwd_qualTrim',
        op = 'qualTrim',
        data_source = "n005",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.15,
        cache_data = TRUE),
    'n007' =
      list(name = 'fwd_trimAffixes',
        op = 'trimAffixes',
        data_source = "n006",
        
        primer_seq = fwd_primer_seq,
        primer_lens = fwd_primer_lens,
        min_score = fwd_min_score,

        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n008' = 
      list(name = 'rev_loadData',
        op = 'loadData',
        data_source = rev_file,
        cache_data = TRUE),
    'n009' =
      list(name = 'rev_basicQC',
        op = 'basicQC',
        data_source = "n008",
        cache_data = FALSE),
    'n010' =
      list(name = 'rev_ambigSeqs',
        op = 'ambigSeqs',
        data_source = "n008",
        threshold = 0.02,
        cache_data = TRUE),
    'n011' =
      list(name = 'rev_primerDimer',
        op = 'primerDimer',
        data_source = "n010",
        threshold = 80,
        cache_data = TRUE),
    'n012' =
      list(name = 'rev_seqLength',
        op = 'seqLength',
        data_source = "n011",
        threshold = min_read_length,
        cache_data = TRUE),
    'n013' =
      list(name = 'rev_qualTrim',
        op = 'qualTrim',
        data_source = "n012",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.15,
        cache_data = TRUE),
    'n014' =
      list(name = 'rev_trimAffixes',
        op = 'trimAffixes',
        data_source = "n013",
        
        primer_seq = rev_primer_seq,
        primer_lens = rev_primer_lens,
        min_score = rev_min_score,
        
        primer_location = 'front',
        front_gaps_allowed = 0,
        cache_data = TRUE),
    'n015' =
      list(name = 'fwd_extractPIDs',
        op = 'extractPIDs',
        data_source = "n007",
        pid_in_which_fragment = fwd_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n016' =
      list(name = 'rev_extractPIDs',
        op = 'extractPIDs',
        data_source = "n014",
        pid_in_which_fragment = rev_pid_in_which_fragment,
        pattern_to_chop_from_names = pattern_to_chop_from_names,
        pid_gaps_allowed = 0,
        cache_data = TRUE),
    'n017' =
      list(name = 'matchPairs',
        op = 'matchPairs',
        data_source = c("fwd" = "n015", "rev" = "n016"),
        cache_data = TRUE),
    'n018' =
      list(name = 'processBadPIDs',
        op = 'processBadPIDs',
        data_source = "n017",
        cache_data = TRUE),
    'n019' =
      list(name = 'mergePEAR',
        op = 'mergePEAR',
        data_source = "n018",
        cache_data = TRUE),
    'n020' =
      list(name = 'merge_qualTrim',
        op = 'qualTrim',
        data_source = "n019",
        avg_qual = 20,
        bad_base_threshold = 10,
        max_bad_bases = 0.05,
        cache_data = TRUE),
    'n021' =
      list(name = 'binSizeCheck',
        op = 'binSizeCheck',
        data_source = "n020",
        min_bin_size = 3,
        cache_data = TRUE),
    'n022' =
      list(name = 'alignBinsMSA',
        op = 'alignBinsMSA',
        bins_to_process = bins_to_process,
        data_source = "n021",
        cache_data = TRUE),
    'n023' =
      list(name = 'buildConsensus',
        op = 'buildConsensus',
        data_source = "n022",
        cache_data = TRUE),
    'n024' =
      list(name = 'primerSeqErr',
        op = 'primerSeqErr',
        data_source = c("fwd" = "n007", "rev" = "n014"),
        cache_data = FALSE),
    'n025' =
      list(name = 'binSeqErr',
        op = 'binSeqErr',
        data_source = c("bin_msa_merged" = "n022", "cons_merged" = "n023", "primer_err" = "n024"),
        cache_data = FALSE),
    'n100' =
      list(name = 'dataTracing',
        op = 'dataTracing',
        data_source = c(
          "fwdReads.1" = "n001", "fwdReads.2" = "n003", "fwdReads.3" = "n004",
          "fwdReads.4" = "n005", "fwdReads.5" = "n006", "fwdReads.6" = "n007",

          "revReads.1" = "n008", "revReads.2" = "n010", "revReads.3" = "n011",
          "revReads.4" = "n012", "revReads.5" = "n013", "revReads.6" = "n014",

          "mergeReads.1" = "n017", "mergeReads.2" = "n018", "mergeReads.3" = "n019",
          "mergeReads.4" = "n020", "mergeReads.5" = "n021", "mergeReads.6" = "n022",
          "mergeReads.7" = "n023"),
      cache_data = FALSE)
  )

  return(list(operation_list = operation_list,
              output_dir = output_dir,
              base_for_names = base_for_names,
              intermediate_reports = intermediate_reports,
              verbosity = verbosity,
              erase_history = erase_history,
              report_type = report_type,
              ncpu = ncpu))
}

print(opt)

config <-
buildConfig(fwd_file = opt$fwd_file,
            fwd_primer_seq = opt$fwd_primer_seq,
            fwd_primer_lens = as.numeric(strsplit(opt$fwd_primer_lens, ',')[[1]]),
            fwd_min_score = as.numeric(opt$fwd_primer_min_score),
            rev_file = opt$rev_file,
            rev_primer_seq = opt$rev_primer_seq,
            rev_primer_lens = as.numeric(strsplit(opt$rev_primer_lens, ',')[[1]]),
            rev_min_score = as.numeric(opt$rev_primer_min_score),
            fwd_pid_in_which_fragment = ifelse(opt$fwd_pid_in_which_fragment == "NULL", 
                                               "NULL", 
                                               as.numeric(opt$fwd_pid_in_which_fragment)),
            rev_pid_in_which_fragment = ifelse(opt$rev_pid_in_which_fragment == "NULL", 
                                               "NULL", 
                                               as.numeric(opt$rev_pid_in_which_fragment)),
            min_read_length = as.numeric(opt$min_read_length),
            output_dir = opt$output_dir,
            base_for_names = opt$base_for_names,
            ncpu = opt$ncpu)

print(config)

unlink(file.path(config$output_dir, config$base_for_names), recursive=T)
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
all_results <- applyOperation(all_results, config, op_number = 'n019') # merge
all_results <- applyOperation(all_results, config, op_number = 'n020') # qualTrim
all_results <- applyOperation(all_results, config, op_number = 'n021') # binSizeCheck
all_results <- applyOperation(all_results, config, op_number = 'n022') # alignBinsMSA
all_results <- applyOperation(all_results, config, op_number = 'n023') # buildConsensus
all_results <- applyOperation(all_results, config, op_number = 'n024') # primerSeqErr
all_results <- applyOperation(all_results, config, op_number = 'n025') # binSeqErr

all_results <- applyOperation(all_results, config, op_number = 'n100') # dataTracing

x <- genReport(all_results, config)

save.image(file.path(opt$output_dir, 
                     opt$base_for_names, 
                     paste(opt$base_for_names, 
                           '.Rdata', 
                           sep = '')))

