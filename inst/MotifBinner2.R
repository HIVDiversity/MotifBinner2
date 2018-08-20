#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option("--non_overlapping", 
            action = "store_true",
            default = FALSE,
            help = "If your forwards and reverse reads do not overlap, include this flag. Also, if they do overlap, but you do not want them to be merged, include this flag."),
make_option("--overlapping", 
            action = "store_true",
            default = FALSE,
            help = "If your forwards and reverse reads overlap and you want them to be merged, include this flag."),
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
make_option("--max_seq",
            default = NULL,
            help = paste("Maximum number of sequences to read from the input files. By default the ",
                         "forward and reverse read files produced by a MiSeq are in the same order, ",
                         "so just reading out the first 'max_seq' sequences from each file should ",
                         "yield 'max_seq' pair-end reads in the dataset.",
                         sep = "")),
make_option("--min_read_length",
            default = 295,
            help = paste("Require reads to be longer than this in the seqLength trimming step. ",
                         "It is recommended to set a strict value for Illumina machines - ",
                         "consider expected read length - 5.",
                         sep = "")),
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
                         sep = "")),
make_option("--merged_read_length",
            default = 100,
            help = paste("Require reads to be longer than this in the seqLength trimming step AFTER building consensus sequences ",
                         sep = "")),
make_option("--header_format",
            default = "MiSeq",
            help = paste("Format of the sequence names. Set to 'SRA' if downloaded from NCBI's SRA. ",
                         sep = "")),
make_option("--front_gaps_allowed",
            default = 0,
            help = paste("Number of non-primer letters allowed before the specified primer sequence in the read. Strongly recommend 0, unless you have a variable length primer. ",
                         sep = "")),
make_option("--erase_history", 
            action = "store_true",
            default = FALSE,
            help = "Delete all data from memory that is not needed for subsequent steps.")
)

opt <- parse_args(OptionParser(option_list = option_list,
  description = "Bin Illumina reads produced with Primer ID approach",
  epilogue = "Example Call:
./MotifBinner2.R --non_overlapping --fwd_file=/fridge/data/walter_nol_2016_11/raw/CAP188_4220_004wpi_v1v3_R1.fastq --fwd_primer_seq=NNNNNNNNCAAAGYCTAAARCCATGTGTA --fwd_primer_lens=29 --fwd_primer_min_score=25 --rev_file=/fridge/data/walter_nol_2016_11/raw/CAP188_4220_004wpi_v1v3_R2.fastq --rev_primer_seq=TGACNNNNNNNNNNNNNNNTACTAATGTTACAATRTGC --rev_primer_lens=4,15,19 --rev_primer_min_score=34 --fwd_pid_in_which_fragment=NULL --rev_pid_in_which_fragment=2 --output_dir=/fridge/data/walter_nol_2016_11 --base_for_names=CAP188_4220_004wpi_v1v3 --ncpu=6 --merged_read_length=250"))

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

print(opt)
if ((!opt$overlapping) & (!opt$non_overlapping)){
  stop('Must specify either --overlapping or --non_overlapping')
} else if (opt$overlapping & opt$non_overlapping){
  stop('Must specify ONLY ONE of --overlapping or --non_overlapping')
} else if (opt$non_overlapping) {
  opt$overlapping <- FALSE
}

fwd_primer_length <- nchar(opt$fwd_primer_seq)
print(strsplit(opt$fwd_primer_lens, ',')[[1]])
fwd_spec_length <- sum(as.numeric(strsplit(opt$fwd_primer_lens, ',')[[1]]))
if (fwd_primer_length != fwd_spec_length){
  stop('Forward primer does not match specified lengths')
}
rev_primer_length <- nchar(opt$rev_primer_seq)
rev_spec_length <- sum(as.numeric(strsplit(opt$rev_primer_lens, ',')[[1]]))
if (rev_primer_length != rev_spec_length){
  stop('Reverse primer does not match specified lengths')
}
if (!file.exists(opt$fwd_file)){stop('Error: Forward file not found - are you sure the file name you specified is correct?')}
if (!file.exists(opt$rev_file)){stop('Error: Reverse file not found - are you sure the file name you specified is correct?')}

print(as.numeric(opt$max_seq))

config <-
buildConfig(overlapping = opt$overlapping,
            fwd_file = opt$fwd_file,
            fwd_primer_seq = opt$fwd_primer_seq,
            fwd_primer_lens = as.numeric(strsplit(opt$fwd_primer_lens, ',')[[1]]),
            fwd_min_score = as.numeric(opt$fwd_primer_min_score),
            rev_file = opt$rev_file,
            rev_primer_seq = opt$rev_primer_seq,
            rev_primer_lens = as.numeric(strsplit(opt$rev_primer_lens, ',')[[1]]),
            rev_min_score = as.numeric(opt$rev_primer_min_score),
            max_seq = as.numeric(opt$max_seq),
            fwd_pid_in_which_fragment = ifelse(opt$fwd_pid_in_which_fragment == "NULL", 
                                               "NULL", 
                                               as.numeric(opt$fwd_pid_in_which_fragment)),
            rev_pid_in_which_fragment = ifelse(opt$rev_pid_in_which_fragment == "NULL", 
                                               "NULL", 
                                               as.numeric(opt$rev_pid_in_which_fragment)),
            min_read_length = as.numeric(opt$min_read_length),
            output_dir = opt$output_dir,
            base_for_names = opt$base_for_names,
            ncpu = opt$ncpu,
            merged_read_length = opt$merged_read_length,
            header_format = opt$header_format,
            front_gaps_allowed = opt$front_gaps_allowed,
            erase_history = opt$erase_history)

dput(config)

unlink(file.path(config$output_dir, config$base_for_names), recursive=T)
all_results <- list()
class(all_results) <- 'allResults'

for (op_num in names(config$operation_list)){
  all_results <-  applyOperation(all_results, config, op_number = op_num)
}

file.copy(from = file.path(all_results$n027_cons_ambigSeqs$config$op_dir, 
                           paste(config$base_for_names, "kept", "cons_ambigSeqs.fastq", sep = "_")),
          to = file.path(config$output_dir, config$base_for_names))

x <- genReport(all_results, config)

save.image(file.path(opt$output_dir, 
                     opt$base_for_names, 
                     paste(opt$base_for_names, 
                           '.Rdata', 
                           sep = '')))
