library(Rcpp)
library(ShortRead)

seq_dat <- readFastq('test_dat.fastq')
#seq_dat <- readFastq('/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_20k_R1.fastq')
seq_dat <- readFastq('/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq')

prefix <- 'TATGGGAYSAAAGYCTMAARCCATGTG'

print(seq_dat@sread)

cat('\n\n\n\n\n\n\n\n\n'); sourceCpp('trimFront.cpp')
trimmed <- trimFront_cpp(as.character(seq_dat@sread),
                        as.character(seq_dat@quality@quality),
                        prefix, c(13, 7, 7))

str(trimmed)

