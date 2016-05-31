library(Rcpp)
library(ShortRead)

seq_dat <- readFastq('CAP256_3100_030wpi_v1v2_20k_R1.fastq')

sourceCpp('trimEnds.cpp')

print(seq_dat)
print_cpp(as.character(seq_dat@sread),
          as.character(seq_dat@id),
          as.character(seq_dat@quality@quality))
