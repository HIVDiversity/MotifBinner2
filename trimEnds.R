library(Rcpp)
library(ShortRead)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

seq_dat <- readFastq('/fridge/data/MotifBinner2_test/raw/CAP256_3100_030wpi_v1v2_R1.fastq')
seq_dat <- readFastq('test_dat.fastq')
seq_dat <- seq_dat[1:4]

sourceCpp('trimEnds.cpp')
prefix <- "ATATATATATATATATATATATATATA"
prefix <- "TATGGGAYSAAAGYCTMAARCCATGTG"

system.time(
trimmed <- 
trimEnds_cpp(as.character(seq_dat@sread),
             as.character(seq_dat@id),
             as.character(seq_dat@quality@quality),
             prefix)
)

trimmed <- ShortReadQ(sread = DNAStringSet(trimmed$sread),
                      quality = BStringSet(trimmed$qual),
                      id = BStringSet(trimmed$id))
print(seq_dat)
print(seq_dat@sread)
print(trimmed@sread)


