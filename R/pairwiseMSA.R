#' Generates a MSA relative to a profile without any inserts into the profile
#'
#' Sequence quality information will be ignored.
#'
#' @param reads Reads either as a ShortRead or DNAStringSet object with the @sread and @id fields.
#' @param ref_seq The profile sequence to use as a reference
#' @export

pairwiseMSA <- function(reads, ref_seq){
  stuff <- pairWiseMSA_cpp(as.character(reads@sread), 
                           as.character(ref_seq),
                           as.character(reads@id))

  alned_seqs <- DNAStringSet(unlist(stuff$mapped))
  insert_map <- list2env(stuff$inserts)
  
  return(list(seqs = alned_seqs,
              insert_map = stuff$inserts))
}
