#include <Rcpp.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <vector>

using namespace seqan;
using namespace Rcpp;

void findPrefixMatch(StringSet<IupacString> haystack,
    StringSet<CharString> id, StringSet<CharString> qual,
    IupacString needle)
{
  int seq_len = 0;
  Align<IupacString, ArrayGaps> align;
  resize(rows(align), 2);
  Row<Align<IupacString, ArrayGaps> >::Type & row0 = row(align, 0);
  Row<Align<IupacString, ArrayGaps> >::Type & row1 = row(align, 1);

  for (unsigned i = 0; i < length(haystack); ++i)
  {
    seq_len = length(haystack[i]);
    assignSource(row(align, 0), infix(haystack[i], 0, std::min(50, seq_len)));
    assignSource(row(align, 1), needle);
    std::cout << "  " << i;
  }
  std::cout << std::endl;
}

// [[Rcpp::export]]

void trimEnds_cpp(CharacterVector r_sread, CharacterVector r_id, 
    CharacterVector r_qual, CharacterVector r_prefix)
{
  StringSet<IupacString> sq_sread;
  StringSet<CharString> sq_id;
  StringSet<CharString> sq_qual;
  StringSet<IupacString> sq_prefix;

  sq_sread = Rcpp::as<std::vector<std::string> >(r_sread);
  sq_id = Rcpp::as<std::vector<std::string> >(r_id);
  sq_qual = Rcpp::as<std::vector<std::string> >(r_qual);
  sq_prefix = Rcpp::as<std::vector<std::string> >(r_prefix);

  std::cout << "length of haystack " << length(sq_sread) << std::endl;
  std::cout << "length of prefix " << length(sq_prefix) << std::endl;
  std::cout << "length of first element of prefix " << length(sq_prefix[0]) << std::endl;

  findPrefixMatch(sq_sread, sq_id, sq_qual, sq_prefix[0]);


//  sq_id = r_id;
//  sq_qual = r_qual;

  std::cout << "inside print_cpp" << std::endl;
}

