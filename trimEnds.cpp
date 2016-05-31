#include <Rcpp.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <vector>

using namespace seqan;
using namespace Rcpp;

// [[Rcpp::export]]

void print_cpp(CharacterVector r_sread, CharacterVector r_id, CharacterVector r_qual)
{
  StringSet<IupacString> sq_sread;
  StringSet<CharString> sq_id;
  StringSet<CharString> sq_qual;

  sq_sread = Rcpp::as<std::vector<std::string> >(r_sread);
  sq_id = Rcpp::as<std::vector<std::string> >(r_id);
  sq_qual = Rcpp::as<std::vector<std::string> >(r_qual);

//  sq_id = r_id;
//  sq_qual = r_qual;

  std::cout << "inside print_cpp" << std::endl;
}

