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
  IupacString sq_sread;
  CharString sq_id;
  CharString sq_qual;

//  sq_sread = r_sread;
//  sq_id = r_id;
//  sq_qual = r_qual;

  std::cout << "inside print_cpp" << std::endl;
}

