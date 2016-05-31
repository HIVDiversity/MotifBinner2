#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]

void print_cpp(CharacterVector sread, CharacterVector seq_id, CharacterVector qual)
{
  std::cout << "inside print_cpp" << std::endl;
}

