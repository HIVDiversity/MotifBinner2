#include <Rcpp.h>
#include <map>
#include <queue>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

std::map<char, int> getConsensusMatrixRowOrder()
{
  std::map<char, int> imap = 
  {
    {'A',  0},
    {'C',  1},
    {'G',  2},
    {'T',  3},
    {'M',  4},
    {'R',  5},
    {'W',  6},
    {'S',  7},
    {'Y',  8},
    {'K',  9},
    {'V', 10},
    {'H', 11},
    {'D', 12},
    {'B', 13},
    {'N', 14},
    {'-', 15},
    {'+', 16},
    {'.', 17}
  };
  return imap;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix gapQualityTweaker(CharacterVector reads, NumericMatrix q_mat)
{
  std::queue<int> last_sizes;
  last_sizes.push(0);
  double current_avg_quality = 0;
  int read_length;
  int jj;
  read_length = reads[1].size();
  bool read_has_begun;
  for (int i = 0; i < 7; ++i)
  {
    read_has_begun = false;
    current_avg_quality = 0;
    while (!last_sizes.empty()){
      last_sizes.pop();
    }
    if (i %2 == 0)
    {
      for (int j = 0; j < read_length; ++j)
      {
        // even number - fwd read
        jj = read_length - 1 - j;
        if (jj < 0) {
          throw std::range_error("fwd read tracing went out of bounds");
        }
        if (reads[i][jj] != '-' & !read_has_begun)
        {
          read_has_begun = true;
          for (int k = 0; k < 5; ++k)
          {
            last_sizes.push(q_mat(i,jj));
          }
          current_avg_quality = q_mat(i,jj);
        }
        if (read_has_begun & reads[i][jj] == '-')
        {
          q_mat(i,jj) = current_avg_quality;
        }
        if (read_has_begun & reads[i][jj] != '-')
        {
          last_sizes.push(q_mat(i,jj));
          current_avg_quality = std::max(0.0,
              current_avg_quality + last_sizes.back()/5.0 - last_sizes.front()/5.0);
          last_sizes.pop();
        }
      } 
    } else {
      for (int j = 0; j < read_length; ++j)
      {
//        std::cout << j << " " << current_avg_quality << " ";
        // odd number - rev read
        if (reads[i][j] != '-' & !read_has_begun)
        {
//          std::cout << "read begins now ";
          read_has_begun = true;
          for (int k = 0; k < 5; ++k)
          {
            last_sizes.push(q_mat(i,j));
          }
          current_avg_quality = q_mat(i,j);
        }
        if (read_has_begun & reads[i][j] == '-')
        {
//          std::cout << " gap detected, updating q_mat ";
          q_mat(i,j) = current_avg_quality;
        }
        if (read_has_begun & reads[i][j] != '-')
        {
          last_sizes.push(q_mat(i,j));
          current_avg_quality = std::max(0.0,
              current_avg_quality + last_sizes.back()/5.0 - last_sizes.front()/5.0);
//          std::cout << " base detected, updating curr_qual ";
//          std::cout << " new qual: " << q_mat(i,j) << " " << last_sizes.back() << 
//                       " old qual: " << last_sizes.front();
          last_sizes.pop();
        }
//        std::cout << std::endl;
      }
    }
  }
  return q_mat;
}

// [[Rcpp::export]]
Rcpp::List score_alignment_positions(CharacterVector reads, NumericMatrix q_mat)
{
  int k = 0;
  char c_let;
  int c_qual;
  int score_mat_row;
  std::map<char, int> score_mat_row_lookup;
  score_mat_row_lookup = getConsensusMatrixRowOrder();
  Rcpp::NumericMatrix score_mat (18, reads[1].size());
  Rcpp::NumericMatrix occurrence_mat (18, reads[1].size());
  for (int j = 0; j < reads[1].size(); ++j)
  {
    for (int i = 0; i < reads.size(); ++i)
    {
      c_let = reads[i][j];
      c_qual = q_mat(i,j);
      score_mat_row = score_mat_row_lookup[c_let];
      score_mat(score_mat_row, j) += c_qual;
      occurrence_mat(score_mat_row, j) ++;
//      std::cout << j << " " << i << " " << score_mat_row << " " << c_qual << std::endl;
    }
  }
//  std::cout << k << std::endl;

  Rcpp::List result;

  result = Rcpp::List::create(
    Rcpp::Named("score_mat") = score_mat,
    Rcpp::Named("occurrence_mat") = occurrence_mat
      );
  return result;
}
