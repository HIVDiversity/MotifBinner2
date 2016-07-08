#include <Rcpp.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

struct dynMatEntry {
  int origin_i;
  int origin_j;
  int Fscore;
};

std::map<char, int> getBaseOrder()
{
  std::map<char, int> imap = 
  {
    {'U',  0},
    {'A',  1},
    {'C',  2},
    {'M',  3},
    {'G',  4},
    {'R',  5},
    {'S',  6},
    {'V',  7},
    {'T',  8},
    {'W',  9},
    {'Y', 10},
    {'H', 11},
    {'K', 12},
    {'D', 13},
    {'B', 14},
    {'N', 15}
  };
  return imap;
}

std::vector<std::vector<int> > getScoreMatrix()
{
  std::vector<std::vector<int> > smat(16, std::vector<int>(16));
  smat = {
      { 1,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1},
      { 0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1},
      { 0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  1,  1},
      { 0,  1,  1,  1,  0,  1,  1,  1,  0,  1,  1,  1,  0,  1,  1,  1},
      { 0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  1,  1,  1,  1},
      { 0,  1,  0,  1,  1,  1,  1,  1,  0,  1,  0,  1,  1,  1,  1,  1},
      { 0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1,  1,  1,  1},
      { 0,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1},
      { 1,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1},
      { 1,  1,  0,  1,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1},
      { 1,  0,  1,  1,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
      { 1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
      { 1,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
      { 1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
      { 1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
      { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1}
    };
  return smat;
}

void print_align_mat(const std::vector<std::vector<dynMatEntry> >& Fmat,
    const std::string & colNames,
    const std::string & rowNames)
{
  std::cout << "The dynamic programming matrix:" << std::endl;
  std::cout << ".";
  for (int i = 0; i != colNames.length(); i++)
  {
    std::cout << std::setw(3) << std::right << colNames[i];
  }
  std::cout << std::endl;
  for (int i = 0; i != Fmat.size(); i++)
  {
    std::cout << rowNames[i];
    for (int j = 0; j != Fmat[0].size(); j ++)
    {
      std::cout << std::setw(3) << std::right << Fmat[i][j].Fscore;
    }
    std::cout << std::endl;
  }
}

std::vector<std::vector<dynMatEntry> > 
initializeFmat(std::vector<std::vector<dynMatEntry> > Fmat,
    bool vseq_start_gap_pen,
    bool hseq_start_gap_pen,
    int vseq_len,
    int hseq_len)
{
  // The scores
  if (!vseq_start_gap_pen)
  {
    for (int j = 0; j != hseq_len; j++)
    { Fmat[0][j].Fscore = 0; }
  } else {
    for (int j = 0; j != hseq_len; j++)
    { Fmat[0][j].Fscore = -j; }
  }

  if (!hseq_start_gap_pen)
  {
    for (int i = 0; i != vseq_len; i++)
    { Fmat[i][0].Fscore = 0; }
  } else {
    for (int i = 0; i != vseq_len; i++)
    { Fmat[i][0].Fscore = -i; }
  }

  // The origins
  Fmat[0][0].origin_i = 0; 
  Fmat[0][0].origin_j = 0; 
  for (int i = 1; i != vseq_len; i++)
  { 
    Fmat[i][0].origin_i = i - 1; 
    Fmat[i][0].origin_j = 0; 
  }
  for (int j = 0; j != hseq_len; j++)
  { 
    Fmat[0][j].origin_i = 0; 
    Fmat[0][j].origin_j = j - 1; 
  }

  return Fmat;
}

std::vector<std::vector<dynMatEntry> > 
fillFmat(std::vector<std::vector<dynMatEntry> > Fmat,
    std::string vseq, 
    std::string hseq,
    int vseq_len,
    int hseq_len)
{
  std::map<char, int> imap;
  imap = getBaseOrder();
  std::vector<std::vector<int> > smat(16, std::vector<int>(16));
  smat = getScoreMatrix();

  // Fill out whole matrix
  int from_above = 0;
  int from_left = 0;
  int from_diag = 0;
  int max;
  int match_score;
  for (int i = 1; i != vseq_len; i++){
    for (int j = 1; j != hseq_len; j++){

      match_score = smat[imap[hseq[j]]][imap[vseq[i]]];
      from_above = Fmat[i - 1][j    ].Fscore - 2;
      from_left  = Fmat[i    ][j - 1].Fscore - 2;
      from_diag  = Fmat[i - 1][j - 1].Fscore + match_score;

      max = (from_above < from_left) ? from_left : from_above;
      max = (max < from_diag) ? from_diag : max;
      Fmat[i][j].Fscore = max;
      if (max == from_diag){
        Fmat[i][j].origin_i = i - 1;
        Fmat[i][j].origin_j = j - 1;
      } else if (max == from_above) {
        Fmat[i][j].origin_i = i - 1;
        Fmat[i][j].origin_j = j;
      } else {
        Fmat[i][j].origin_i = i;
        Fmat[i][j].origin_j = j - 1;
      }
    }
  }
  return Fmat;
}

std::pair<int, int> 
findBacktraceStart(std::vector<std::vector<dynMatEntry> > Fmat,
    bool hseq_end_gap_pen,
    bool vseq_end_gap_pen,
    int vseq_len,
    int hseq_len)
{
  // Find backtrace start position
  int i; int j;
  if (hseq_end_gap_pen and vseq_end_gap_pen)
  { // bottom right corner
    i = vseq_len - 1;
    j = hseq_len - 1;
  }
  int max;
  if (!hseq_end_gap_pen and vseq_end_gap_pen)
  { // last column highest score
    j = hseq_len - 1;
    max = -10000;
    for (int k = 0; k != vseq_len; k++)
    {
      if (Fmat[k][j].Fscore > max)
      {
        i = k;
        max = Fmat[k][j].Fscore;
      }
    }
  }
  if (hseq_end_gap_pen and !vseq_end_gap_pen)
  { // last row highest score
    i = vseq_len - 1;
    max = -10000;
    for (int m = 0; m != hseq_len; m++)
    {
      if (Fmat[i][m].Fscore > max)
      {
        j = m;
        max = Fmat[i][m].Fscore;
      }
    }
  }
  if (!hseq_end_gap_pen and !vseq_end_gap_pen)
  { // highest score of either last row or column
    int k = vseq_len - 1;
    max = -10000;
    for (int m = 0; m != hseq_len; m++)
    {
      if (Fmat[k][m].Fscore > max)
      {
        i = k;
        j = m;
        max = Fmat[k][m].Fscore;
      }
    }

    int m = hseq_len - 1;
    max = -10000;
    for (int k = 0; k != vseq_len; k++)
    {
      if (Fmat[k][m].Fscore > max)
      {
        i = k;
        j = m;
        max = Fmat[k][m].Fscore;
      }
    }
  }

  std::cout << "The starting position is: (" << i << ", " << j << "). The bottom left position is (" <<
    vseq_len - 1 << ", " << hseq_len - 1 << ")" << std::endl;
  return {i, j};
}

std::pair<std::string, std::string> 
fullBacktrace(int i, int j,
  int vseq_len, int hseq_len,
  std::string vseq, std::string hseq,
  std::vector<std::vector<dynMatEntry> > Fmat
  ) 
{
  // Backtrace
  std::string vseq_align;
  std::string hseq_align;

  // Navigate from bottom right to true start
  if (i != (vseq_len - 1) and j != (hseq_len - 1))
  {
    throw std::range_error("Back tracing must start in either the last row OR column - invalid scores?");
  }

  if (i != vseq_len - 1)
  {
    // starting in NOT in last row
    // navigate up
    // gap in hseq to base in vseq
    for (int k = vseq_len - 1; k != i; k--)
    {
      vseq_align += vseq[k];
      hseq_align += '-';
    }
  }
  if (j != hseq_len - 1)
  {
    // starting in NOT in last col
    // navigate left
    // gap in vseq to base in hseq
    for (int m = hseq_len - 1; m != j; m--)
    {
      vseq_align += '-';
      hseq_align += hseq[m];
    }
  }

  // Normal Backtrace
  while (i > 0 and j > 0)
  {
    if(Fmat[i-1][j-1].Fscore >= Fmat[i][j-1].Fscore and 
       Fmat[i-1][j-1].Fscore >= Fmat[i-1][j].Fscore)
    { // diag
      vseq_align += vseq[i];
      hseq_align += hseq[j];
      i--;
      j--;
      continue;
    }
    if (Fmat[i-1][j].Fscore >= Fmat[i][j-1].Fscore)
    { // from_above
      vseq_align += vseq[i];
      hseq_align += '-';
      i--;
      continue;
    }
    // from left
    vseq_align += '-';
    hseq_align += hseq[j];
    j--;
  }
  std::cout << "tracing endpoint: (" << i << ", " << j << ")" << std::endl;
  print_align_mat(Fmat, hseq, vseq);

  // Navigate from true end to the top left

  if (i > 0){
    for (int k = i; k != 0; k--)
    { // ended NOT in first ROW
      hseq_align += '-';
      vseq_align += vseq[k];
    }
  }
  if (j > 0){
    for (int m = j; m != 0; m--)
    { // ended NOT in first COLUMN
      hseq_align += hseq[m];
      vseq_align += '-';
    }
  }
  std::cout << "vseq: " << vseq_align << std::endl;
  std::cout << "hseq: " << hseq_align << std::endl;
  return {vseq_align, hseq_align};
}

std::vector<std::string> align(std::string vseq, std::string hseq,
    bool vseq_start_gap_pen,
    bool hseq_start_gap_pen,
    bool hseq_end_gap_pen,
    bool vseq_end_gap_pen)
{
  vseq = '-'+vseq;
  hseq = '-'+hseq;
  std::cout << "Seqs to align:" << std::endl;
  std::cout << "vseq: " << vseq << std::endl;
  std::cout << "hseq: " << hseq << std::endl;
  std::vector<std::string> alignment(2);

  int vseq_len = vseq.length();
  int hseq_len = hseq.length();

  std::vector<std::vector<dynMatEntry> > Fmat(vseq_len, std::vector<dynMatEntry>(hseq_len));
  Fmat = initializeFmat(Fmat, vseq_start_gap_pen, hseq_start_gap_pen, vseq_len, hseq_len);
  Fmat = fillFmat(Fmat, vseq, hseq, vseq_len, hseq_len);

  std::pair<int, int> backtrace_start;
  backtrace_start = findBacktraceStart(Fmat, hseq_end_gap_pen, vseq_end_gap_pen, vseq_len, hseq_len);

  std::pair<std::string, std::string> aligned_pair;
  aligned_pair = fullBacktrace(backtrace_start.first, 
    backtrace_start.second, 
    vseq_len, hseq_len, 
    vseq, hseq, Fmat);

  //reverse the alignments
  std::string vseq_align_fwd = aligned_pair.first;
  std::string hseq_align_fwd = aligned_pair.second;
  for (int i = aligned_pair.first.length()-1; i != -1; i--)
  {
    vseq_align_fwd[aligned_pair.first.length() - i -1] = aligned_pair.first[i];
    hseq_align_fwd[aligned_pair.first.length() - i -1] = aligned_pair.second[i];
  }

  alignment[0] = vseq_align_fwd;
  alignment[1] = hseq_align_fwd;
  return alignment;
}

// [[Rcpp::export]]
Rcpp::List trimFront_cpp(CharacterVector r_sread, CharacterVector r_qual,
    CharacterVector r_primer, int pref_len, int pid_len, int suf_len, int verbosity)
{
  std::vector<std::string> alignment(2);
  std::string read_start;
  for (int i=0; i!= r_sread.size(); i++)
  {
    read_start = Rcpp::as<std::string>(r_sread[i]);
    read_start = read_start.substr(0, 20+pref_len+pid_len+suf_len);

    alignment = align(read_start,
                      Rcpp::as<std::string>(r_primer),
                      false, false, false, false);

    if (verbosity > 1)
    {
      std::cout << "Sequence Number: " << i << std::endl;
      std::cout << alignment[0] << std::endl;
      std::cout << alignment[1] << std::endl;
      std::cout << "*****************************" << std::endl;
      std::cout << "\n\n\n\n\n\n\n\n\n\n\n" << std::endl;
    }
  }
  Rcpp::List trim_result;

  trim_result = Rcpp::List::create(Rcpp::Named("sread") = alignment);

  return trim_result;
}

