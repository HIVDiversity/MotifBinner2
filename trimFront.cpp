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

struct fragmentStore {
  std::string read_segment;
  int read_gaps;
  int read_bases;
  std::string prefix_segment;
  int prefix_gaps;
  int prefix_bases;
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

//  std::cout << "The starting position is: (" << i << ", " << j << "). The bottom right position is (" <<
//    vseq_len - 1 << ", " << hseq_len - 1 << ")" << std::endl;
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

//  // Wrong Normal Backtrace
//  while (i > 0 and j > 0)
//  {
//    if(Fmat[i-1][j-1].Fscore >= Fmat[i][j-1].Fscore and 
//       Fmat[i-1][j-1].Fscore >= Fmat[i-1][j].Fscore)
//    { // diag
//      vseq_align += vseq[i];
//      hseq_align += hseq[j];
//      i--;
//      j--;
//      continue;
//    }
//    if (Fmat[i-1][j].Fscore >= Fmat[i][j-1].Fscore)
//    { // from_above
//      vseq_align += vseq[i];
//      hseq_align += '-';
//      i--;
//      continue;
//    }
//    // from left
//    vseq_align += '-';
//    hseq_align += hseq[j];
//    j--;
//  }
//  std::cout << "tracing endpoint: (" << i << ", " << j << ")" << std::endl;
//  print_align_mat(Fmat, hseq, vseq);

  // Retrieve 'Stored' Backtrace
  int new_i;
  int new_j;
  while (i > 0 and j > 0)
  {
//    std::cout << "Current: (" << i << ", " << j << ")" << 
//      "Next: (" << Fmat[i][j].origin_i << ", " << 
//      Fmat[i][j].origin_j << ")" << std::endl;
//    std::cout << "vseq: " << vseq_align << std::endl;
//    std::cout << "hseq: " << hseq_align << std::endl;
    if (Fmat[i][j].origin_i != i and Fmat[i][j].origin_j != j)
    { //diag
      vseq_align += vseq[i];
      hseq_align += hseq[j];
    } else if (Fmat[i][j].origin_i != i)
    { //from_above
      vseq_align += vseq[i];
      hseq_align += '-';
    } else if (Fmat[i][j].origin_j != j)
    { //from_left
      vseq_align += '-';
      hseq_align += hseq[j];
    } else {
      throw std::range_error("Back tracing inifite loop issues");
    }
    new_i = Fmat[i][j].origin_i;
    new_j = Fmat[i][j].origin_j;
    i = new_i;
    j = new_j;
  }

//  std::cout << "tracing endpoint: (" << i << ", " << j << ")" << std::endl;
//  print_align_mat(Fmat, hseq, vseq);
  
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
//  std::cout << "vseq: " << vseq_align << std::endl;
//  std::cout << "hseq: " << hseq_align << std::endl;
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
//  std::cout << "Seqs to align:" << std::endl;
//  std::cout << "vseq: " << vseq << std::endl;
//  std::cout << "hseq: " << hseq << std::endl;
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
    CharacterVector r_primer, std::vector<int> prefix_lens) // int pref_len, int pid_len, int suf_len, int verbosity)
{

  std::vector<std::string> alignment(2);
  std::string read_start;

  int read_bases;
  int prefix_bases;
  int read_gaps;
  int prefix_gaps;

  int prev_read_bases;
  int prev_prefix_bases;
  int prev_read_gaps;
  int prev_prefix_gaps;

  int cur_prefix_fragment;
  int next_prefix_break;
  int cur_fragment_start_in_aln;
  std::string read_segment;
  std::string prefix_segment;

  std::vector<int> read_front_gaps(r_sread.size());
  std::vector<int> prefix_front_gaps(r_sread.size());

  std::vector<std::vector<int> > fragment_ints(r_sread.size(), 
      std::vector<int>(prefix_lens.size()*4));
  std::vector<std::vector<std::string> > fragment_strings(r_sread.size(), 
      std::vector<std::string>(prefix_lens.size()*2));

  for (int i=0; i!= r_sread.size(); i++)
  {
    read_start = Rcpp::as<std::string>(r_sread[i]);
    read_start = read_start.substr(0, 20+Rcpp::as<std::string>(r_primer[0]).length());

    alignment = align(read_start,
                      Rcpp::as<std::string>(r_primer),
                      false, false, false, false);

    if (alignment[0].length() != alignment[1].length())
    {
      throw std::range_error("Alignments must be of equal lengths");
    }

    std::string& read = alignment[0]; // reference for easier name
    std::string& prefix = alignment[1]; // reference for easier name
    read_gaps = 0; read_bases = 0;
    prefix_gaps = 0; prefix_bases = 0;
    cur_prefix_fragment = 0;
    next_prefix_break = prefix_lens[0];

    read_front_gaps[i] = -1; prefix_front_gaps[i] = -1;

    for (int j = 0; j < read.size(); ++j)
    {
      // basic tracking
      if (read[j] == '-'){ ++read_gaps; } else { ++read_bases;}
      if (prefix[j] == '-'){ ++prefix_gaps; } else { ++prefix_bases; }
      
      // once per alignment
      if (read_front_gaps[i] == -1 and read_bases == 1) {read_front_gaps[i] = read_gaps;}
      if (prefix_front_gaps[i] == -1 and prefix_bases == 1) 
      {
        prefix_front_gaps[i] = prefix_gaps;
        cur_fragment_start_in_aln = j;
        if (read[j] == '-')
        {
          prev_read_gaps = read_gaps - 1;
          prev_read_bases = read_bases;
        } else {
          prev_read_gaps = read_gaps;
          prev_read_bases = read_bases - 1;
        }

        prev_prefix_gaps = prefix_gaps;
        prev_prefix_bases = prefix_bases - 1;
      }

      // once for each prefix segment
      if (prefix_bases >= next_prefix_break)
      {
        read_segment = read.substr(cur_fragment_start_in_aln, j - cur_fragment_start_in_aln);
        prefix_segment = prefix.substr(cur_fragment_start_in_aln, j - cur_fragment_start_in_aln);
        cur_fragment_start_in_aln = j;

        fragment_ints[i][cur_prefix_fragment*4 + 0] = read_gaps - prev_read_gaps;
        fragment_ints[i][cur_prefix_fragment*4 + 1] = read_bases - prev_read_bases;
        fragment_ints[i][cur_prefix_fragment*4 + 2] = prefix_gaps - prev_prefix_gaps;
        fragment_ints[i][cur_prefix_fragment*4 + 3] = prefix_bases - prev_prefix_bases;

        fragment_strings[i][cur_prefix_fragment*2 + 0] = read_segment;
        fragment_strings[i][cur_prefix_fragment*2 + 1] = prefix_segment;
        
        cur_prefix_fragment ++;
        if (cur_prefix_fragment < prefix_lens.size()){
          next_prefix_break += prefix_lens[cur_prefix_fragment];
        } else {
          next_prefix_break += 10000;
        }
        
        prev_read_gaps = read_gaps;
        prev_read_bases = read_bases;
        prev_prefix_gaps = prefix_gaps;
        prev_prefix_bases = prefix_bases;
      }
    }
  }

  for (int i = 0; i != read_front_gaps.size(); ++i)
  {

    std::cout << "Sequence Number: " << i << std::endl;
    std::cout << "Read Front Gaps: " << read_front_gaps[i] << ";   Prefix Front Gaps: " 
      << prefix_front_gaps[i] << std::endl;
    for (int j = 0; j != prefix_lens.size(); ++j)
    {
        std::cout << "Read Segment: " << fragment_strings[i][j*2+0] << "   -   Read gaps: " 
          << fragment_ints[i][j*4+0] << "; Read bases: " << fragment_ints[i][j*4+1] << std::endl;
        std::cout << "Pref Segment: " << fragment_strings[i][j*2+1] << "   -   Prefix gaps: " 
          << fragment_ints[i][j*4+2] << "; Prefix bases: " << fragment_ints[i][j*4+3] << std::endl;
    }
    std::cout << "*****************************" << std::endl;
    std::cout << "\n\n" << std::endl;
  }



  Rcpp::List trim_result;

  trim_result = Rcpp::List::create(Rcpp::Named("sread") = alignment);

  return trim_result;
}

