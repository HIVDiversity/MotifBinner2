#include <Rcpp.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>
#include <time.h>

// [[Rcpp::plugins("cpp11")]]

using namespace Rcpp;

struct dynMatEntry {
  int origin_i;
  int origin_j;
  int Fscore;
};

struct alignResult {
  std::string vseq;
  std::string hseq;
  int score;
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

  // The origins aka the arrows
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

alignResult align(std::string vseq, std::string hseq,
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

  alignResult result;
  result.vseq = vseq_align_fwd;
  result.hseq = hseq_align_fwd;
  result.score = Fmat[backtrace_start.first][backtrace_start.second].Fscore;

  return result;
}

// [[Rcpp::export]]
Rcpp::List trimFront_cpp(CharacterVector r_sread, CharacterVector r_qual,
    CharacterVector r_primer, std::vector<int> prefix_lens) // int pref_len, int pid_len, int suf_len, int verbosity)
{

  alignResult alignment;
  std::string read_start;
  std::string full_read;

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
  std::string read_qual_segment;
  std::string prefix_segment;
  Rcpp::NumericVector scores(r_sread.size());

  std::vector<int> read_front_gaps(r_sread.size());
  std::vector<int> prefix_front_gaps(r_sread.size());
  std::vector<std::string> rest_of_read(r_sread.size());
  std::vector<std::string> rest_of_read_qual(r_sread.size());

  Rcpp::NumericMatrix prefix_fragment_gaps(r_sread.size(), prefix_lens.size());
  Rcpp::NumericMatrix prefix_fragment_bases(r_sread.size(), prefix_lens.size());
  Rcpp::NumericMatrix read_fragment_gaps(r_sread.size(), prefix_lens.size());
  Rcpp::NumericMatrix read_fragment_bases(r_sread.size(), prefix_lens.size());

  Rcpp::CharacterMatrix prefix_fragments(r_sread.size(), prefix_lens.size());
  Rcpp::CharacterMatrix read_fragments(r_sread.size(), prefix_lens.size());
  Rcpp::CharacterMatrix read_qual_fragments(r_sread.size(), prefix_lens.size());

  double loop_start = 0;
  double align_time = 0;
  double total_loop_time = 0;

  std::string gapped_qual;
  std::string full_qual;

  for (int i=0; i!= r_sread.size(); i++)
  {
    gapped_qual.clear();
    loop_start = clock();
    full_read = Rcpp::as<std::string>(r_sread[i]);
    full_qual = Rcpp::as<std::string>(r_qual[i]);
    read_start = full_read.substr(0, 20+Rcpp::as<std::string>(r_primer[0]).length());

    alignment = align(read_start,
                      Rcpp::as<std::string>(r_primer),
                      false, false, false, false);
    align_time += (clock() - loop_start);

    if (alignment.vseq.length() != alignment.hseq.length())
    {
      throw std::range_error("Alignments must be of equal lengths");
    }

    std::string& read = alignment.vseq; // reference for easier name
    std::string& prefix = alignment.hseq; // reference for easier name
    scores[i] = alignment.score;
    read_gaps = 0; read_bases = 0;
    prefix_gaps = 0; prefix_bases = 0;
    cur_prefix_fragment = 0;
    next_prefix_break = prefix_lens[0];

    read_front_gaps[i] = -1; prefix_front_gaps[i] = -1;

//    if (i % 100  == 0) {
//      std::cout << i << " ";
//    }
//    std::cout << "Sequence Number: " << i << std::endl;
//    std::cout << "Read Front Gaps: " << read_front_gaps[i] << ";   Prefix Front Gaps: " 
//      << prefix_front_gaps[i] << std::endl;
//    std::cout << read << std::endl;
//    std::cout << prefix << std::endl;

    for (int j = 0; j < read.size(); ++j)
    {
      // basic tracking
      if (read[j] == '-'){ 
        ++read_gaps; 
        gapped_qual += '!';
      } else { 
        gapped_qual += full_qual[read_bases];
        ++read_bases;
      }
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
        read_segment = read.substr(cur_fragment_start_in_aln, j - cur_fragment_start_in_aln + 1);
        read_qual_segment = gapped_qual.substr(cur_fragment_start_in_aln, j - cur_fragment_start_in_aln + 1);
        prefix_segment = prefix.substr(cur_fragment_start_in_aln, j - cur_fragment_start_in_aln + 1);
        cur_fragment_start_in_aln = j+1;

        read_fragment_gaps(i, cur_prefix_fragment) = read_gaps - prev_read_gaps;
        read_fragment_bases(i, cur_prefix_fragment) = read_bases - prev_read_bases;
        prefix_fragment_gaps(i, cur_prefix_fragment) = prefix_gaps - prev_prefix_gaps;
        prefix_fragment_bases(i, cur_prefix_fragment) = prefix_bases - prev_prefix_bases;

        read_fragments(i, cur_prefix_fragment) = read_segment;
        prefix_fragments(i, cur_prefix_fragment) = prefix_segment;
        read_qual_fragments(i, cur_prefix_fragment) = read_qual_segment;

//        std::cout << "Read Segment: " << read_fragments[i][cur_prefix_fragment] << "   -   Read gaps: " 
//          << read_fragment_gaps[i][cur_prefix_fragment] << "; Read bases: " << 
//          read_fragment_bases[i][cur_prefix_fragment] << std::endl;
//        std::cout << "Pref Segment: " << prefix_fragments[i][cur_prefix_fragment] << "   -   Prefix gaps: " 
//          << prefix_fragment_gaps[i][cur_prefix_fragment] << "; Prefix bases: " << 
//          prefix_fragment_bases[i][cur_prefix_fragment] << std::endl;
        
        cur_prefix_fragment ++;
        if (cur_prefix_fragment < prefix_lens.size()){
          next_prefix_break += prefix_lens[cur_prefix_fragment];
        } else {
          rest_of_read[i] = full_read.substr(read_bases, std::string::npos); //full_read.size()-read_bases);
          rest_of_read_qual[i] = full_qual.substr(read_bases, std::string::npos); 
          break; // loop through bases in alignment
        }
        
        prev_read_gaps = read_gaps;
        prev_read_bases = read_bases;
        prev_prefix_gaps = prefix_gaps;
        prev_prefix_bases = prefix_bases;
      } // if transition to next prefix fragment
    } // loop over bases of read
//    std::cout << rest_of_read[i] << std::endl;
//    std::cout << "*****************************" << std::endl;
//    std::cout << "\n\n" << std::endl;
    total_loop_time += (clock() - loop_start);
  }

  std::cout << std::endl;
  std::cout << "Total loop time: " << total_loop_time << std::endl;
  std::cout << "Total align time: " << align_time << std::endl;
  std::cout << "Percentage spent on align: " << align_time/total_loop_time << std::endl;

  Rcpp::List trim_result;

  trim_result = Rcpp::List::create(
    Rcpp::Named("scores") = scores,

    Rcpp::Named("read_front_gaps") = read_front_gaps,
    Rcpp::Named("prefix_front_gaps") = prefix_front_gaps,
    Rcpp::Named("rest_of_read") = rest_of_read,
    Rcpp::Named("rest_of_read_qual") = rest_of_read_qual,

    Rcpp::Named("read_fragments") = read_fragments,
    Rcpp::Named("read_qual_fragments") = read_qual_fragments,
    Rcpp::Named("prefix_fragments") = prefix_fragments,

    Rcpp::Named("read_fragment_gaps") = read_fragment_gaps,
    Rcpp::Named("read_fragment_bases") = read_fragment_bases,
    Rcpp::Named("prefix_fragment_gaps") = prefix_fragment_gaps,
    Rcpp::Named("prefix_fragment_bases") = prefix_fragment_bases
      );

  return trim_result;
}

// [[Rcpp::export]]
Rcpp::List transfer_gaps_cpp(CharacterVector aligned_read, CharacterVector r_qual, NumericVector gap_only_cols)
{
  int read_bases;
  std::string read;
  std::string gapped_qual;
  std::string gapped_read;
  std::string qual;
  Rcpp::CharacterVector all_gapped_qual(aligned_read.size());
  Rcpp::CharacterVector all_gapped_read(aligned_read.size());
  for (int i=0; i!= aligned_read.size(); i++)
  {
    read_bases = 0;
    read = Rcpp::as<std::string>(aligned_read[i]);
    qual = Rcpp::as<std::string>(r_qual[i]);
    gapped_qual.clear();
    gapped_read.clear();
    for (int j = 0; j < read.size(); ++j)
    {
      if (std::find(gap_only_cols.begin(), gap_only_cols.end(), j) != gap_only_cols.end()){
        if (read[j] != '-'){
          throw std::range_error("Non-gap letter in gap only column");
        }
      } else if (read[j] == '-'){ 
        gapped_qual += '!';
        gapped_read += '-';
      } else { 
        gapped_qual += qual[read_bases];
        gapped_read += read[j];
        ++read_bases;
      }
    }
    all_gapped_qual[i] = gapped_qual;
    all_gapped_read[i] = gapped_read;
  }

  Rcpp::List result;

  result = Rcpp::List::create(
    Rcpp::Named("reads") = all_gapped_read,
    Rcpp::Named("quals") = all_gapped_qual
      );
  return result;
}


// [[Rcpp::export]]
Rcpp::List score_alignment_positions(CharacterVector reads, NumericMatrix q_mat)
{
  std::map<char, int> score_mat_row_lookup = 
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

  int k = 0;
  char c_let;
  int c_qual;
  int score_mat_row;
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
