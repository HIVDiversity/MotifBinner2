#include <Rcpp.h>
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <seqan/align.h>
#include <vector>

using namespace seqan;
using namespace Rcpp;

// Extend SeqAn by a user-define scoring matrix.
namespace seqan {

// We have to create a new specialization of the ScoringMatrix_ class
// for the DNA alphabet.  For this, we first create a new tag.
struct UserDefinedMatrix {};
// Then, we specialize the class ScoringMatrix_ for the Iupac alphabet.
template <>
struct ScoringMatrixData_<int, Iupac, UserDefinedMatrix>
{
  enum
  {
    VALUE_SIZE = ValueSize<Iupac>::VALUE,
    TAB_SIZE = VALUE_SIZE * VALUE_SIZE
  };

  static inline int const * getData()
  {
    // The user defined data table.
    static int const _data[TAB_SIZE] =
    {
       0, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
      -1,  0, -1,  0, -1,  0, -1,  0, -1,  0, -1,  0, -1,  0, -1,  0,
      -1, -1,  0,  0, -1, -1,  0,  0, -1, -1,  0,  0, -1, -1,  0,  0,
      -1,  0,  0,  0, -1,  0,  0,  0, -1,  0,  0,  0, -1,  0,  0,  0,
      -1, -1, -1, -1,  0,  0,  0,  0, -1, -1, -1, -1,  0,  0,  0,  0,
      -1,  0, -1,  0,  0,  0,  0,  0, -1,  0, -1,  0,  0,  0,  0,  0,
      -1, -1,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,
      -1,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,
       0, -1, -1, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0, -1,  0, -1,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0, -1,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0, -1, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0
    };
    return _data;
  }
};
}  // namespace seqan

template <typename TString>
std::string convert_to_string(TString x, char filler, int num_fill)
{
  std::string y;
  resize(y, length(x)+num_fill+1);
  if (num_fill > 0){
    for (int i = 0; i != num_fill; ++i)
    {
      y[i] = filler;
    }
  }
  for (int i = 0; i != length(x); ++i)
  {
    y[i+num_fill] = x[i];
  }
  return y;
}

//template <typename TGap>
//std::string convert_to_string_gap(TGap x)
//{
//  std::string z;                                                                                    
//  resize(z, length(x));
//  typedef Iterator<TGap> TGapIterator;                                                        
//  TGapIterator it = begin(x), itEnd = end(x);                                                 
//  int i = 0;                                                                                        
//  for(; it != itEnd; ++it)                                                                          
//  {                                                                                                 
//    Iupac c = isGap(it) ? gapValue<Iupac>() : *it;                                                  
//    z[i] = c;                                                                                       
//    i++;                                                                                            
//  }                                                                                                 
//  return z;
//}

Rcpp::List findPrefixMatch(StringSet<IupacString> haystack,
  StringSet<CharString> id, StringSet<CharString> qual,
  IupacString needle, int pref_len, int pid_len, int suf_len)
{
  typedef int TValue;
  typedef Score<TValue, ScoreMatrix<Iupac, Default> > TScoringScheme;
  int const gapOpenScore = -1;
  int const gapExtendScore = -1;
  //StringSet<IupacString> trim_haystack(length(haystack));
  //StringSet<CharString> trim_qual(length(qual));
  //
  std::vector<std::string> pref(length(haystack));
  std::vector<std::string> pid(length(haystack));
  std::vector<std::string> suf(length(haystack));
  std::vector<std::string> read(length(haystack));
  std::vector<std::string> read_qual(length(qual));
  std::vector<std::string> trim_id(length(id));

  TScoringScheme scoringScheme(gapExtendScore, gapOpenScore);

  for (unsigned i = 0; i < ValueSize<Iupac>::VALUE; ++i)
  {
    for (unsigned j = 0; j < ValueSize<Iupac>::VALUE; ++j)
    {
      setScore(scoringScheme, Iupac(i), Iupac(j), i * j);
    }
  }
  setDefaultScoreMatrix(scoringScheme, UserDefinedMatrix());

  // iteration trackers
  int seq_len = 0;
  int gap_search_range;
  int haystack_gaps;
  int needle_gaps;
    
  int aln_pref_start;
  int aln_pid_start;
  int aln_suf_start;
  int aln_read_start;
  int gaps_before_read;

  // result storage
  std::vector<int> scores(length(haystack));
  std::vector<int> trim_spots(length(haystack));
  std::vector<int> gaps_at_front_of_read(length(haystack));
    
  // alginment data structures
  Align<IupacString, ArrayGaps> align;
  resize(rows(align), 2);
  Row<Align<IupacString, ArrayGaps> >::Type & row0 = row(align, 0);
  Row<Align<IupacString, ArrayGaps> >::Type & row1 = row(align, 1);

  for (unsigned i = 0; i < length(haystack); ++i)
  {
    seq_len = length(haystack[i]);
    assignSource(row(align, 0), infix(haystack[i], 0, std::min(80, seq_len)));
    assignSource(row(align, 1), needle);

    int score = globalAlignment(align, scoringScheme, AlignConfig<true, false, false, true>(), LinearGaps());
    scores[i] = score;

    haystack_gaps = 0;
    needle_gaps = 0;
    aln_pref_start = -1;
    aln_pid_start = -1;
    aln_suf_start = -1;
    aln_read_start = -1;
    gaps_before_read = -1;

    for (unsigned j = 0; j != length(row1); ++j)
    {
      if (row0[j] == '-')
      { haystack_gaps++; }
      if (row1[j] == '-')
      { needle_gaps++; }

      if (pid_len > 0)
      {
        if (aln_pref_start == -1 and (j - needle_gaps) == 1)
        { aln_pref_start = j; }
        if (aln_pid_start == -1 and (j - needle_gaps) == pref_len)
        { aln_pid_start = j; }
        if (aln_suf_start == -1 and (j - needle_gaps) == pref_len + pid_len)
        { aln_suf_start = j; }
        if (aln_read_start == -1 and (j - needle_gaps) == pref_len + pid_len + suf_len)
        { aln_read_start = j; 
          gaps_before_read = haystack_gaps; }
      } else {
        if (aln_pref_start == -1 and (j - needle_gaps) == 1)
        { aln_pref_start = j; 
          gaps_before_read = haystack_gaps; }
      }
    }

//    // old trim_spot calc - keep for now to test against.
//    for (unsigned j = 0; j != length(row1); ++j)
//    {
//      if (row1[j] != '-')
//      {
//        trim_spots[i] = j;
//        break;
//      }
//    }

//    // front gaps - do not incorporate with other tracking
//    //   it gets too complex and leads to many unnecessary
//    //   if evaluations
//    gaps_at_front_of_read[i] = 0;
//    gap_search_range = std::min(length(needle)+trim_spots[i], length(row0));
//    for (unsigned j = trim_spots[i]; j <= gap_search_range; ++j)
//    {
//      if (row0[j] == '-')
//      {
//        gaps_at_front_of_read[i] = j - trim_spots[i] + 1;
//      } else {
//        break;
//      }
//    }

    std::cout << "pref" << std::endl;
    pref[i] = convert_to_string(infix(row0, 5, 10), '-', 0);
//    pref[i] = convert_to_string(infix(row0, aln_pref_start, aln_pid_start), '-', 0);
//    if (pid_len > 0)
//    {
//    std::cout << "pid" << std::endl;
//      pid[i] = convert_to_string(infix(row0, aln_pid_start, aln_suf_start), '-', 0);
//    }
//    if (suf_len > 0)
//    {
//    std::cout << "suf" << std::endl;
//      suf[i] = convert_to_string(infix(row0, aln_suf_start, aln_read_start), '-', 0);
//    }
    std::cout << "read" << std::endl;
    read[i] = convert_to_string(infix(haystack[i], aln_read_start - gaps_before_read, 
          length(haystack[i])-1), '-', 0);
    read_qual[i] = convert_to_string(infix(qual[i], aln_read_start - gaps_before_read, 
          length(qual[i])-1), '!', 0);
    trim_id[i] = toCString(id[i]);  //convert_to_string(id[i], '-', 0);


    std::cout << "Aligned Prefix " << pref[i] << std::endl;
    std::cout << "Sequence " << i << std::endl;
    std::cout << align;
    std::cout << "Score = " << scores[i] << ";  Trim Spot = " << trim_spots[i] << "; Last Read Gap = " << gaps_at_front_of_read[i] << std::endl;
    std::cout << " --------------------- " << std::endl;
    std::cout << std::endl;

  }
  std::cout << std::endl;
  std::cout << "number of alignments performed " << length(haystack) << std::endl;
  std::cout << "number of alignments results obtained " << length(pref) << std::endl;

//  // compute best possible score
//  assignSource(row(align, 0), needle);
//  assignSource(row(align, 1), needle);
//  int score = globalAlignment(align, Score<int, Simple>(0, -1, -1), AlignConfig<true, false, false, true>(), LinearGaps());
//  std::cout << "best alignment: prefix to prefix" << std::endl;
//  std::cout << score << std::endl;

  return Rcpp::List::create(Rcpp::Named("sread") = read,
                            Rcpp::Named("id") = trim_id,
                            Rcpp::Named("qual") = read_qual,
                            Rcpp::Named("score") = scores,
                            Rcpp::Named("trim_spot") = trim_spots,
                            Rcpp::Named("gaps_at_front_of_read") = gaps_at_front_of_read);
}

// [[Rcpp::export]]

Rcpp::List trimEnds_cpp(CharacterVector r_sread, CharacterVector r_id, 
    CharacterVector r_qual, CharacterVector r_prefix, 
    int pref_len, int pid_len, int suf_len)
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

  Rcpp::List trim_result;
  trim_result = findPrefixMatch(sq_sread, sq_id, // haystack, id
      sq_qual, sq_prefix[0], // qual, needle
      pref_len, pid_len, suf_len); //pref_len, pid_len, suf_len

  return trim_result;
}
