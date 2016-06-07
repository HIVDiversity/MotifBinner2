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

Rcpp::List findPrefixMatch(StringSet<IupacString> haystack,
  StringSet<CharString> id, StringSet<CharString> qual,
  IupacString needle)
{
  typedef int TValue;
  typedef Score<TValue, ScoreMatrix<Iupac, Default> > TScoringScheme;
  int const gapOpenScore = -1;
  int const gapExtendScore = -1;
  //StringSet<IupacString> trim_haystack(length(haystack));
  //StringSet<CharString> trim_qual(length(qual));
  std::vector<std::string> trim_haystack(length(haystack));
  std::vector<std::string> trim_qual(length(qual));
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

  int seq_len = 0;
  int gap_search_range;
  std::vector<int> scores(length(haystack));
  std::vector<int> trim_spots(length(haystack));
  std::vector<int> gaps_at_front_of_read(length(haystack));

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

    for (unsigned j = 0; j != length(row1); ++j)
    {
      if (row1[j] != '-')
      {
        trim_spots[i] = j;
        break;
      }
    }

    gaps_at_front_of_read[i] = 0;
    gap_search_range = std::min(length(needle)+trim_spots[i], length(row0));
    for (unsigned j = trim_spots[i]; j <= gap_search_range; ++j)
    {
      if (row0[j] == '-')
      {
        gaps_at_front_of_read[i] = j - trim_spots[i] + 1;
      } else {
        break;
      }
    }

    trim_haystack[i] = convert_to_string(infix(haystack[i], trim_spots[i], length(haystack[i])-1),
        '-', std::max(0, gaps_at_front_of_read[i]));
    trim_qual[i] = convert_to_string(infix(qual[i], trim_spots[i], length(qual[i])-1),
        '!', std::max(0, gaps_at_front_of_read[i]));
    trim_id[i] = toCString(id[i]);  //convert_to_string(id[i], '-', 0);
//    std::cout << "Sequence " << i << std::endl;
//    std::cout << align;
//    std::cout << "Score = " << scores[i] << ";  Trim Spot = " << trim_spots[i] << "; Last Read Gap = " << gaps_at_front_of_read[i] << std::endl;
//    std::cout << " --------------------- " << std::endl;
//    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "number of alignments performed " << length(haystack) << std::endl;
  std::cout << "number of alignments results obtained " << length(trim_haystack) << std::endl;
  assignSource(row(align, 0), needle);
  assignSource(row(align, 1), needle);

  int score = globalAlignment(align, Score<int, Simple>(0, -1, -1), AlignConfig<true, false, false, true>(), LinearGaps());
//  std::cout << "best alignment: prefix to prefix" << std::endl;
//  std::cout << score << std::endl;

  return Rcpp::List::create(Rcpp::Named("sread") = trim_haystack,
                            Rcpp::Named("id") = trim_id,
                            Rcpp::Named("qual") = trim_qual,
                            Rcpp::Named("score") = scores,
                            Rcpp::Named("trim_spot") = trim_spots,
                            Rcpp::Named("gaps_at_front_of_read") = gaps_at_front_of_read);
}

// [[Rcpp::export]]

Rcpp::List trimEnds_cpp(CharacterVector r_sread, CharacterVector r_id, 
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

  Rcpp::List trim_result;
  trim_result = findPrefixMatch(sq_sread, sq_id, sq_qual, sq_prefix[0]);

  return trim_result;
}

