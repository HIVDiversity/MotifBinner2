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

void findPrefixMatch(StringSet<IupacString> haystack,
  StringSet<CharString> id, StringSet<CharString> qual,
  IupacString needle)
{
  typedef int TValue;
  typedef Score<TValue, ScoreMatrix<Iupac, Default> > TScoringScheme;
  int const gapOpenScore = -1;
  int const gapExtendScore = -1;

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
  Align<IupacString, ArrayGaps> align;
  resize(rows(align), 2);
  Row<Align<IupacString, ArrayGaps> >::Type & row0 = row(align, 0);
  Row<Align<IupacString, ArrayGaps> >::Type & row1 = row(align, 1);

  for (unsigned i = 0; i < length(haystack); ++i)
  {
    seq_len = length(haystack[i]);
    assignSource(row(align, 0), infix(haystack[i], 0, std::min(50, seq_len)));
    assignSource(row(align, 1), needle);

    int score = globalAlignment(align, scoringScheme, AlignConfig<true, false, false, true>(), LinearGaps());
    std::cout << i << " = " << score << ";  " ;
  }
  std::cout << std::endl;
  assignSource(row(align, 0), needle);
  assignSource(row(align, 1), needle);

  int score = globalAlignment(align, Score<int, Simple>(0, -1, -1), AlignConfig<true, false, false, true>(), LinearGaps());
  std::cout << "best alignment: prefix to prefix" << std::endl;
  std::cout << score << std::endl;
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

