/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#include <iosfwd>  // ostream
#include <stddef.h>  // size_t
#include <string>
#include <vector>

namespace args {
  using namespace std;
  typedef unsigned uint;

  struct seq_set_info {  // misc info about a sequence set
    size_t num;  // number of seqs
    size_t len;  // total bp
    double gc;   // GC composition

    seq_set_info(const vector<vector<uint> > & seqs);
  };

  const uint alphsize = 4;

  extern string seq_file;  // input file: sequences to be analyzed
  extern string mat_file;  // input file: count matrices to use
  extern vector<string> bg_files;  // background sequence files
  extern uint shuffles;  // number of randomization tests to perform
  extern double pthresh;  // P-value threshold
  extern double hit_thresh;  // score threshold for individual motif hits
  extern bool shuf_nuc;  // randomize mononucleotides?
  extern bool shuf_di;  // randomize dinucleotides?
  extern bool shuf_mat;  // randomize matrix columns?
  extern bool mask_lower;  // mask lowercase letters?
  extern bool verbose;  // print scores of significant motifs for each sequence
  extern double pseudocount;  // add to matrix entries
  extern uint random_seed;

  void parse(int argc, char **argv);  // parse the arguments using getopt
  void print(ostream & strm, uint mat_num, const seq_set_info & seq_info, const vector<seq_set_info> & bg_seq_info);
}
