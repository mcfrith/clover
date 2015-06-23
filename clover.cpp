/*
Clover: Cis-eLement OVERrepresentation

Examine whether cis-element motifs are significantly overrepresented in
a set of sequences.

Inputs:
1. fasta-format file of sequences
2. file of matrices
3. Optional but recommended file of background sequences
4. number of random comparisons to perform

Outputs:
1. A raw score for each motif
2. P-values of this score based on randomization tests
3. Locations of significant motifs in the sequences
*/

/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>  // accumulate
#include <string>
#include <vector>
#include <cmath>  // log
#include <libgen.h>  // basename

#include "MCFgen.hpp"  // die
#include "MCFbio.hpp"  // count_oligos
#include "markov.hpp"
#include "args.hpp"

// how many digits needed to write the number in base 10?
// what's the coolest way to do this?
unsigned num_digits(unsigned x)
{
  unsigned n = 1;
  while (x /= 10) ++n;
  return n;
}

// DNA-to-number translator that only recognizes uppercase
inline unsigned nolower_translator(char c)
{
  switch(c) {
  case 'A':
    return 0u;
  case 'C':
    return 1u;
  case 'G':
    return 2u;
  case 'T':
    return 3u;
  default:
    return 4u;
  }
}

// randomly shuffle elements of a vector,
// leaving elements equal to 'm' in place
// Fisher-Yates shuffle
template<class T> void shuffle_except(std::vector<T> & v, const T & m)
{
  for (unsigned i = v.size(); i > 0; --i) {
    if (v[i-1] == m) continue;
    unsigned j;
    do {
      j = rand(i);
    } while (v[j] == m);
    std::swap(v[i-1], v[j]);
  }
}

// randomly shuffle the rows of a matrix, using Fisher-Yates shuffle
// this would be easier if I used a vector of vectors
template<class T> void shuffle_matrix_cols(mcf::matrix<T> & m)
{
  for (unsigned r = m.rows(); r > 0; --r) {
    unsigned s = rand(r);
    if (s != r-1)  // swap_ranges requires non-overlapping ranges
      std::swap_ranges(m[r-1], m[r], m[s]);
  }
}

namespace clover {
  using namespace std;
  using mcf::matrix;
  using mcf::markov;
  using args::alphsize;
  typedef unsigned uint;
  typedef vector<matrix<double> > matvec;

  struct result {
    uint motif;  // index of motif
    double raw_score;
    vector<double> pvalues;
    vector<double> seq_scores;  // score for each sequence

    bool operator> (const result & r) const
    { return raw_score > r.raw_score; }
  };

  struct hit {
    uint motif;  // index of motif
    uint strand;  // 0 or 1: palindromes??
    uint location;
    double score;

    hit(uint m, uint st, uint l, double sc)
      : motif(m), strand(st), location(l), score(sc) {}

    bool operator< (const hit & h) const
    { return location < h.location; }
  };

  vector<matrix<float> > ss_motifs;  // single stranded motifs
  vector<matvec> ds_motifs;  // double strand motifs
  vector<string> motif_names;
  vector<vector<uint> > seqs;
  vector<string> seq_names;
  vector<vector<double> > base_probs;  // residue abundances per sequence
  vector<double> dp;  // dynamic programming matrix used in combine_scores
  vector<result> results;  // motif raw scores & pvalues
  vector<vector<hit> > hits;  // motif locations in each sequence

  // (I think uint(-1) == ~0u == UINT_MAX)
  double scan_seq(const vector<uint> & seq, const matvec & motif, const vector<double> & base_probs, uint seqnum = ~0u, uint motnum = ~0u);
  void init_dp(uint num_seqs);
  double combine_scores(const vector<double> & scores);
  void get_hits();
  void get_raw_scores();
  void rand_test(const vector<vector<uint> > & myseqs, const vector<vector<double> > & b_probs, vector<matvec> & motifs, vector<uint> & losses);
  void copy_masks(const vector<uint> & source, vector<uint> & dest);
  void get_base_probs(const vector<uint> & seq, vector<double> & probs);
  double get_gc_composition(const vector<vector<uint> > & seqs);
  void bg_fragment(const vector<vector<uint> > & bg_seqs, vector<uint> & frag, uint len, const vector<uint> & frag_num, uint frag_tot);
  void shuffle_bgseq(const vector<vector<uint> > & bg_seqs);
  void shuffle_mono();
  void shuffle_di();
  void shuffle_mat();
  void get_ss_motifs();
  void fifth_column(const matrix<float> & m1, matrix<double> & m2);
  void get_ds_motifs(const vector<matrix<float> > & ss_mots, vector<matvec> & ds_mots);
  void get_motifs();
  void get_seqs(const string & filename, vector<vector<uint> > & myseqs, vector<string> & names);
  bool is_significant(const result & r);
  void print_hits(uint wn, uint ws);
  void print_per_sequence(uint wn);
  void print_results();
}

// a result is significant if all of its pvalues are significant
bool clover::is_significant(const result & r)
{
  for (uint p = 0; p < r.pvalues.size(); ++p)
    if (r.pvalues[p] > args::pthresh && r.pvalues[p] < 1 - args::pthresh)
      return false;

  return true;
}

// get the score for one motif in one sequence:
// prob(sequence | motif) / prob(sequence | null model)
// default arguments shouldn't be repeated, apparently
double clover::scan_seq(const vector<uint> & seq, const matvec & motif, const vector<double> & b_probs, uint seqnum, uint motnum)
{
  double tot_score = 0;

  // iterate over motif orientations:
  for (uint m = 0; m < motif.size(); ++m) {
    matrix<double> pssm(motif[m]);
    if (pssm.rows() > seq.size())
      continue;

    // incorporate base probs into PSSM:
    for (uint r = 0; r < pssm.rows(); ++r)
      for (uint c = 0; c < alphsize; ++c)
	pssm[r][c] /= b_probs[c];

    // finally, scan the PSSM against the sequence:
    double score = 0;
    const uint posns = seq.size() - pssm.rows() + 1;
    for (vector<uint>::const_iterator n = seq.begin(); n != seq.begin() + posns; ++n) {
      double s = 1;
      for (uint k = 0; k < pssm.rows(); ++k)
	s *= pssm[k][*(n+k)];
      score += s;
      if (seqnum != ~0u && log(s) >= args::hit_thresh)
	hits[seqnum].push_back(hit(motnum, m, n - seq.begin(), s));
    }

    tot_score += score / posns / motif.size();
  }

  return tot_score;
}

// initialize the dynamic programming matrix used in combine_scores
void clover::init_dp(uint num_seqs)
{
  dp.resize(num_seqs+1);
}

// Combine separate sequence scores into an overall raw score.
// METHOD:
// Define dp[i][j] = average of all products of i elements
// from among the first j elements of scores.
// Hence dp[i][S] = average of all products of i elements from scores.
// Return the average of dp[i][S]. (Excluding i=0.)
double clover::combine_scores(const vector<double> & scores)
{
  const uint score_num = scores.size();

  // dp_i_j   means: dp[i][j]
  // dp_i_j1  means: dp[i][j-1]
  // dp_i1_j1 means: dp[i-1][j-1]

  for (uint j = 1; j <= score_num; ++j) {
    double s = scores[j-1];

    // set the boundary condition: dp[0][j-1] = 1
    double dp_i_j1 = 1;

    // set the boundary condition: dp[j][j-1] = 0
    dp[j] = 0;

    // At this point, the dp vector holds the values: dp[1...j][j-1].

    for (uint i = 1; i <= j; ++i) {
      double dp_i1_j1 = dp_i_j1;
      /****/ dp_i_j1  = dp[i];
      double dp_i_j   = (i * s * dp_i1_j1 + (j-i) * dp_i_j1) / j;
      dp[i]           = dp_i_j;
    }

    // At this point, the dp vector holds the values: dp[1...j][j].
  }

  // At this point, the dp vector holds the values: dp[1...S][S].

  double comb_score = 0;
  for (uint i = 1; i <= score_num; ++i)
    comb_score += dp[i];

  return comb_score / score_num;
}

void clover::get_hits()
{
  cerr << "Getting motif locations..." << endl;
  hits.resize(seqs.size());

  for (uint m = 0; m < ds_motifs.size(); ++m)
    if (is_significant(results[m]))  // assume results finished & not shuffled!
      for (uint s = 0; s < seqs.size(); ++s)
	scan_seq(seqs[s], ds_motifs[m], base_probs[s], s, m);
}

void clover::get_raw_scores()
{
  cerr << "Getting raw scores..." << endl;

  for (uint m = 0; m < ds_motifs.size(); ++m) {
    for (uint s = 0; s < seqs.size(); ++s)
      results[m].seq_scores.push_back(scan_seq(seqs[s], ds_motifs[m], base_probs[s]));
    results[m].raw_score = combine_scores(results[m].seq_scores);
  }
}

// conduct one randomization test on each motif
void clover::rand_test(const vector<vector<uint> > & myseqs, const vector<vector<double> > & b_probs, vector<matvec> & motifs, vector<uint> & losses)
{
  for (uint m = 0; m < motifs.size(); ++m) {
    vector<double> scores;
    for (uint s = 0; s < myseqs.size(); ++s)
      scores.push_back(scan_seq(myseqs[s], motifs[m], b_probs[s]));
    double raw = combine_scores(scores);
    if (raw >= results[m].raw_score)
      ++losses[m];
  }
}

// copy mask characters from one sequence to another
void clover::copy_masks(const vector<uint> & source, vector<uint> & dest)
{
  for (uint i = 0; i < source.size(); ++i)
    if (source[i] == alphsize)
      dest[i] = alphsize;
}

// get background abundances of residues in the sequence
void clover::get_base_probs(const vector<uint> & seq, vector<double> & probs)
{
  probs.clear();  // just in case
  vector<size_t> counts;  // no Laplace's rule
  mcf::count_residues(seq, counts, alphsize);
  const size_t zero = 0;
  const double tot = accumulate(counts.begin(), counts.end(), zero);
  for (uint i = 0; i < alphsize; ++i)
    probs.push_back(counts[i] / tot);
}

void clover::shuffle_mono()
{
  cerr << "Randomizing mononucleotides..." << endl;

  vector<uint> losses(ds_motifs.size(), 0u);
  vector<vector<uint> > r_seqs(seqs);

  for (uint r = 0; r < args::shuffles; ++r) {
    for (uint s = 0; s < r_seqs.size(); ++s)
      shuffle_except(r_seqs[s], alphsize);
    rand_test(r_seqs, base_probs, ds_motifs, losses);
  }

  for (uint m = 0; m < ds_motifs.size(); ++m)
    results[m].pvalues.push_back(losses[m] / (double)args::shuffles);
}

void clover::shuffle_di()
{
  // 1) Preprocessing - create Markov models
  vector<markov> models;
  for (uint s = 0; s < seqs.size(); ++s) {
    vector<uint> dinuc_counts(16, 0);  // no Laplace's rule here
    mcf::count_oligos(seqs[s], dinuc_counts, 2, alphsize);
    models.push_back(markov(1, alphsize, dinuc_counts));
  }

  cerr << "Randomizing dinucleotides..." << endl;

  vector<uint> losses(ds_motifs.size(), 0u);

  // 2) Repeatedly generate Markov sequences and score them:
  for (uint r = 0; r < args::shuffles; ++r) {
    vector<vector<uint> > r_seqs(seqs.size());
    vector<vector<double> > b_probs(seqs.size());

    for (uint s = 0; s < seqs.size(); ++s) {
      models[s].generate(r_seqs[s], seqs[s].size());
      copy_masks(seqs[s], r_seqs[s]);
      get_base_probs(r_seqs[s], b_probs[s]);
    }

    rand_test(r_seqs, b_probs, ds_motifs, losses);
  }

  for (uint m = 0; m < ds_motifs.size(); ++m)
    results[m].pvalues.push_back(losses[m] / (double)args::shuffles);
}

void clover::shuffle_mat()
{
  cerr << "Randomizing matrix columns..." << endl;

  vector<uint> losses(ds_motifs.size(), 0u);
  vector<matrix<float> > r_ss_motifs(ss_motifs);

  for (uint r = 0; r < args::shuffles; ++r) {
    for (uint m = 0; m < r_ss_motifs.size(); ++m)
      shuffle_matrix_cols(r_ss_motifs[m]);
    vector<matvec> r_ds_motifs;
    get_ds_motifs(r_ss_motifs, r_ds_motifs);
    rand_test(seqs, base_probs, r_ds_motifs, losses);
  }

  for (uint m = 0; m < ds_motifs.size(); ++m)
    results[m].pvalues.push_back(losses[m] / (double)args::shuffles);
}

// get a fragment of bg sequences with same length as fg sequence,
// avoiding masked bases
void clover::bg_fragment(const vector<vector<uint> > & bg_seqs, vector<uint> & frag, uint len, const vector<uint> & frag_num, uint frag_tot)
{
  uint b;  // which background sequence
  uint r = rand(frag_tot);
  for (b = 0; b < bg_seqs.size(); ++b) {
    if (frag_num[b] > r)
      break;
    r -= frag_num[b];
  }
  assert(b != bg_seqs.size());

  vector<uint>::const_iterator p;
  const uint posns = bg_seqs[b].size() - len + 1;
  do {
    p = bg_seqs[b].begin() + rand(posns);
  } while (find(p, p + len, alphsize) != p + len);  // check for masked bases

  frag.assign(p, p + len);
}

void clover::shuffle_bgseq(const vector<vector<uint> > & bg_seqs)
{
  cerr << "Analyzing background sequences..." << endl;

  // 1) Preprocessing - how many ways can each fg seq fit in each bg seq?
  vector<vector<uint> > frag_nums;
  for (vector<vector<uint> >::const_iterator i = seqs.begin();
       i != seqs.end(); ++i) {
    frag_nums.push_back(vector<uint>());
    for (vector<vector<uint> >::const_iterator j = bg_seqs.begin();
	 j != bg_seqs.end(); ++j) {
      uint frags = 0;
      uint r = 0;
      for (vector<uint>::const_iterator n = j->begin(); n != j->end(); ++n)
	if (*n == alphsize)
	  r = 0;
	else if (++r >= i->size())
	  ++frags;
      frag_nums.back().push_back(frags);
    }
  }

  vector<uint> frag_tots;
  for (vector<vector<uint> >::const_iterator f = frag_nums.begin();
       f != frag_nums.end(); ++f) {
    uint tot = accumulate(f->begin(), f->end(), 0u);
    if (tot == 0)
      mcf::die("Can't get fragments of control sequences to match all target sequences.");
    frag_tots.push_back(tot);
  }

  vector<uint> losses(ds_motifs.size(), 0u);

  // 2) Repeatedly pick random fragments of the bg sequence set and score them:
  for (uint r = 0; r < args::shuffles; ++r) {
    vector<vector<uint> > r_seqs(seqs.size());
    vector<vector<double> > b_probs(seqs.size());

    for (uint s = 0; s < seqs.size(); ++s) {
      bg_fragment(bg_seqs, r_seqs[s], seqs[s].size(), frag_nums[s], frag_tots[s]);
      copy_masks(seqs[s], r_seqs[s]);
      get_base_probs(r_seqs[s], b_probs[s]);
    }

    rand_test(r_seqs, b_probs, ds_motifs, losses);
  }

  for (uint m = 0; m < ds_motifs.size(); ++m)
    results[m].pvalues.push_back(losses[m] / (double)args::shuffles);
}

// ***** Functions to read input files: *****

// read matrices from a file and partially preprocess them
// Preprocessing steps:
// 1. normalize rows (counts -> probs)
// (2. multiply first row by 0.5: 50% prob of each orientation)
void clover::get_ss_motifs()
{
  ifstream file(args::mat_file.c_str());
  if (!file) mcf::die("Sorry, couldn't open " + args::mat_file);

  const vector<float> pseudos(alphsize, args::pseudocount);
  matrix<float> matf;
  string title;

  while (mcf::get_simple_pssm(file, matf, title, alphsize)) {
    if (matf.rows() == 0) mcf::die("Empty matrix not allowed: " + title);
    motif_names.push_back(title);
    mcf::normalize_pssm(matf, pseudos);
    //    transform(matf[0], matf[1], matf[0], bind2nd(multiplies<float>(), 0.5f));
    ss_motifs.push_back(matf);
  }

  if (!file.eof())  // catches some but not all errors
    mcf::die("Sorry, couldn't understand the matrix file " + args::mat_file);
}

// add a column to a matrix, with all cells = zero
void clover::fifth_column(const matrix<float> & m1, matrix<double> & m2)
{
  for (uint r = 0; r < m1.rows(); ++r) {
    copy(m1[r], m1[r+1], m2[r]);
    m2[r][alphsize] = 0;
  }
}

// finish preprocessing matrices
// 3. add an extra column of all zeros (for 'n' bases)
// 4. make reverse complements
void clover::get_ds_motifs(const vector<matrix<float> > & ss_mots, vector<matvec> & ds_mots)
{
  for (uint m = 0; m < ss_mots.size(); ++m) {
    matrix<float> matf(ss_mots[m]);
    matrix<double> matd(matf.rows(), alphsize+1);
    ds_mots.push_back(matvec());
    fifth_column(matf, matd);
    ds_mots.back().push_back(matd);
    //continue;
    if (!matf.is_rotate180()) {  // if the motif isn't palindromic
      matf.rotate180();  // reverse complement
      fifth_column(matf, matd);
      ds_mots.back().push_back(matd);
    }
  }
}

// read motifs from a file
void clover::get_motifs()
{
  get_ss_motifs();
  get_ds_motifs(ss_motifs, ds_motifs);
}

void clover::get_seqs(const string & filename, vector<vector<uint> > & myseqs, vector<string> & names)
{
  ifstream file(filename.c_str());
  if (!file) mcf::die("Sorry, couldn't open file " + filename);
  unsigned (*translator)(char) = (args::mask_lower ? nolower_translator : mcf::DNA_to_number);

  while (1 == 1) {
    myseqs.push_back(vector<uint>(0));  // push back an empty vector
    string n;
    if (!mcf::get_fasta(file, myseqs.back(), n, translator))
      break;
    names.push_back(n);
  }

  myseqs.pop_back();
}

// ***** Functions to print output: *****

// print locations of motifs in sequences
void clover::print_hits(uint wn, uint ws)
{
  assert(seqs.size() == hits.size());

  string::size_type max_seq_len = 0;
  for (uint s = 0; s < seqs.size(); ++s)
    max_seq_len = max(max_seq_len, seqs[s].size());
  // coordinate field width: at least 3
  const uint wc = max(num_digits(max_seq_len), 3u);  // The u is needed

  cout << "*** Motif Instances with Score >= " << args::hit_thresh << ":\n\n"
       << setw(wn) << "Motif"
       << "  " << setw(2*wc+3) << "Location"
       << "  " << "Strand"
       << "  " << setw(ws) << "Sequence"
       << "  " << "Score" << "\n\n";

  for (uint s = 0; s < seqs.size(); ++s) {
    if (hits[s].empty())
      continue;
    cout << '>' << seq_names[s] << endl;
    sort(hits[s].begin(), hits[s].end());
    for (vector<hit>::const_iterator h = hits[s].begin(); h != hits[s].end(); ++h) {
      const uint motif_width = ss_motifs[h->motif].rows();
      string site;
      for (unsigned int k = 0; k < motif_width; ++k)
	site += mcf::number_to_DNA(seqs[s][h->location + k]);

      cout << setw(wn) << motif_names[h->motif].c_str();
      cout.setf(ios::right);  // right justification
      cout << "  " << setw(wc) << h->location + 1
	   << " - " << setw(wc) << h->location + motif_width;
      cout.setf(ios::left);  // left justification
      if (ds_motifs[h->motif].size() == 1)
	cout << "  " << "  p   ";
      else
	cout << "  " << (h->strand == 0 ? "  +   " : "  -   ");
      cout << "  " << setw(ws) << site.c_str()
	   << "  " << log(h->score) << '\n';
    }
    cout << endl;
  }
}

// print per sequence motif scores
void clover::print_per_sequence(uint wn)
{
  cout << "*** Per Sequence Scores of Significant Motifs:\n\n";
  for (uint s = 0; s < seqs.size(); ++s) {
    cout << '>' << seq_names[s] << endl;
    for (vector<result>::const_iterator r = results.begin(); r != results.end(); ++r)
      if (is_significant(*r))
	cout << setw(wn) << motif_names[r->motif].c_str()
	     << "  " << log(r->seq_scores[s]) << '\n';
    cout << endl;
  }
}

void clover::print_results()
{
  const uint wd = 10;  // data field width
  // name field width: (won't be pretty if names have tabs)
  string::size_type wn = string("Motif").size();
  uint ws = 8;  // site field width: at least 8
  for (uint m = 0; m < ss_motifs.size(); ++m)
    if (is_significant(results[m])) {
      wn = max(wn, motif_names[m].size());
      ws = max(ws, ss_motifs[m].rows());
    }

  for (uint m = 0; m < results.size(); ++m)
    results[m].motif = m;
  sort(results.begin(), results.end(), greater<result>());

  cout.setf(ios::left);  // left justification

  cout << "*** Over- and Under-represented Motifs:\n\n";
  cout << setw(wn) << "Motif" << "  " << setw(wd) << "Raw score";
  if (args::shuf_nuc || args::shuf_di || args::shuf_mat || !args::bg_files.empty())
    cout << " " << "P-value from randomizing";
  cout << '\n' << setw(wn) << "" << "  " << setw(wd) << "";
  if (args::shuf_nuc)
    cout << " " << setw(wd) << "Mononucs";
  if (args::shuf_di)
    cout << " " << setw(wd) << "Dinucs";
  if (args::shuf_mat)
    cout << " " << setw(wd) << "Matrix";
  for (uint i = 0; i < args::bg_files.size(); ++i) {
    // is there no better way to feed a C++ string into basename?
    vector<char> v(args::bg_files[i].begin(), args::bg_files[i].end());
    v.push_back(0);
    string name(basename(&v[0]));
    if (name.size() > wd) name.resize(wd);
    cout << " " << setw(wd) << name.c_str();
  }
  cout << '\n';

  for (vector<result>::const_iterator r = results.begin(); r != results.end(); ++r) {
    if (!is_significant(*r))
      continue;
    cout << setw(wn) << motif_names[r->motif].c_str()
	 << "  " << setw(wd) << log(r->raw_score);
    for (uint p = 0; p < r->pvalues.size(); ++p)
      cout << " " << setw(wd) << r->pvalues[p];
    cout << '\n';
  }
  cout << endl;

  print_hits(wn, ws);

  if (args::verbose)
    print_per_sequence(wn);
}

int main(int argc, char **argv)
{
  using namespace std;

  cout << "Clover: Cis-eLement OVERrepresentation\n"
       << "Compiled on " __DATE__ "\n" << endl;

  args::parse(argc, argv);  // parse command line

  clover::get_motifs();
  clover::get_seqs(args::seq_file, clover::seqs, clover::seq_names);
  if (clover::seqs.empty())
    mcf::die("No sequences read.");
  args::seq_set_info seq_info(clover::seqs);

  clover::base_probs.resize(clover::seqs.size());
  for (unsigned s = 0; s < clover::seqs.size(); ++s)
    clover::get_base_probs(clover::seqs[s], clover::base_probs[s]);

  clover::init_dp(clover::seqs.size());
  clover::results.resize(clover::ds_motifs.size());
  srand(args::random_seed);

  clover::get_raw_scores();

  if (args::shuf_nuc)
    clover::shuffle_mono();
  if (args::shuf_di)
    clover::shuffle_di();
  if (args::shuf_mat)
    clover::shuffle_mat();

  vector<args::seq_set_info> bg_info;
  for (unsigned i = 0; i < args::bg_files.size(); ++i) {
    vector<vector<unsigned> > bg_seqs;
    vector<string> junk;
    cerr << "Reading background sequences..." << endl;
    clover::get_seqs(args::bg_files[i], bg_seqs, junk);
    bg_info.push_back(args::seq_set_info(bg_seqs));
    clover::shuffle_bgseq(bg_seqs);
  }

  clover::get_hits();

  cerr << endl;
  cout.precision(3);
  cout << "Sequence file: " << args::seq_file << " (" << seq_info.num << " sequences, "
       << seq_info.len << " bp, " << seq_info.gc * 100 << "% C+G)\n"
       << "Motif file: " << args::mat_file << " (" << clover::ds_motifs.size() << " motifs)\n"
       << endl;
  clover::print_results();
  args::print(cout, clover::ds_motifs.size(), seq_info, bg_info);
}
