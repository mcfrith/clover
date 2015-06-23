/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#include <ostream>
#include <unistd.h>  // getopt
#include "MCFgen.hpp"  // die, tostring
#include "MCFbio.hpp"  // count_residues
#include "args.hpp"

namespace args {
  string seq_file;
  string mat_file;
  vector<string> bg_files;
  uint shuffles = 1000;
  double pthresh = 0.01;
  double hit_thresh = 6;
  bool shuf_nuc = false;
  bool shuf_di = false;
  bool shuf_mat = false;
  bool mask_lower = false;
  bool verbose = false;
  double pseudocount = 0.375;
  uint random_seed = 1;  // if srand() isn't called at all, seed=1
}

args::seq_set_info::seq_set_info(const vector<vector<uint> > & seqs)
  : num(seqs.size())
{
  len = 0;
  for (uint s = 0; s < seqs.size(); ++s)
    len += seqs[s].size();

  vector<size_t> counts;  // no Laplace's rule
  for (size_t s = 0; s < seqs.size(); ++s)
    mcf::count_residues(seqs[s], counts, alphsize);
  const size_t zero = 0;
  gc = double(counts[1] + counts[2]) / accumulate(counts.begin(), counts.end(), zero);
}

void args::parse(int argc, char **argv)
{
  using mcf::tostring;

  const string doc = "\
INTRODUCTION\n\
\n\
Clover is a program for identifying functional sites in DNA sequences. If you\n\
give it a set of DNA sequences that share a common function, it will compare\n\
them to a library of sequence motifs (e.g. transcription factor binding\n\
patterns), and identify which if any of the motifs are statistically\n\
overrepresented in the sequence set.\n\
\n\
REQUIRED INPUT\n\
\n\
Clover requires two inputs: a file of DNA sequences in FASTA format, and a file\n\
of sequence motifs in a FASTA-like format described below. Any non-alphabetic\n\
characters in the sequences are ignored, and any alphabetic characters except A,\n\
C, G, T (uppercase or lowercase) are converted to 'n' and excluded from matching\n\
motifs. The motif file should look like this:\n\
\n\
>TATA\n\
0 0 0 10\n\
10 0 0 0\n\
0 0 0 10\n\
10 0 0 0\n\
>E-box\n\
1 20 1 1\n\
(etc)\n\
\n\
Each motif begins with a title line containing the character '>' followed by the\n\
motif's name. Subsequent lines represent successive positions of the motif, from\n\
5' to 3', and the columns contain counts of A, C, G, and T, respectively,\n\
observed at each position. These numbers typically come from from an alignment\n\
of several binding sites for a transcription factor.\n\
\n\
RAW SCORES AND P-VALUES\n\
\n\
Clover will compare each motif in turn to the sequence set, and calculate a \"raw\n\
score\" indicating how strongly the motif is present in the sequence set. Raw\n\
scores by themselves are hard to interpret, so Clover provides options (which we\n\
recommend you use) to determine the statistical significance of the raw scores.\n\
Four ways of determining statistical significance are available. The first\n\
involves providing Clover with one or more files of background DNA sequences.\n\
Each background file should contain sequences in FASTA format, with total length\n\
much greater than the target sequence set. For each background set, Clover will\n\
repeatedly extract random fragments matched by length to the target sequences,\n\
and calculate raw scores for these fragments. The proportion of times that the\n\
raw score of a fragment set exceeds or equals the raw score of the target set,\n\
e.g. 0.02, is called a P-value. The P-value indicates the probability that the\n\
motif's presence in the target set can be explained just by chance. For each\n\
motif, a separate P-value is calculated for each background file.\n\
\n\
The second way of determining statistical significance is to repeatedly shuffle\n\
the letters within each target sequence, and use these shuffled sequence sets as\n\
controls. P-values are calculated as above. The third way is to create random\n\
sequences with the same dinucleotide compositions as each target sequence. The\n\
fourth way is to shuffle the motif matrices, and obtain control raw scores by\n\
comparing the shuffled motifs to the target sequences. When shuffling a motif,\n\
the counts of A, C, G and T within each position are not shuffled, but the\n\
positions are shuffled among one another.\n\
\n\
ADVICE\n\
\n\
In our experience to date, the use of background sequence sets works best.\n\
However, it is necessary to choose the background sets carefully: they should\n\
ideally come from the same taxonomic group as the target sequences, and have\n\
similar repetitive element and GC content. We like to cover our bases by using\n\
multiple background sets, e.g. for human target sequences, we might use a human\n\
chromosome, a set of human CpG islands, and a set of human gene upstream regions\n\
as backgrounds. The methods that randomize nucleotides and dinucleotides suffer\n\
from predicting motifs that lie in Alus and other common repetitive elements to\n\
be significant. You should avoid including orthologous sequences from closely\n\
related species, e.g. human and mouse, as that will artefactually boost the\n\
significance of motifs in these sequences.\n\
\n\
OUTPUT\n\
\n\
Clover prints details for statistically significant motifs (all P-values <= some\n\
threshold, by default 0.01), and then finds instances of these motifs in the\n\
sequences. Motif instances are scored using the standard log likelihood ratio\n\
method:\n\
  score = log[ prob(sequence|motif) / prob(sequence|random) ]\n\
Details are printed for motif instances with score >= some threshold, by default\n\
6.\n\
\n\
OPTIONS\n\
\n\
In addition to the two required inputs, there several options for modifying\n\
Clover's behavior:\n\
 -h  Help: print documentation. You already know this one.\n\
 -r  Number of randomized/control raw scores to calculate for comparison with\n\
     each target raw score.\n\
 -t  P-value threshold: only print results for motifs whose P-values don't\n\
     exceed this amount.\n\
 -u  Score threshold for printing locations of significant motifs. This\n\
     parameter doesn't affect raw score and P-value calculations, just which\n\
     motif instances get printed.\n\
 -n  Perform sequence (nucleotide) shuffles.\n\
 -d  Perform dinucleotide randomizations.\n\
 -m  Perform motif shuffles.\n\
 -l  Mask (convert to 'n') any lowercase letters in the target and sequences\n\
     (and background sequences, if any). Lowercase letters are often used to\n\
     indicate repetitive elements.\n\
 -v  Verbose: print per-sequence scores for significant motifs. When calculating\n\
     a motif's raw score, preliminary scores are first obtained for the motif\n\
     compared to each sequence, and these are then combined to form the overall\n\
     raw score. The -v option causes these preliminary scores to be displayed.\n\
 -p  Pseudocount to add to each entry of the motif matrices. Pseudocounts are a\n\
     widely used technique, with a theoretical underpinning in Bayesian\n\
     statistics, for estimating underlying frequencies from a limited number of\n\
     counts. If your matrices contain probabilities rather than counts, you\n\
     should probably set the pseudocount to zero.\n\
 -s  Seed for random number generator (default = 1).\n\
\n\
Example usage:\n\
  clover -t 0.05 mymotifs myseqs.fa background1.fa background2.fa\n\
\n\
Good luck finding those motifs!\n\
";

  const string usage =
    "Usage summary: clover [options] mymotifs myseqs.fa [BGfiles]\n"
    "Options:\n"
    "-h  help: print documentation\n"
    "-r  number of randomizations (" + tostring(shuffles) + ")\n"
    "-t  P-value threshold (" + tostring(pthresh) + ")\n"
    "-u  motif score threshold (" + tostring(hit_thresh) + ")\n"
    "-n  randomize sequences (nucleotides)\n"
    "-d  randomize sequences (dinucleotides)\n"
    "-m  randomize motifs\n"
    "-l  filter lowercase letters in sequences\n"
    "-v  verbose: print scores of significant motifs for each sequence\n"
    "-p  pseudocount to add to all matrix elements (" + tostring(pseudocount) + ")\n"
    "-s  seed for random number generator (" + tostring(random_seed) + ")\n"
    ;

  int c;

  while ((c = getopt(argc, argv, "hr:t:u:ndmlvp:s:")) != -1)
    switch (c) {
    case 'h':
      cout << doc << endl;
      exit(0);
    case 'r':
      shuffles = atoi(optarg);
      break;
    case 't':
      pthresh = atof(optarg);
      break;
    case 'u':
      hit_thresh = atof(optarg);
      break;
    case 'n':
      shuf_nuc = true;
      break;
    case 'd':
      shuf_di = true;
      break;
    case 'm':
      shuf_mat = true;
      break;
    case 'l':
      mask_lower = true;
      break;
    case 'v':
      verbose = true;
      break;
    case 'p':
      pseudocount = atof(optarg);
      break;
    case 's':
      random_seed = atoi(optarg);
      break;
    case '?':
      mcf::die("\n" + usage);  // "invalid option" message is printed by getopt
    }

  if (optind + 2 > argc)  // there should be 2 more non-option arguments
    mcf::die("Error: motif and sequence files required\n\n" + usage);

  mat_file = argv[optind++];
  seq_file = argv[optind++];

  while (optind < argc)
    bg_files.push_back(argv[optind++]);
}

// print values of command line arguments - always do this!
void args::print(ostream & stream, uint mat_num, const seq_set_info & seq_info, const vector<seq_set_info> & bg_info)
{
  /*  stream
    << "Sequence file: " << seq_file << " (" << seq_info.num << " sequences, "
    << seq_info.len << " bp, " << seq_info.gc * 100 << "% C+G)\n"
    << "Motif file: " << mat_file << " (" << mat_num << " motifs)\n"
    ;*/
  if (!bg_files.empty()) {
    stream << "Background sequence files:\n";
    for (uint s = 0; s < bg_files.size(); ++s)
      stream << "  " << bg_files[s] << " (" << bg_info[s].num << " sequences, "
	     << bg_info[s].len << " bp, " << bg_info[s].gc * 100 << "% C+G)\n";
  }
  stream
    << "Randomizations:         " << shuffles << '\n'
    << "P-value threshold:      " << pthresh << '\n'
    << "Motif score threshold:  " << hit_thresh << '\n'
    << "Randomize nucleotides   " << (shuf_nuc ? "ON\n" : "OFF\n")
    << "Randomize dinucleotides " << (shuf_di ? "ON\n" : "OFF\n")
    << "Randomize motifs        " << (shuf_mat ? "ON\n" : "OFF\n")
    << "Lowercase filtering     " << (mask_lower ? "ON\n" : "OFF\n")
    << "Pseudocount:            " << pseudocount << '\n'
    << "Random seed:            " << random_seed << '\n'
    << flush;
}
