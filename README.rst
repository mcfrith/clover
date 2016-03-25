Clover: Cis-eLement OVERrepresentation
======================================

Introduction
------------

Clover is a program for identifying functional sites in DNA
sequences. If you give it a set of DNA sequences that share a common
function, it will compare them to a library of sequence motifs
(e.g. transcription factor binding patterns), and identify which if
any of the motifs are statistically overrepresented in the sequence
set.

FAQ
---

| Q: What does it mean if a motif gets a negative raw score, but a low
  P-value?
| A: This means that the motif itself is not overrepresented, but a
  similar motif is. Or, equivalently, it means that the motif is
  overrepresented, but your motif pattern is slightly inaccurate.

| Q: What is meant by the error message "Can't get fragments of
  control sequences to match all target sequences"?
| A: This means that your background sequences are too short, or your
  target sequences are too long. Clover tries to find random fragments
  of the background sequences matched by length to the target
  sequences. If any of the target sequences are longer than all the
  background sequences, this is impossible and you will get this
  message. To fix this, either shorten your target sequences or use
  longer background sequences. (This problem can also arise if the
  sequences are fragmented by ambiguous 'N' nucleotides, e.g. from
  repeat masking. If the background sequences are chosen carefully,
  repeat masking is not necessary and not recommended.)

| Q: How can I make Clover run faster?  
| A: You can split your motif library into several smaller libraries,
  and use them in parallel.

Publication
-----------

Martin C Frith, Yutao Fu, Liqun Yu, Jiang-Fan Chen, Ulla Hansen,
Zhiping Weng (2004). Detection of functional DNA motifs via
statistical over-representation. Nucleic Acids Research 32(4):1372-81.

Setup
-----

Using the command line, go into the clover directory and type "make".

Get a Motif Library
-------------------

To use Clover, you will need a library of sequence motifs, which you
can get from e.g. `JASPAR <http://jaspar.genereg.net/>`_.  You can use
jaspar2clover to convert from JASPAR's matrix_only.txt format to
Clover's format.

Required Input
--------------

Clover requires two inputs: a file of DNA sequences in FASTA format,
and a file of sequence motifs in a FASTA-like format described
below. Any non-alphabetic characters in the sequences are ignored, and
any alphabetic characters except A, C, G, T (uppercase or lowercase)
are converted to 'n' and excluded from matching motifs. The motif file
should look like this::

    >TATA
    0 0 0 10
    10 0 0 0
    0 0 0 10
    10 0 0 0
    >E-box
    1 20 1 1
    (etc)

Each motif begins with a title line containing the character '>'
followed by the motif's name. Subsequent lines represent successive
positions of the motif, from 5' to 3', and the columns contain counts
of A, C, G, and T, respectively, observed at each position. These
numbers typically come from an alignment of several binding sites for
a transcription factor.

Raw Scores and P-values
-----------------------

Clover will compare each motif in turn to the sequence set, and
calculate a "raw score" indicating how strongly the motif is present
in the sequence set. Raw scores by themselves are hard to interpret,
so Clover provides options (which we recommend you use) to determine
the statistical significance of the raw scores. Four ways of
determining statistical significance are available. The first involves
providing Clover with one or more files of background DNA
sequences. Each background file should contain sequences in FASTA
format, with total length much greater than the target sequence
set. For each background set, Clover will repeatedly extract random
fragments matched by length to the target sequences, and calculate raw
scores for these fragments. The proportion of times that the raw score
of a fragment set exceeds or equals the raw score of the target set,
e.g. 0.02, is called a P-value. The P-value indicates the probability
that the motif's presence in the target set can be explained just by
chance. For each motif, a separate P-value is calculated for each
background file.

The second way of determining statistical significance is to
repeatedly shuffle the letters within each target sequence, and use
these shuffled sequence sets as controls. P-values are calculated as
above. The third way is to create random sequences with the same
dinucleotide compositions as each target sequence. The fourth way is
to shuffle the motif matrices, and obtain control raw scores by
comparing the shuffled motifs to the target sequences. When shuffling
a motif, the counts of A, C, G and T within each position are not
shuffled, but the positions are shuffled among one another.

Advice
------

In our experience to date, the use of background sequence sets works
best. However, it is necessary to choose the background sets
carefully: they should ideally come from the same taxonomic group as
the target sequences, and have similar repetitive element and GC
content. We like to cover our bases by using multiple background sets,
e.g. for human target sequences, we might use a human chromosome, a
set of human CpG islands, and a set of human gene upstream regions as
backgrounds. The methods that randomize nucleotides and dinucleotides
suffer from predicting motifs that lie in Alus and other common
repetitive elements to be significant. You should avoid including
orthologous sequences from closely related species, e.g. human and
mouse, as that will artefactually boost the significance of motifs in
these sequences.

Output
------

Clover prints details for statistically significant motifs (all
P-values <= some threshold, by default 0.01), and then finds instances
of these motifs in the sequences. Motif instances are scored using the
standard log likelihood ratio method::

    score = log[ prob(sequence|motif) / prob(sequence|random) ]

Details are printed for motif instances with score >= some threshold,
by default 6.

Options
-------

In addition to the two required inputs, there several options for
modifying Clover's behavior:

-h  Help: print documentation.

-r  Number of randomized/control raw scores to calculate for
    comparison with each target raw score.

-t  P-value threshold: only print results for motifs whose P-values
    don't exceed this amount.

-u  Score threshold for printing locations of significant motifs. This
    parameter doesn't affect raw score and P-value calculations, just
    which motif instances get printed.

-n  Perform sequence (nucleotide) shuffles.

-d  Perform dinucleotide randomizations.

-m  Perform motif shuffles.

-l  Mask (convert to 'n') any lowercase letters in the target
    sequences (and background sequences, if any). Lowercase letters
    are often used to indicate repetitive elements.

-v  Verbose: print per-sequence scores for significant motifs. When
    calculating a motif's raw score, preliminary scores are first
    obtained for the motif compared to each sequence, and these are
    then combined to form the overall raw score. The -v option causes
    these preliminary scores to be displayed.

-p  Pseudocount to add to each entry of the motif
    matrices. Pseudocounts are a widely used technique, with a
    theoretical underpinning in Bayesian statistics, for estimating
    underlying frequencies from a limited number of counts. If your
    matrices contain probabilities rather than counts, you should
    probably set the pseudocount to zero.

-s  Seed for the random number generator (default = 1).

-z  Which DNA strand(s) to analyze: 1=forward, 2=both (default = 2).

Example usage::

    clover -t 0.05 mymotifs myseqs.fa background1.fa background2.fa

**Good luck finding those motifs!**
