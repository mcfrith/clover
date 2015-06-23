// Class for generating random sequences from a Markov model
/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

// The sequences consist of integers between 0 and alphsize-1, inclusive.
// Oligonucleotides are encoded in a fairly obvious way:
// E.g. the oligonucleotide (1, 3, 2) is encoded as
// 2 + 3 * alph_size + 1 * alph_size^2
// This code gives the index of the oligo_counts array corresponding to each oligo
// oligo_counts must contain a^(o+1) elements

#ifndef MCF_MARKOV_H
#define MCF_MARKOV_H

#include <vector>

namespace mcf {
  class markov {
  public:
    // constructor
    markov(unsigned o,  // order of Markov model
	   unsigned a,  // alphabet size
	   const std::vector<unsigned> oligo_counts);
    // counts of all a^(o+1) oligos of length o+1

    // generate a sequence of a given length, APPEND it to seq
    void generate(std::vector<unsigned> & seq, unsigned length) const;

  private:
    unsigned order;  // order of the Markov model
    unsigned alphsize;  // size of alphabet
    unsigned num_prefixes;  // alphsize ^ order (used internally)
    std::vector<std::vector<float> > cond_probs;
  };
}

#endif
