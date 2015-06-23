/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#include <algorithm>  // min
#include <numeric>  // accumulate
#include <cassert>
#include "MCFgen.hpp"  // int_pow
#include "markov.hpp"

// constructor
// assumes oligo_counts has a^(o+1) elements!!
mcf::markov::markov(unsigned o, unsigned a, const std::vector<unsigned> oligo_counts)
  : order(o), alphsize(a), num_prefixes(int_pow(a, o)), cond_probs(o+1)
{
  assert(oligo_counts.size() == num_prefixes * alphsize);

  // copy oligo_counts into cond_probs.back():
  cond_probs.back().insert(cond_probs.back().begin(), oligo_counts.begin(), oligo_counts.end());

  for (std::vector<std::vector<float> >::reverse_iterator k = cond_probs.rbegin(); k < cond_probs.rend(); ++k) {
    if (k+1 != cond_probs.rend())
      reserve_or_die(*(k+1), k->size() / alphsize);
    for (std::vector<float>::iterator p = k->begin(); p < k->end(); p += alphsize) {
      float total = std::accumulate(p, p + alphsize, 0.0f);
      /* warning: the probabilities might not sum exactly to 1!
         this might cause problems elsewhere */
      for (std::vector<float>::iterator i = p; i < p + alphsize; ++i) {
        *i /= total;
      }
      if (k+1 != cond_probs.rend())
        (k+1)->push_back(total);
    }
  }
}

// generate a sequence and APPEND it to seq
void mcf::markov::generate(std::vector<unsigned> & seq, unsigned length) const
{
  reserve_or_die(seq, seq.size() + length);

  unsigned end = std::min(order, length);
  unsigned prefix = 0;
  unsigned i;

  for (i = 0u; i < end; ++i) {
    std::vector<float>::const_iterator p = cond_probs[i].begin() + prefix;
    int base = random_choice(p, p + alphsize, 1.0f) - p;
    seq.push_back(base);
    prefix = (prefix + base) * alphsize;
  }

  for (; i < length; ++i) {
    std::vector<float>::const_iterator p = cond_probs[order].begin() + prefix;
    int base = random_choice(p, p + alphsize, 1.0f) - p;
    seq.push_back(base);
    prefix = (prefix + base) % num_prefixes * alphsize;
  }
}
