/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#include "MCFbio.hpp"

void mcf::count_residues(
const std::vector<unsigned> & seq, std::vector<size_t> & counts, unsigned alphsize)
{
  if (counts.size() < alphsize)
    counts.resize(alphsize);

  for (std::vector<unsigned>::const_iterator n = seq.begin(); n < seq.end(); ++n) {
    unsigned i = *n;
    if (i < alphsize)
      ++counts[i];
  }
}
