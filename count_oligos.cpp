/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#include <cassert>
#include "MCFgen.hpp"  // int_pow
#include "MCFbio.hpp"

void mcf::count_oligos(
const std::vector<unsigned> & seq, std::vector<unsigned> & counts, unsigned oli_len, unsigned alphsize)
{
  const unsigned oli_num = int_pow(alphsize, oli_len);
  assert(counts.size() == oli_num);
  unsigned oli = 0u;  // code of current oligo
  unsigned len = 0u;  // length of current oligo

  for (std::vector<unsigned>::const_iterator n = seq.begin(); n < seq.end(); ++n) {
    unsigned i = *n;
    if (i < alphsize) {  // is it a normal residue?
      oli = (oli * alphsize) % oli_num + i;
      ++len;
    } else
      len = 0u;  // reset oligo length
    if (len >= oli_len)  // have we reached desired oligo length?
      ++counts[oli];  // count this oligo
  }
}
