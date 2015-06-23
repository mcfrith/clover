// Read a position specific weight/count/score matrix in eFASTA format
/* Written by Martin C Frith */
/* I intend that anyone who finds this code useful be free to use,
   modify, or redistribute it without any restrictions. Naturally, I
   hope it will not be used for evil purposes. */

#include <istream>
#include "MCFbio.hpp"

// Industrial strength PSSM reader
// (alphsize could be figured out from the number of entries per line...)
std::istream & mcf::get_simple_pssm(
  std::istream & strm,
  matrix<float> & mat,
  std::string & title,
  unsigned alphsize)  // not allowed to repeat the default argument?
{
  std::string t;
  matrix<float> m(0, alphsize);
  char c = 0;
  bool titflag = false;  // have we read a title line yet?

  while (strm >> c) {
    if (c == '>') {
      if (titflag || m.rows() != 0) {
	strm.unget();
	break;
      } else {
	std::getline(strm, t);
	titflag = true;
      }
    } else if (c == '#') {  // skip comments
      std::string junk;
      std::getline(strm, junk);
    } else {
      strm.unget();
      std::vector<double> v;
      for (unsigned i = 0; i < alphsize; ++i) {
	double d;
	strm >> d;
	v.push_back(d);
      }
      if (!strm)
	return strm;  // failed to read alphsize doubles
      m.push_row(v.begin());
    }
  }

  // if reached EOF but read something, clear the stream state:
  if (strm.eof() && (titflag || m.rows() != 0))
    strm.clear();

  if (strm) {
    title = t;
    mat = m;
  }
  return strm;
}
