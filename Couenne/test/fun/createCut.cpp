#include <OsiCuts.hpp>

#define COUENNE_EPS 1e-7

int createCut (OsiCuts &cs,
	       double rhs, int sign, 
	       int i1, double c1,
	       int i2, double c2,
	       int i3, double c3,
	       bool is_global) {

  int nterms = 0;

  if (i1 >= 0) nterms++; else c1 = 0;
  if (i2 >= 0) nterms++; else c2 = 0;
  if (i3 >= 0) nterms++; else c3 = 0;

  if (!nterms) // nonsense cut
    return 0;

  double    *coeff = new double [nterms]; 
  int       *index = new int    [nterms];
  OsiRowCut *cut   = new OsiRowCut;

  if (i1 >= 0) {coeff [0] = c1; index [0] = i1;}
  if (i2 >= 0) {coeff [1] = c2; index [1] = i2;}
  if (i3 >= 0) {coeff [2] = c3; index [2] = i3;}

  if (sign <= 0) cut -> setUb (rhs);
  if (sign >= 0) cut -> setLb (rhs);

  cut -> setRow (nterms, index, coeff);

  delete [] coeff;
  delete [] index;

  cut -> setGloballyValid (is_global); // global?
  cs.insert (cut);
  delete cut;

  return 1;
}
