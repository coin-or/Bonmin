/*
 * Name:    SdpCutGen.hpp
 * Author:  Pietro Belotti
 * Purpose: Cut separator for PSD matrices of the form X\succeq x x^T
 *
 *
 * (C) Pietro Belotti. This file is published under the Common Public License
 */

#ifndef SDPCUTGEN_HPP
#define SDPCUTGEN_HPP

#include <CglCutGenerator.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))


/// Class to separate SDP cuts

class SdpCutGen: public CglCutGenerator {

private:

  int  n_;        /// order of the matrix
  int *indices_;  /// indices to be tracked
  double  *b_;    /// value of original variables
  double **Q_;    /// value of RLT variables Xij = xi * xj

public:

  /// constructor
  SdpCutGen  (int, double *, double **);

  /// copy constructor
  SdpCutGen  (const SdpCutGen &);

  /// destructor
  ~SdpCutGen ();

  /// clone 
  SdpCutGen *clone () const
  {return new SdpCutGen (*this);}

  /// the main cut generator
  void generateCuts (const OsiSolverInterface &,  OsiCuts &, 
		     const CglTreeInfo = CglTreeInfo ()) const;

  /// separation procedures
  void separateEV  (const OsiSolverInterface &, OsiCuts &) const;
  void separateLU  (const OsiSolverInterface &, OsiCuts &) const;
  void separateBK  (const OsiSolverInterface &, OsiCuts &) const;
  void separateTRM (const OsiSolverInterface &, OsiCuts &) const;

  /// play with eigenvalues/vectors 
  void eigenPlay (OsiCuts &, int n, const double *, 
		  int m, double *vector, double *value) const;

  /// insert a SDP cut a' X b >= 0 given vectors a and b and size of the matrix
  void genSDPcut (OsiCuts &, double *, double *) const;
};


/// wrapper for Lapack's Fortran routine to compute eigenvalues/vectors
void dsyevx_wrapper (int, double *, int &, double * &, double * &);

/// partially decompose matrix (leave all 0 pivots at lower right
/// submatrix)
extern "C" {
  double **partialLU (double **, int);
}

#endif
