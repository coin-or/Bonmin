#ifndef SDPCUTGEN_HPP
#define SDPCUTGEN_HPP

#include <CglCutGenerator.hpp>

#define indexQ(i,j,n) ((n) + (i) * (2*(n)-1-(i)) / 2 + (j))

class SdpCutGen: public CglCutGenerator {

private:

  /// order of the matrix
  int n_; 

  /// indices to be tracked
  int     *indices_;

  /// value of original variables
  double  *b_;

  /// value of RLT variables Xij = xi * xj
  double **Q_;

public:

  /// constructor
  SdpCutGen  (int n, double *b, double **Q):
    n_ (n) {

    b_ = new double   [n];
    Q_ = new double * [n];
    for (int i=n; i--;)
      Q_ [i] = new double [n];

    for (int i=n; i--;) {
      b_ [i] = b [i];
      for (int j=n; j--;) 
	Q_ [i] [j] = Q [i] [j];
    }
  }

  /// copy constructor
  SdpCutGen  (const SdpCutGen &rhs):
    n_ (rhs.n_) {

    b_ = new double   [n_];
    Q_ = new double * [n_];
    for (int i=n_; i--;)
      Q_ [i] = new double [n_];

    for (int i=n_; i--;) {
      b_ [i] = rhs.b_ [i];
      for (int j=n_; j--;) 
	Q_ [i] [j] = rhs.Q_ [i] [j];
    }
  }

  /// destructor
  ~SdpCutGen () {
    free (b_); 
    while (n_--) 
      free (Q_ [n_]);
    free (Q_);
  }

  /// clone 
  SdpCutGen *clone () const
  {return new SdpCutGen (*this);}

  /// the main cut generator
  void generateCuts (const OsiSolverInterface &, 
		     OsiCuts &, 
		     const CglTreeInfo = CglTreeInfo ()) const;

  void separateEV  (const OsiSolverInterface &, OsiCuts &) const;
  void separateLU  (const OsiSolverInterface &, OsiCuts &) const;
  void separateBK  (const OsiSolverInterface &, OsiCuts &) const;
  void separateTRM (const OsiSolverInterface &, OsiCuts &) const;

  /// play with eigenvalues/vectors 
  void eigenPlay (OsiCuts &, int n, int m, double *vector, double *value) const;

  /// insert a SDP cut a' X a >= 0 given vector a and size of the matrix
  void genSDPcut (OsiCuts &, double *, double *) const;

  ///
  void dsyevx_wrapper (int, double *, int &, double *, double *) const;
};


/// partially decompose matrix (leave all 0 pivots at lower right
/// submatrix)
extern "C" {
  double **partialLU (double **, int);
}
#endif
