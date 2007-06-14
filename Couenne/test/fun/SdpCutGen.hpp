#ifndef SDPCUTGEN_HPP
#define SDPCUTGEN_HPP

#include <CglCutGenerator.hpp>

class SdpCutGen: public CglCutGenerator {

private:

  int n_;
  double *b_;
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
};

#endif
