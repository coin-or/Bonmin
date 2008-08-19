#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>

#include <stdlib.h>

#include "CoinTime.hpp"
#include "BonminConfig.h"
#include "BonCouenneInterface.hpp"
#include "BonIpoptSolver.hpp"

#include "CoinHelperFunctions.hpp"
#include "BonCouenneSetup.hpp"

#include "BonCbc.hpp"

#include "CbcCutGenerator.hpp"
#include "CouenneProblem.hpp"
#include "CouenneCutGenerator.hpp"

#include "asl.h"
#include "nlp.h"
#include "getstub.h"
#include "r_opn.hd" // for N_OPS 
#include "opcode.hd"

int pint;

static keyword keywds[] = { /* must be sorted */
  { const_cast<char*>("barrier"), 
    I_val,  
    &pint,
    (const char *) "non faccio nulla"}
};

extern Option_Info Oinfo;

//extern 
//Option_Info* TheOInfo;
//extern ASL* asl;

using namespace Bonmin;

int main (int argc, char *argv[]) {

  using namespace Ipopt;

  //  char * pbName = NULL;
  //  double time_start = CoinCpuTime();

  Bonmin::Bab bb;
  bb.setUsingCouenne (true);

  CouenneSetup bonmin;
  bonmin.InitializeCouenne (argv);

  //

  SmartAsl *aslfg = new SmartAsl;
  aslfg -> asl = readASLfg (argv);

  int NumberOfVariables = bonmin.couennePtr () -> Problem () -> nOrigVars ();

  typedef struct {char *msg; int code, wantsol;} Sol_info;

  static Sol_info solinfo [] = {
    {"Optimal", 0, 1}
  };

  static SufDecl suftab [] = {
    {(const char *) "newlb", 0, ASL_Sufkind_var | ASL_Sufkind_real, 0},
    {(const char *) "newub", 0, ASL_Sufkind_var | ASL_Sufkind_real, 0}};

  char oinf = 0;

  //Option_Info *TheOInfo = (Option_Info*) &oinf;

  Sol_info *Si;

  suf_declare_ASL (aslfg -> asl, suftab, 2); //sizeof (suftab) / sizeof (Sufdecl));

  // add an AMPL suffix
  SufDesc* vnewLb = suf_get_ASL(aslfg -> asl, "newlb", ASL_Sufkind_var);
  SufDesc* vnewUb = suf_get_ASL(aslfg -> asl, "newub", ASL_Sufkind_var);

  double *x = new double [NumberOfVariables];

  CoinFillN (x, NumberOfVariables, 0.);

  write_sol_ASL (aslfg -> asl, (const char *) "solved!", 0, 0, &Oinfo);

  delete [] x;

  vnewLb -> u.r = (real*)M1zapalloc_ASL(&aslfg -> asl->i, NumberOfVariables * sizeof(real));
  vnewUb -> u.r = (real*)M1zapalloc_ASL(&aslfg -> asl->i, NumberOfVariables * sizeof(real));

  const double
    *newL = bonmin.couennePtr () -> Problem () -> Lb (),
    *newU = bonmin.couennePtr () -> Problem () -> Ub ();

  // return variable integrality
  for(int i = 0; i < NumberOfVariables; i++) {
    vnewLb->u.r[i] = newL [i];
    vnewUb->u.r[i] = newU [i];
  }

  return 0;
}
