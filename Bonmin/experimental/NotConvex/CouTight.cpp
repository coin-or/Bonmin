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

static int empty_int = 0;

static keyword keywds[] = { /* must be sorted */
  KW(const_cast<char*>("empty"), 
     I_val, 
     empty_int, 
     const_cast<char*>("nothing")), 
};

extern Option_Info Oinfo;

using namespace Bonmin;

int main (int argc, char *argv[]) {

  using namespace Ipopt;

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
    {
      const_cast<char*>("Optimal"), 
      0, 
      1
    }
  };

  static SufDecl suftab [] = {
    {const_cast<char*>("newlb"), 0, 
     ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_output, 0},
    {const_cast<char*>("newub"), 0, 
     ASL_Sufkind_var | ASL_Sufkind_real | ASL_Sufkind_output, 0}};

  char oinf = 0;

  Sol_info *Si;

  suf_declare_ASL (aslfg -> asl, suftab, sizeof (suftab) / sizeof (SufDecl));

  // add an AMPL suffix
  SufDesc* vnewLb = suf_get_ASL(aslfg -> asl, "newlb", ASL_Sufkind_var);
  SufDesc* vnewUb = suf_get_ASL(aslfg -> asl, "newub", ASL_Sufkind_var);

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

  write_sol_ASL (aslfg -> asl, const_cast<char*>("solved!"), 0, 0, &Oinfo);
  return 0;
}
