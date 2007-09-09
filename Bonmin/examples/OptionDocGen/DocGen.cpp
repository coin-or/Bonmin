// (C) Copyright Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  03/17/2006


#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#endif
#include <iomanip>
#include <fstream>

#include "CoinTime.hpp"

#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "MyTMINLP.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"


int main (int argc, char *argv[])
{
  using namespace Ipopt;
  using namespace Bonmin;
  SmartPtr<MyTMINLP> tminlp = new MyTMINLP;
  

  BonminSetup bonmin;
  bonmin.initializeOptionsAndJournalist();
  //Now initialize from tminlp
  bonmin.initialize(GetRawPtr(tminlp));

  std::ofstream of("Table.tex");
  bonmin.roptions()->writeLatexOptionsTable(of, Bonmin::RegisteredOptions::BonminCategory);
  return 0;
}

