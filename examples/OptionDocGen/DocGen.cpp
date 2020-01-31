// (C) Copyright Carnegie Mellon University 2006, 2007
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// P. Bonami, Carnegie Mellon University
//
// Date :  03/17/2006


#include <iomanip>
#include <fstream>

#include "CoinPragma.hpp"
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

  std::ofstream of("options_list_bonmin_content.tex");
  bonmin.roptions()->writeLatexHtmlDoc(of, Bonmin::RegisteredOptions::BonminCategory);
  of.close();
  of.open("options_list_ipopt_content.tex");
  bonmin.roptions()->writeLatexHtmlDoc(of, Bonmin::RegisteredOptions::IpoptCategory);
  of.close();
  of.open("options_list_filter_content.tex");
  bonmin.roptions()->writeLatexHtmlDoc(of, Bonmin::RegisteredOptions::FilterCategory);
  of.close();

  of.open("options_table.tex");
  bonmin.roptions()->writeLatexOptionsTable(of, Bonmin::RegisteredOptions::BonminCategory);
  of.close();

  of.open("bonmin.opt");
  bonmin.roptions()->writeBonminOpt(of, Bonmin::RegisteredOptions::BonminCategory);
  return 0;
}

