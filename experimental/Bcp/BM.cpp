// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006, 2007 
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors :
// Laszlo Ladanyi, International Business Machines Corporation
// Pierre Bonami, Carnegie Mellon University,
//
// Date : 03/15/2006
#include <sstream>

/* To get the cumulative time spent on a processor just use a gawk command
   like this below. Look at the output first; probably the process id needs
   to be prepended to the regexp and the procid may also change the $7 to
   some other word.
   gawk -e 'BEGIN {t=0} /^BCP_lp: Time spent in this node:/ {t+=$7} END {print t}' outfile
*/

#include "BM.hpp"

using namespace std;

//#############################################################################

int main(int argc, char* argv[])
{
    CoinError::printErrors_ = true;
    BM_init user_init;
    int retcode = -1;
#if 1
    retcode = bcp_main(argc, argv, &user_init);
#else
    try {
      retcode = bcp_main(argc, argv, &user_init);
    }
    catch(Bonmin::TNLPSolver::UnsolvedError &E) {
      //      E.writeDiffFiles();
      E.printError(std::cerr);
   }
   catch(Bonmin::OsiTMINLPInterface::SimpleError &E) {
     std::cerr<<E.className()<<"::"<<E.methodName()
              <<std::endl
              <<E.message()<<std::endl;
   }
   catch(CoinError &E) {
     std::cerr<<E.className()<<"::"<<E.methodName()
              <<std::endl
              <<E.message()<<std::endl;
   }
   catch (Ipopt::OPTION_INVALID &E)
   {
    std::cerr<<"Ipopt exception : "<<E.Message()<<std::endl;
   }
   catch(...) {
     std::cerr<<" unrecognized exception"<<std::endl;
     throw;
   }
#endif

    return retcode;
}

//#############################################################################

template <>
void BCP_parameter_set<BM_par>::create_keyword_list() {
    // Create the list of keywords for parameter file reading
    keys.push_back(make_pair(BCP_string("BM_DisregardPriorities"),
			     BCP_parameter(BCP_CharPar, DisregardPriorities)));
    keys.push_back(make_pair(BCP_string("BM_PrintBranchingInfo"),
			     BCP_parameter(BCP_CharPar, PrintBranchingInfo)));
    keys.push_back(make_pair(BCP_string("BM_UsePseudoCosts"),
			     BCP_parameter(BCP_IntPar, UsePseudoCosts)));
    keys.push_back(make_pair(BCP_string("BM_DecreasingSortInSetupList"),
			     BCP_parameter(BCP_IntPar, DecreasingSortInSetupList)));
    keys.push_back(make_pair(BCP_string("BM_PreferHighCombinationInBranching"),
			     BCP_parameter(BCP_IntPar, PreferHighCombinationInBranching)));
    keys.push_back(make_pair(BCP_string("BM_NumNlpFailureMax"),
			     BCP_parameter(BCP_IntPar, NumNlpFailureMax)));
    keys.push_back(make_pair(BCP_string("BM_NL_filename"),
			     BCP_parameter(BCP_StringPar, NL_filename)));
    keys.push_back(make_pair(BCP_string("BM_IpoptParamfile"),
			     BCP_parameter(BCP_StringPar, IpoptParamfile)));
}

/****************************************************************************/

template <>
void BCP_parameter_set<BM_par>::set_default_entries() {
    set_entry(DisregardPriorities, false);
    set_entry(PrintBranchingInfo, true);
    set_entry(UsePseudoCosts, 1);
    set_entry(DecreasingSortInSetupList, 1);
    set_entry(PreferHighCombinationInBranching, 0);
    set_entry(NumNlpFailureMax, 5);
    set_entry(NL_filename, "");
    set_entry(IpoptParamfile, "");
}

//#############################################################################

BCP_lp_user *
BM_init::lp_init(BCP_lp_prob& p)
{
    return new BM_lp;
}

/****************************************************************************/

BCP_tm_user *
BM_init::tm_init(BCP_tm_prob& p,
                 const int argnum, const char * const * arglist)
{
    BM_tm* tm = new BM_tm;

    if (argnum == 2) {
	tm->par.read_from_file(arglist[1]);
    } else if (argnum == 1) {
	// work with defaults
    } else {
	tm->par.read_from_arglist(argnum, arglist);
    }

    tm->readIpopt();

    return tm;
}

/****************************************************************************/

BCP_user_pack *
BM_init::packer_init(BCP_user_class* p)
{
    return new BM_pack;
}

/****************************************************************************/

BM_stats::~BM_stats()
{
  // LACI: It would be nice to also print the process ID here, but I
  // wasn't sure how to get it...

  printf("Stats: #NodesSol = %d #SbSol = %d #Fixed = %d #SbDone = %d SumInd = %d SumPos = %e\n", numberNodeSolves_, numberSbSolves_, numberFixed_, numberStrongBranching_, sumStrongBranchingListIndices_, sumStrongBranchingListPositions_);
}
