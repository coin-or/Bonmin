// (C) Copyright International Business Machines Corporation and Carnegie Mellon University 2006 
// All Rights Reserved.
// This code is published under the Common Public License.
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
    BM_init user_init;
    return bcp_main(argc, argv, &user_init);
}

//#############################################################################

template <>
void BCP_parameter_set<BM_par>::create_keyword_list() {
    // Create the list of keywords for parameter file reading
    keys.push_back(make_pair(BCP_string("BranchOnSos"),
			     BCP_parameter(BCP_CharPar, BranchOnSos)));
    keys.push_back(make_pair(BCP_string("PureBranchAndBound"),
			     BCP_parameter(BCP_CharPar, PureBranchAndBound)));
    keys.push_back(make_pair(BCP_string("PrintBranchingInfo"),
			     BCP_parameter(BCP_CharPar, PrintBranchingInfo)));
    keys.push_back(make_pair(BCP_string("CombinedDistanceAndPriority"),
			     BCP_parameter(BCP_CharPar,
					   CombinedDistanceAndPriority)));
    keys.push_back(make_pair(BCP_string("SosWithLowPriorityMoreImportant"),
			     BCP_parameter(BCP_CharPar,
					   SosWithLowPriorityMoreImportant)));
    keys.push_back(make_pair(BCP_string("VarWithLowPriorityMoreImportant"),
			     BCP_parameter(BCP_CharPar,
					   VarWithLowPriorityMoreImportant)));
    keys.push_back(make_pair(BCP_string("NumNlpFailureMax"),
			     BCP_parameter(BCP_IntPar, NumNlpFailureMax)));
    keys.push_back(make_pair(BCP_string("WarmStartStrategy"),
			     BCP_parameter(BCP_IntPar, WarmStartStrategy)));
    keys.push_back(make_pair(BCP_string("NL_filename"),
			     BCP_parameter(BCP_StringPar, NL_filename)));
    keys.push_back(make_pair(BCP_string("IpoptParamfile"),
			     BCP_parameter(BCP_StringPar, IpoptParamfile)));
}

/****************************************************************************/

template <>
void BCP_parameter_set<BM_par>::set_default_entries() {
    set_entry(BranchOnSos, true);
    set_entry(PureBranchAndBound, false);
    set_entry(PrintBranchingInfo, true);
    set_entry(CombinedDistanceAndPriority, true);
    set_entry(SosWithLowPriorityMoreImportant, true);
    set_entry(VarWithLowPriorityMoreImportant, true);
    set_entry(NumNlpFailureMax, 5);
    set_entry(WarmStartStrategy, WarmStartFromRoot);
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
