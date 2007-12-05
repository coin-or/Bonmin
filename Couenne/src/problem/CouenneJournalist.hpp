// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptAlg.hpp 988 2007-06-01 21:57:27Z andreasw $
//
// Author:  Andreas Waechter           IBM    2007-12-04

// This file is a wrapper for the Journalist from the Ipopt project.
// The only thing it adds over the original Journalist class is that
// the names are easier to reach, and that the categories are given
// real names.

#ifndef CouenneJournalist_hpp
#define CouenneJournalist_hpp

#include "IpJournalist.hpp"

typedef Ipopt::SmartPtr<Ipopt::Journalist> JnlstPtr;

const Ipopt::EJournalCategory J_BRANCHING(Ipopt::J_USER1);
const Ipopt::EJournalCategory J_BOUNDTIGHENING(Ipopt::J_USER2);
const Ipopt::EJournalCategory J_CONVEXIFYING(Ipopt::J_USER3);
const Ipopt::EJournalCategory J_PROBLEM(Ipopt::J_USER4);

#endif
