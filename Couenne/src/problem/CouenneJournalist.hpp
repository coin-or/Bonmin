// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: CouenneJournalist.hpp 988 2007-06-01 21:57:27Z andreasw $
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
typedef Ipopt::SmartPtr<const Ipopt::Journalist> ConstJnlstPtr;

const Ipopt::EJournalCategory J_BRANCHING       (Ipopt::J_USER1);
const Ipopt::EJournalCategory J_BOUNDTIGHTENING (Ipopt::J_USER2);
const Ipopt::EJournalCategory J_CONVEXIFYING    (Ipopt::J_USER3);
const Ipopt::EJournalCategory J_PROBLEM         (Ipopt::J_USER4);
const Ipopt::EJournalCategory J_NLPHEURISTIC    (Ipopt::J_USER5);
const Ipopt::EJournalCategory J_DISJCUTS        (Ipopt::J_USER6);
const Ipopt::EJournalCategory J_REFORMULATE     (Ipopt::J_USER7);

#endif
