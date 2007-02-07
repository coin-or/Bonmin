/*
 * Name:    exprOpp.C
 * Author:  Pietro Belotti
 * Purpose: definition of the opposite -f(x) of a function
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprOpp.h>


// find bounds of -x given bounds on x

void exprOpp::getBounds (expression *&lb, expression *&ub) {

    expression *lba, *uba;
    argument_ -> getBounds (lba, uba);

    lb = new exprOpp (uba);
    ub = new exprOpp (lba);
  }


// differentiation

inline expression *exprOpp::differentiate (int index) 
{return new exprOpp (argument_ -> differentiate (index));}


// printing

void exprOpp::print (std::ostream& out) const
{exprUnary::print (out, "-", PRE);}
