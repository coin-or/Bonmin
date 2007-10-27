/*
 * Name:    CouObjStats.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: Statistics on CouenneObjects (branches created, where, when, on which variable)
 *
 * (C) Carnegie-Mellon University, 2006-07.
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUOBJSTATS_HPP
#define COUOBJSTATS_HPP

#include <vector>

class CouObjStats {

public:

  CouObjStats () {}

  CouObjStats (const CouObjStats &) {}

  CouObjStats *clone () 
  {return new CouObjStats (*this);}

  ~CouObjStats () {}

protected:

  static int Index_;

  // data for all calls
  std::vector <double> infeasibility_;

  std::vector <double> lb_;
  std::vector <double> ub_;
  std::vector <double> brpoint_;
  std::vector <int>    index_;
};

#endif
