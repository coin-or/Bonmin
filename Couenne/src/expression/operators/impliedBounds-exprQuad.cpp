/*
 * Name:    impliedBounds-exprQuad.cpp
 * Author:  Pietro Belotti
 * Purpose: inferring bounds on independent variables of an exprQuad
 *          given bounds on the auxiliary varible
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <exprQuad.hpp>

/// implied bound processing for quadratic form upon change in lower-
/// and/or upper bound of w, whose index is wind

bool exprQuad::impliedBound (int wind, CouNumber *l, CouNumber *u, t_chg_bounds *chg) {

  return false;
}
