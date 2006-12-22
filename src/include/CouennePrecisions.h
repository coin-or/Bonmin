/*
 * Name:    CouennePrecisions.h
 * Author:  Pietro Belotti
 * Purpose: constants for evaluation procedures
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_PRECISION_H
#define COUENNE_PRECISION_H

#include <CouenneTypes.h>
#include <math.h>

const CouNumber COUENNE_EPS       = 1e-100;
const CouNumber COUENNE_EPS_SIMPL = 1e-200;
const CouNumber COUENNE_INFINITY  = 1e+300;

inline int FELINE_round (double x) {
  return (int) (ceil (x - 0.5));
}

#endif
