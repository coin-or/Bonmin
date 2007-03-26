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

#define COUENNE_EPS       1e-7 /* keep it at least 1e-7 or strange things happen */
#define COUENNE_EPS_SIMPL 1e-20
#define COUENNE_INFINITY  1e+50

#define FELINE_round(x) ((int) (floor ((x)+0.5)))

#endif
