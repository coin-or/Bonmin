/*
 * Name:    CouennePrecisions.hpp
 * Author:  Pietro Belotti
 * Purpose: constants for evaluation procedures
 *
 * (C) Carnegie-Mellon University, 2006-08. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_PRECISIONS_HPP
#define COUENNE_PRECISIONS_HPP

#include <math.h>

// must be >= 1e-7
#define COUENNE_EPS           1.e-7

// for integrality check
#define COUENNE_EPS_INT       1.e-9

// for simplification
#define COUENNE_EPS_SIMPL     1.e-20

// for bounds
#define COUENNE_INFINITY      1.e+50

// for cuts, ensures stability and scaling
#define COU_MAX_COEFF     1.e+9

// for cuts, ditto
#define COU_MIN_COEFF     1.e-9

// rounds to nearest integer
#define COUENNE_round(x) ((int) (floor ((x) + 0.5)))

#define MAX_BOUND 1.e45

#endif
