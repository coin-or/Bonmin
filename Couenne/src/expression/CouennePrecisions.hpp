/*
 * Name:    CouennePrecisions.hpp
 * Author:  Pietro Belotti
 * Purpose: constants for evaluation procedures
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_PRECISIONS_HPP
#define COUENNE_PRECISIONS_HPP

#include <math.h>

/* keep it at least 1e-7, or strange things happen */
#define COUENNE_EPS           1e-7

#define COUENNE_EPS_SIMPL     1e-20

#define COUENNE_INFINITY      1.0e+50

#define COU_MAX_COEFF     1.0e+12

#define COUENNE_round(x) ((int) (floor ((x) + 0.5)))

#endif
