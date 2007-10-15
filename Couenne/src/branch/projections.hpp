/*
 * Name:    projections.hpp
 * Authors: Pietro Belotti, Carnegie Mellon University
 * Purpose: tools for projecting points on lines/planes
 *
 * (C) Carnegie-Mellon University, 2006. 
 * This file is licensed under the Common Public License (CPL)
 */

#include "CouennePrecisions.hpp"


/** Compute projection of point (x0, y0) on the segment defined by
 *  line ax + by + c <>= 0 (sign provided by parameter sign) and
 *  bounds [lb, ub] on x. Return distance from segment, 0 if satisfied
 */

CouNumber project (CouNumber a,   CouNumber b, CouNumber c, 
		   CouNumber x0,  CouNumber y0, 
		   CouNumber lb,  CouNumber ub, 
		   int sign,
		   CouNumber *xp = NULL, CouNumber *yp = NULL);

/** Compute projection of point (x0, y0) on the segment defined by two
 *  points (x1,y1), (x2, y2) -- sign provided by parameter
 *  sign. Return distance from segment, 0 if on it.
 */

CouNumber projectSeg (CouNumber x0,  CouNumber y0, 
		      CouNumber x1,  CouNumber y1, 
		      CouNumber x2,  CouNumber y2, 
		      int sign,
		      CouNumber *xp = NULL, CouNumber *yp = NULL);
