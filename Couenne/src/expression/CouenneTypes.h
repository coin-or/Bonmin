/*
 * Name:    CouenneTypes.h
 * Author:  Pietro Belotti
 * Purpose: define number types used throughout the code
 *
 * (C) 2006 Pietro Belotti, Carnegie Mellon University.
 * This file is licensed under the Common Public License (CPL).
 */

#ifndef COUENNE_TYPES_H
#define COUENNE_TYPES_H

enum nodeType       {CONST=0, VAR, UNARY, N_ARY, COPY, AUX, EMPTY};
enum linearity_type {ZERO=0, CONSTANT, LINEAR, QUADRATIC, NONLINEAR};
enum pos            {PRE=0, POST, INSIDE, NONE};
enum con_sign       {COUENNE_EQ, COUENNE_LE, COUENNE_GE, COUENNE_RNG};
enum opt_sense      {MAXIMIZE, MINIMIZE};
enum conv_type      {CURRENT_ONLY, UNIFORM_GRID, AROUND_CURPOINT};

enum expr_type {/*COU_EXPRAUX,  COU_EXPRCLONE, COU_EXPRCOPY, */
                COU_EXPRESSION, /***** variables, constants **************/
		COU_EXPRCONST, COU_EXPRVAR, COU_EXPRLBOUND, COU_EXPRUBOUND, 
		/*COU_EXPRIVAR, */ 
		COU_EXPROP,     /***** n-ary operators *******************/
		COU_EXPRSUB,  COU_EXPRSUM, COU_EXPRGROUP, 
		COU_EXPRMIN,  COU_EXPRMUL, COU_EXPRPOW, COU_EXPRMAX, COU_EXPRDIV, 
		/*COU_EXPRBDIV,  COU_EXPRBMUL,*/ 
		COU_EXPRUNARY,  /***** unary operators *******************/
		COU_EXPRCOS,  COU_EXPRABS,
		COU_EXPREXP,  COU_EXPRINV,   COU_EXPRLOG,    
		COU_EXPROPP,   COU_EXPRSIN
		};

enum convexity {NONCONVEX, CONVEX, CONCAVE, AFFINE};

typedef double CouNumber;
typedef CouNumber (*unary_function) (CouNumber);

#endif
