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

/** type of a node in an expression tree */
enum nodeType       {CONST=0, VAR, UNARY, N_ARY, COPY, AUX, EMPTY};

/** linearity of an expression, as returned by the method Linearity() */
enum linearity_type {ZERO=0, CONSTANT, LINEAR, QUADRATIC, NONLINEAR};

/** position where the operator should be printed when printing the expression
 *
 *  For instance, it is INSIDE for exprSum, exprMul, exprDiv, while it
 *  is PRE for exprLog, exprSin, exprExp...
 **/
enum pos            {PRE=0, POST, INSIDE, NONE};

/** sign of constraint */
enum con_sign       {COUENNE_EQ, COUENNE_LE, COUENNE_GE, COUENNE_RNG};

/** optimization sense of an objective function */
enum opt_sense      {MAXIMIZE, MINIMIZE};

/** position and number of convexification cuts added for a lower
    convex (upper concave) envelope */
enum conv_type      {CURRENT_ONLY, UNIFORM_GRID, AROUND_CURPOINT};

/** code returned by the method expression::code() */
enum expr_type {/*COU_EXPRAUX,  COU_EXPRCLONE, COU_EXPRCOPY, */
                COU_EXPRESSION, /***** variables, constants **************/
		COU_EXPRCONST, COU_EXPRVAR, COU_EXPRLBOUND, COU_EXPRUBOUND, 
		/*COU_EXPRIVAR, */ 
		COU_EXPROP,     /***** n-ary operators *******************/
		COU_EXPRSUB,  COU_EXPRSUM, COU_EXPRGROUP, COU_EXPRQUAD,
		COU_EXPRMIN,  COU_EXPRMUL, COU_EXPRPOW, COU_EXPRMAX, COU_EXPRDIV, 
		/*COU_EXPRBDIV,  COU_EXPRBMUL,*/ 
		COU_EXPRUNARY,  /***** unary operators *******************/
		COU_EXPRCOS,  COU_EXPRABS,
		COU_EXPREXP,  COU_EXPRINV,   COU_EXPRLOG,    
		COU_EXPROPP,   COU_EXPRSIN
		};

/** convexity type of an expression */
enum convexity {UNSET, NONCONVEX, CONVEX, CONCAVE, AFFINE};

/** (un)changed bounds */

#define UNCHANGED 0
#define CHANGED   1
#define EXACT     2

/** status of lower/upper bound of a variable, to be checked/modified
    in bound tightening */
typedef struct chg_bounds {
  char lower;
  char upper;
} t_chg_bounds;

/** main number type in Couenne */
typedef double CouNumber;

/** unary function, used in all exprUnary */
typedef CouNumber (*unary_function) (CouNumber);

#endif
