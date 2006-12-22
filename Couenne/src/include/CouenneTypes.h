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
enum linearity_type {CONSTANT = 0, LINEAR, QUADRATIC, NONLINEAR};
enum pos            {PRE=0, POST, INSIDE, NONE};
enum con_sign       {COUENNE_EQ, COUENNE_LE, COUENNE_GE, COUENNE_RNG};
enum opt_sense      {MAXIMIZE, MINIMIZE};
enum conv_type      {CURRENT_ONLY, UNIFORM_GRID, AROUND_CURPOINT};

typedef double CouNumber;
typedef CouNumber       (*unary_function)   (CouNumber);

#endif
