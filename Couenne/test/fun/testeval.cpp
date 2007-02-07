//
// File: testeval.cpp
// Author: Pietro Belotti
// Purpose: test function evaluation library
//
// (C) Pietro Belotti, 2006, Carnegie Mellon University
// This file is distributed under the Common Public License (CPL)
//

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

#include <stdio.h>
#include <iostream>

#include <math.h>

//#include <expressions.h>
#include <expression.h>
#include <exprMul.h>
#include <exprSin.h>
#include <exprSum.h>
#include <exprSub.h>
#include <exprExp.h>
#include <exprPow.h>
#include <exprLog.h>
#include <CouenneTypes.h>
#include <CouenneProblem.h>
#include <rootQ.h>

/*
void test (expression * &root_) {
  expression *rs = root_ -> simplify ();
  if (!rs) rs =  root_;
  else      delete root_;
  //  rs -> print (std::cout);
  //  printf ("\n");
  root_ = rs;
}
*/


int main (int argc, char **argv) {

  CouNumber variables [1000] = { 1.5,    3.3, 45.6};
  CouNumber l_bounds  [1000] = {-1.0,    1.3,  5.6};
  CouNumber u_bounds  [1000] = { 3.0,    3.4, 55.6};

  /*
  Qroot qq;
  for (int i=100000; i--;)
    for (int k = 12; k-- > 1;) {
      if (i==1) 
        printf ("root %3d: %.6f\n", k, qq (k));
      qq (k);
    }
  */

  /*
  exprVar var0 (0);
  exprVar var1 (1);
  exprVar var2 (2);

  exprConst k1 (1.4);
  exprConst k2 (1.0);
  exprConst k3 (9);

  expression *al1 [] = {&k1, &var0};  
  expression *al2 [] = {&k2, &var1};
  expression *al3 [] = {&k3, &var2};

  exprOp divi (al1, 2, &oDiv);
  exprOp mult (al2, 2, &oMul);
  exprOp sum2 (al3, 2, &oAdd);

  expression *al0 [] = {&divi, &mult, &sum2};

  exprOp root (al0, 3, &oAdd);
  */

  exprVar *vx = new exprVar (0), 
          *vy = new exprVar (1);

  expression **al0 = new expression * [2];
  al0 [0] = new exprCopy (vx);
  al0 [1] = new exprCopy (vy);

  exprConst *c13 = new exprConst (13);

  expression **al1 = new expression * [2];
  al1 [0] = c13;
  al1 [1] = new exprSub (al0, 2);

  exprMul *mult1 = new exprMul (al1, 2);

  exprSin *sin1 = new exprSin (mult1);

  expression **al2 = new expression * [2];
  al2 [0] = new exprCopy (vx);
  al2 [1] = new exprCopy (vx);

  exprMul *mult2 = new exprMul (al2, 2);

  exprExp *exp1 = new exprExp (new exprCopy (vx));

  expression **al3 = new expression * [2];
  al3 [0] = new exprCopy (vy);
  al3 [1] = exp1;

  exprMul *mult3 = new exprMul (al3, 2);

  expression **al4 = new expression * [2];
  al4 [0] = mult2;
  al4 [1] = mult3;

  //  exprSum *root = new exprSum (al4, 2);

  exprSum *sum4 = new exprSum (al4, 2);

  expression **al5 = new expression * [2];
  al5 [0] = sin1;
  al5 [1] = sum4;

  exprMul *mult4 = new exprMul (al5, 2);

  exprConst *c4 = new exprConst (4);
  expression **al6 = new expression * [2];
  al6 [0] = new exprCopy (vx);
  al6 [1] = c4;

  exprPow *pow1 = new exprPow (al6, 2);

  exprConst *c1 = new exprConst (1);

  expression **al7 = new expression * [2];
  al7 [0] = c1; 
  al7 [1] = pow1;

  exprSum *sum5 = new exprSum (al7, 2);
  exprLog *log1 = new exprLog (sum5);

  exprConst *c3 = new exprConst (3);
  expression **al8 = new expression * [2];
  al8 [0] = c3;
  al8 [1] = log1;

  exprMul *mult5 = new exprMul (al8, 2);

  expression **al9 = new expression * [2];
  al9 [0] = mult4;
  al9 [1] = mult5;

  exprSum *root = new exprSum (al9, 2);

  expression::update (variables, l_bounds, u_bounds);

  //  CouenneProblem *p = new CouenneProblem ();

  /*
  expression *root2 = 
    new exprSum 
     (new exprMul (new exprClone (vx), new exprClone (vy)),
      new exprPow (new exprConst (4), new exprClone (vx)));

  expression *con0 = 
    new exprSum 
     (new exprSin (new exprSum (new exprClone (vy), new exprExp (new exprClone (vy)))),
      new exprInv (new exprClone (vx)));

  expression *con1 = 
    new exprSum 
     (new exprMul (new exprClone (vx), new exprClone (vy)),
      new exprLog (new exprClone (vx)));

  expression *con2 = 
    new exprSub
     (new exprExp (new exprClone (vy)),
      new exprPow (new exprClone (vx), new exprConst (4)));

  p -> addObjective (root2, "max");

  for (int i=0; i<100; i++) {
    p -> addLEConstraint (new exprClone (con0), new exprConst (i*3));
    p -> addEQConstraint (new exprClone (con1), new exprConst (i*4 % 3));
    p -> addGEConstraint (new exprClone (con2), new exprConst (2 + 0.01*i));
  }
  */

  ///////////////////////////////////////////////////////////

  //  expression *con1 = new exprPow (new exprConst (3), new exprClone (vy));
  /*
  expression *con1 = new exprSum (new exprMul (new exprClone (vy), new exprClone (vx)),
				  new exprLog (new exprClone (vx)));

  expression *cons = con1 -> simplify ();

  if (cons) 
    con1 = cons;

  p -> addLEConstraint (con1, new exprConst (4));
  //  p -> addGEConstraint (con1, new exprConst (1));

  printf ("-------------------Original:\n");
  p -> print (std::cout);

  p -> standardize ();

  printf ("-------------------Standardized:\n");
  p -> print (std::cout);

  p -> convexify ();

  printf ("-------------------Convexified:\n");
  p -> print (std::cout);
  */

  /*
  struct rusage u;

  getrusage (RUSAGE_SELF, &u);
  double now = u.ru_utime.tv_sec + 1e-6 * u.ru_utime.tv_sec;

  int neval = 0;

  for (int i=0; i<10000; i++) {

    for (register int j = p -> nCons (); j--;) {
      expression ** coeffs = p -> Con (j) -> Coeff ();
      for (register int k = p -> Con (j) -> nTerms (); k--;neval++)
	(*(coeffs [k])) ();
    }
  }

  getrusage (RUSAGE_SELF, &u);
  fprintf (stderr, "time = %.3f for %d evaluations in full matrix update\n", 
	   now = (double) u.ru_utime.tv_sec + 1.0e-6 * u.ru_utime.tv_usec - now, neval);
  */
  /*
  for (int i=0; i<4; i++) {

    printf ("(%.6f, %.6f) (%.6f, %.6f)\n", 
	    l_bounds [0], u_bounds [0], l_bounds [1], u_bounds [1]);

    for (register int j = 0; j < p -> nCons (); j++) {
      expression ** coeffs = p -> Con (j) -> Coeff ();
      for (register int k = 0; k < p -> Con (j) -> nTerms (); k++) {
	printf (" %+.6f x_%d  ", 
		(*(coeffs [k])) (), 
		p -> Con (j) -> Indices () [k]);
	(*(coeffs [k])) ();
      }

      CouNumber rhs = 0;

      if (p -> Con (j) -> Rhs ())
	rhs = (*(p -> Con (j) -> Rhs ())) ();

      printf (" >=< %.6f\n", rhs);
    }

    l_bounds [0] += (u_bounds [0] - l_bounds [0]) / 10;
    u_bounds [0] -= (u_bounds [0] - l_bounds [0]) / 10;

    l_bounds [1] += (u_bounds [1] - l_bounds [1]) / 5;
    u_bounds [1] -= (u_bounds [1] - l_bounds [1]) / 100;
  }

  return 0;

  delete p;
  */

  (*vx) ();
  (*vy) ();

  //  printf ("Root\n");

  //  root -> print (std::cout); printf ("\n");

  //  printf ("Derivatives\n1: ");
  /*

  expression *root_x, *root_y, *rxs, *rys;

  root_x = root -> differentiate (0);
  root_y = root -> differentiate (1);

  //  root_x -> print (std::cout); printf ("\n2: ");
  //  root_y -> print (std::cout); printf ("\n");

  //  printf ("Simplified\n");

  //  printf ("First order:\n");

  test (root_x);
  test (root_y);

  //  printf ("\nSecond order:\n");

  expression *root_xx = root_x -> differentiate (0);
  expression *root_xy = root_x -> differentiate (1);
  expression *root_yx = root_y -> differentiate (0);
  expression *root_yy = root_y -> differentiate (1);

  //  printf ("\n xx: ");
  //  root_xx -> print (std::cout); printf ("\n xy: ");
  //  root_xy -> print (std::cout); printf ("\n yx: ");
  //  root_yx -> print (std::cout); printf ("\n yy: ");
  //  root_yy -> print (std::cout); printf ("\n\n");

  test (root_xx);
  test (root_xy);
  test (root_yx);
  test (root_yy);

  //  root_x -> print (std::cout); printf ("\n2: ");
  //  root_y -> print (std::cout); printf ("\n");
  */

  //  /*return (sin (c1*(x-y)) * (x*x + y*exp (x)) + c2 * log (c3 + pow (x,c4))); */

  //  CouNumber *x = variables;
  //  CouNumber *y = variables + 1;

  root -> print (std::cout); printf ("\n");

  if (1)
  for (register int i=1000000; i--;) {

    /*
    CouNumber c1 = 13, c2 = 3, c3 = 1, c4 = 4;

    *variables += 0.3;
    variables [1] += 0.2;
    variables [2] += 0.1;

    (*vx) (variables);
    (*vy) (variables);
    */

    //    printf ("\r%14.5f", root (variables)); fflush (stdout);

    // std::cout << "calc expr: " << root (variables) << std::endl;

    //    printf ("%20.4f %20.4f\n",
    //    (*root) (variables)
      //	    (sin (c1*(*x-*y)) * (*x**x + *y*exp (*x)) + c2 * log (c3 + pow (*x,c4))));

    (*root) ();

	/*
    (*root_x) ();
    (*root_y) ();

    (*root_xx) ();
    (*root_xy) ();
    (*root_yx) ();
    (*root_yy) ();
	*/

    //    if ((*root_x) () < 2.55);
    //    if ((*root_y) () < 2.55);
  }

  delete root;
  delete vx;
  delete vy;
  /*
  delete root_x;
  delete root_y;

  delete root_xx;
  delete root_yx;
  delete root_xy;
  delete root_yy;
  */
}
