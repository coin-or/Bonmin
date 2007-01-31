/*
 * Name:    exprAux.h
 * Author:  Pietro Belotti
 * Purpose: definition of the auxiliary variable class (used in
 *          standardization and convexification)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#ifndef COUENNE_EXPRAUX_H
#define COUENNE_EXPRAUX_H

#include <iostream>
#include <exprVar.h>
#include <CouenneTypes.h>
#include <CglCutGenerator.hpp>


// expression base class

class exprAux: public exprVar {

 protected:

  // The expression associated with this auxiliary variable
  expression *image_;

  // bounds determined by the associated expression's bounds and the
  // relative constraint's bound
  expression *lb_;
  expression *ub_;

 public:

  // Node type
  inline enum nodeType Type () 
    {return AUX;}

  // Constructor
  exprAux (expression *image, int index): 
    exprVar (index),
    image_  (image)

    {image_ -> getBounds (lb_, ub_);}

  // Destructor
  ~exprAux () {
    delete image_; 
    delete lb_; 
    delete ub_;
  }

  // copy constructor
  exprAux (const exprAux &e):
    exprVar (e.Index ()),
    image_  (e.Image () -> clone ())

    {image_ -> getBounds (lb_, ub_);}

  // cloning method
  virtual exprAux *clone ()
    {return new exprAux (*this);}

  // set lower bound
  void setLB (expression *lb) 
    {if (lb_) delete lb_; lb_ = lb;}

  // set upper bound
  void setUB (expression *ub) 
    {if (ub_) delete ub_; ub_ = ub;}

  // Bound get
  expression *Lb () {return lb_;}
  expression *Ub () {return ub_;}

  // I/O
  void print (std::ostream &out)
    {out << "w_" << varIndex_;}

  // The expression associated with this auxiliary variable
  inline expression *Image () const
    {return image_;}

  // Null function for evaluating the expression
  inline CouNumber operator () () 
    {return (currValue_ = (*image_) ());}

  // Differentiation
  inline expression *differentiate (int index) 
    {return image_ -> differentiate (index);}

  // Dependence on variable set
  inline bool dependsOn (int *indices, int num) 
    {return image_ -> dependsOn (indices, num);}

  // Get a measure of "how linear" the expression is:
  //
  // 0: a constant
  // 1: linear
  // 2: quadratic
  // 3: nonlinear non-quadratic
  inline int Linearity ()
    {return image_ -> Linearity ();}

  // Get lower and upper bound of an expression (if any)
  inline void getBounds (expression *&lb, expression *&ub) 
    {image_ -> getBounds (lb, ub);}

  // generate cuts for expression associated with this auxiliary
  void generateCuts (const OsiSolverInterface &si, 
		     OsiCuts &cs, const CouenneCutGenerator *cg)
    {/*printf ("----------------Generating cut for "); 
    print (std::cout);  printf (" := ");
    image_ -> print (std::cout); printf("\n");
    int j = cs.sizeRowCuts ();*/
    image_ -> generateCuts (this, si, cs, cg);
    /*for (;j < cs.sizeRowCuts ();j++)
      cs.rowCutPtr (j) -> print ();*/}

  // generate equality between *this and *w
  //  void generateCuts (exprAux *w, const OsiSolverInterface &si, 
  //		     OsiCuts &cs, const CouenneCutGenerator *cg);
};

#endif
