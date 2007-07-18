/*
 * Name:    depGraph.hpp
 * Author:  Pietro Belotti
 * Purpose: class for manipulating dependencies between variables
 *
 * This file is licensed under the Common Public License (CPL)
 */

#ifndef DEPGRAPH_H
#define DEPGRAPH_H

#include <vector>
#include <set>

#include <CouenneTypes.h>
#include <expression.hpp>
#include <exprAux.hpp>
#include <CouenneProblemElem.hpp>


/// structure for comparing nodes in the dependence graph
struct compNode {
  inline bool operator () (const DepNode *n0, const DepNode *n1) const;
};


/// vertex of a dependence graph. Contains variable and its forward
/// star (all variables it depends on)

class DepNode {

protected:

  /// index of variable associated with node
  int index_;

  /// index nodes on which this one depends (forward star in
  /// dependence graph)
  std::set <DepNode *, compNode> *depList_;

  /// order in which this variable should be updated, evaluated, etc.
  int order_;

public:

  /// empty constructor
  /*DepNode ():
    index_   (-1),
    depList_ (NULL),
    order_   (-1) {printf ("=============== called null destructor\n");}*/

  /// fictitious constructor: only fill in index (such object is used
  /// in find() and then discarded)
  DepNode  (int ind):
    index_   (ind),
    depList_ (new std::set <DepNode *, compNode>),
    order_   (-1) {}

  /// destructor
  ~DepNode () 
  {if (depList_) delete depList_;}

  /// return index of this variable
  inline int Index () const 
  {return index_;}

  /// return index of this variable
  inline int Order () const 
  {return order_;}

  /// return all variables it depends on
  inline std::set <DepNode *, compNode> *DepList () const 
  {return depList_;}

  /// does this variable depend on variable with index xi?
  bool depends (int , bool = false) const;

  /// assign numbering to all nodes of graph
  void createOrder (DepGraph *);

  /// debugging procedure
  void print (int = 0) const;
};


/// structure for comparing nodes

bool compNode::operator () (const DepNode *n0, const DepNode *n1) const
{return (n0 -> Index () < n1 -> Index ());}


/// Dependence graph. Shows dependence of auxiliary variable on other
/// (auxiliary and/or original) variables

class DepGraph {

protected:

  /// set of variable nodes
  std::set <DepNode *, compNode> vertices_;

  /// counter to assign numbering to all nodes
  int counter_;

public:

  /// constructor
  DepGraph  (): counter_ (0) {}

  /// destructor
  ~DepGraph () {
    for (std::set <DepNode *, compNode>::iterator i = vertices_.begin ();
	 i != vertices_.end (); i++) 
      delete (*i);
  }

  /// node index counter
  inline int &Counter () 
  {return counter_;}

  /// insert new variable if new
  void insert (exprVar *);

  /// insert new auxiliary if new
  void insert (exprAux *);

  /// does w depend on x?
  bool depends (int, int, bool = false);

  /// assign numbering to all nodes of graph
  void createOrder ();

  /// debugging procedure
  void print ();

  /// search for node in vertex set
  DepNode *lookup (int index);
};

#endif
