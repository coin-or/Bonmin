/*
 * Name:    depGraph.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for manipulating dependencies between variables
 *
 * (C) Carnegie-Mellon University, 2007. 
 * This file is licensed under the Common Public License (CPL)
 */

#include <cstdlib>
#include "depGraph.hpp"

//#define DEBUG

// Methods of the class DepNode ///////////////////////////////////////////////////////////

/// does this variable depend on variable with index xi?
bool DepNode::depends (int xi, bool recursive, 
		       std::set <DepNode *, compNode> *already_visited) const {

  // check if any node of the forward star has index xi
  for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin (); 
       i != depList_ -> end (); ++i) {

    if (!already_visited || 
	(already_visited -> find (*i) == already_visited -> end ())) {

      if (((*i) -> Index () == xi) || // check current node
	(recursive && 
	 ((*i) -> depends (xi, recursive, already_visited)))) // check deplist recursively
      { 
#ifdef DEBUG
	printf ("%d <- ", (*i) -> Index ()); 
	fflush (stdout);
#endif
	return true;
      } else {
	if (already_visited) {
	  already_visited -> insert (*i);
	  /*printf ("checked [%d]: ", (*i) -> Index ()); 
	  for (std::set <DepNode *, compNode>::iterator j = already_visited -> begin ();
	       j != already_visited -> end (); ++j)
	    printf ("%d ", (*j) -> Index ());
	    printf ("\n");*/
	}
      }
    }
  }

  return false;
}


/// assign numbering to all nodes of graph
void DepNode::createOrder (DepGraph *g) {

  if (order_ != -1) return;

  if (order_ == -2) {

    printf ("detected cycle in creating order, exiting\n");
    exit (-1);
  }

  order_ = -2;

  for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin();
       i != depList_ -> end (); ++i)
    if ((*i) -> Order () == -1)
      (*i) -> createOrder (g);

  if (order_ == -2)
    order_ = g -> Counter () ++;
}


/// debugging procedure
void DepNode::print (int indent, bool descend) const {

  printf ("%d ", index_); 
  if (order_ >= 0) printf ("[%d]", order_); 
  fflush (stdout);

  if (depList_ -> size () > 0) {
    printf ("("); fflush (stdout);

    for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin();
	 i != depList_ -> end (); ++i)
      if (descend)
	(*i) -> print (indent + 1, descend);
      else printf ("%d ", (*i) -> Index ());

    printf (") "); fflush (stdout);
  }
}


// Methods of the class DepGraph ////////////////////////////////////////////////////////////


/// insert new variable if new
void DepGraph::insert (exprVar *var) {

  DepNode *el = new DepNode (var -> Index ());
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el); 

  if (i == vertices_ . end ())
    vertices_.insert (el);
  else delete el;
}


/// insert new auxiliary if new
void DepGraph::insert (exprAux *aux) {

  if (!aux) return;

  DepNode *el = new DepNode (aux -> Index ());
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el); 

  if (i == vertices_ . end ()) {
    vertices_.insert (el);
    aux      -> Image () -> fillDepSet (el -> DepList (), this);
  } else {
    aux -> Image () -> fillDepSet ((*i) -> DepList (), this);
    delete el;
  }
}


/// erase element from graph
void DepGraph::erase (exprVar *var) {

  DepNode *el = new DepNode (var -> Index ());
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el); 

  if (i != vertices_ . end ())
    vertices_.erase (i);
  delete el;
}

/// does w depend on x?
bool DepGraph::depends (int wi, int xi, bool recursive) {

  DepNode *el = new DepNode (wi);
  std::set <DepNode *, compNode>::iterator i = vertices_. find (el);
  delete el;

  std::set <DepNode *, compNode> already_visited;

  if (i != vertices_. end ())               // if such element is in the set
    return (*i) -> depends (xi, recursive, &already_visited); // then search it
  else return false;
}


/// assign numbering to all nodes of graph
void DepGraph::createOrder () {

  for (std::set <DepNode *, compNode>::iterator i = vertices_. begin();
       i != vertices_. end (); ++i)
    (*i) -> createOrder (this);
}


/// debugging procedure
void DepGraph::print (bool descend) {

  printf ("Dependence graph: \n");
  for (std::set <DepNode *, compNode>::iterator i = vertices_. begin();
       i != vertices_. end (); ++i) {
    (*i) -> print (0, descend);
    printf ("\n");
  }
}


/// search for node in vertex set
DepNode *DepGraph::lookup (int index) {

  DepNode *el = new DepNode (index), *ret;
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el);

  ret = (i == vertices_.end ()) ? NULL : (*i);

  delete el;
  return ret;
}
