/*
 * Name:    depGraph.cpp
 * Author:  Pietro Belotti
 * Purpose: methods for manipulating dependencies between variables
 *
 * This file is licensed under the Common Public License (CPL)
 */

#include <depGraph.hpp>


// Methods of the class DepNode ///////////////////////////////////////////////////////////


/// does this variable depend on variable with index xi?
bool DepNode::depends (int xi, bool recursive) const {

  // check if any node of the forward star has index xi
  for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin (); 
       i != depList_ -> end (); i++) {

    if ((*i) -> Index () == xi) 
      return true;

    if (recursive) {

      // if recursive, propagate the search on elements of the
      // forward star

      std::set <DepNode *, compNode> *gen2 = (*i) -> DepList ();

      for (std::set <DepNode *, compNode>::iterator j = gen2 -> begin (); 
	   j != gen2 -> end (); j++)
	if ((*j) -> Index () == xi) 
	  return true;
    }
  }
  return false;
}


/// assign numbering to all nodes of graph
void DepNode::createOrder (DepGraph *g) {

  if (order_ != -1) return;

  //printf ("[%d ", index_);

  for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin();
       i != depList_ -> end (); i++)
    if ((*i) -> Order () == -1)
      (*i) -> createOrder (g);

  if (order_ == -1) 
    order_ = g -> Counter () ++;
  //printf ("->%d] ", order_);
}


/// debugging procedure
void DepNode::print (int indent) const {

  printf ("%d ", index_); 
  if (order_ >= 0) printf ("[%d]", order_); 
  fflush (stdout);

  if (depList_ -> size () > 0) {
    printf ("("); fflush (stdout);

    for (std::set <DepNode *, compNode>::iterator i = depList_ -> begin();
	 i != depList_ -> end (); i++)
      (*i) -> print (indent + 1);

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
  } else aux -> Image () -> fillDepSet ((*i) -> DepList (), this);
}


/// does w depend on x?
bool DepGraph::depends (int wi, int xi, bool recursive) {

  DepNode *el = new DepNode (wi);
  std::set <DepNode *, compNode>::iterator i = vertices_. find (el);
  delete el;

  if (i != vertices_. end ())               // if such element is in the set
    return (*i) -> depends (xi, recursive); // then search it
  else return false;
}


/// assign numbering to all nodes of graph
void DepGraph::createOrder () {

  for (std::set <DepNode *, compNode>::iterator i = vertices_. begin();
       i != vertices_. end (); i++) {
    //printf ("creating order, vertex %d (%d) ", 
    //	    (*i) -> Index (), (*i) -> Order ());
    (*i) -> createOrder (this);
    //printf ("---> %d\n", (*i) -> Order ());
  }
}


/// debugging procedure
void DepGraph::print () {

  printf ("------------------------------ dependence graph\n");
  for (std::set <DepNode *, compNode>::iterator i = vertices_. begin();
       i != vertices_. end (); i++) {
    (*i) -> print ();
    printf ("\n");
  }
  printf ("------------------------------\n");
}


/// search for node in vertex set
DepNode *DepGraph::lookup (int index) {

  DepNode *el = new DepNode (index), *ret;
  std::set <DepNode *, compNode>::iterator i = vertices_ . find (el);

  ret = (i == vertices_.end ()) ? NULL : (*i);

  delete el;
  return ret;
}
