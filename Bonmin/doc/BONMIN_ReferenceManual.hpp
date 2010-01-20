
// Main bonmin page

/*! \mainpage bonmin reference manual
<p>
bonmin is an open source code for solving general MINLP (Mixed Integer
Non-Linear Programming) problems.
It is distributed on <A HREF='http://www.coin-or.org'> COIN-OR </A> under the
<A HREF='http://www.opensource.org/licenses/cpl.php'> CPL (Common Public
License) </A>. bonmin is a C++ code</p>

<p>
You can use the link at the top of this page to access the documentation
of the different elements of the code.
*/

 In addition you will find here:
<ul>
<li> <A HREF='../Bonmin_UserManual/'> A link </A> to the html version of
bonmin Users manual,
 <li> \subpage Code_structure "Code Structure" briefly presents the code and how it is interfaced to
 <A HREF='http://www.coin-or.org/Cbc'> Cbc </A> and 
 <A HREF='http://projects.coin-or.org/Ipopt'> Ipopt </A> 
 <li> \subpage Example "Example" presents the C++ example provided and explains
 how to interface bonmin directly with a C++ code,
 <li> \subpage Extra "Debugging" presents some functionality of bonmin which
 may be usefull for debugging.
 <li> \subpage Authors, a list of the person who have contributed to this code
 <li> <A HREF='
</ul>
</p>
 */

/*! \page Authors
  Pierre Bonami (Carnegie Mellon University, supported by IBM)<br>
  Carl D. Laird (Carnegie Mellon University)<br>
  Andreas Waechter, project leader (IBM)<br>
  <br>
  have contributed to this code.<br>
  <br>
  We would also like to acknoweldge  <br>
  L.T. Biegler (Carnegie Mellon University)  <br>
  A.R. Conn (IBM)  <br>
  G. Cornuejols (LIF Marseille, Carnegie Mellon University)  <br>
  I.E. Grossmann (Carnegie Mellon University)  <br>
  J. Lee (IBM)  <br>
  A. Lodi (IBM, Universita di Bologna)  <br>
  F. Margot (Carnegie Mellon University)  <br>
  N. Sawaya (Carnegie Mellon University)  <br>
  for their advices, help and other contributions.
*/

// Code structure
/*! \page Code_structure Brief description of code structure
 */


//Debuging stuff
/*!  * \page Extra stuff usefull for debugging */


//Example
/*! \page Example The C++ example
 * \section Coding the NLP for Ipopt
 * \section Transforming the NLP into a MINLP
 */




//Directory layout.

/*! \dir Apps
 * contains the sources to build the bonmin executable program. */
/*! \dir CbcModelForHybrid
* contains version of some of the Cbc files which have been modified to be able to handle MINLP's.*/
/*! \dir Doc
* contains this documentation and the user's manual.*/
/*! \dir IpoptInterface
* contains the source for an 
* <A HREF='http://www.coin-or.org/projects.html#OSI'>  OsiSolverInterface</A> 
* to <A HREF='http://projects.coin-or.org/Ipopt'> Ipopt </A> 
*(OsiInterface's are made for linear problems so it is not truly an OsiInterface but it implements the sufficient for doing a branch-and-bound) */


/*! \dir OaInterface
 * Contains the necessary elements for interfacing outer-approximation inside 
 * <A HREF='http://www.coin-or.org/Cbc'> Cbc </A> using an IpoptInterface-type 
 * interface.*/
