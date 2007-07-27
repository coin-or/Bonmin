/*
 * Name:    rootQ.hpp
 * Author:  Pietro Belotti
 * Purpose: find roots of polynomial Q^k(x) (see Liberti and Pantelides, 2003)
 *
 * (C) Pietro Belotti. This file is licensed under the Common Public License (CPL)
 */

#include <CouenneTypes.h>
#include <map>

/// Find roots of polynomial $Q^k(x) = \sum_{i=1}^{2k} i x^{i-1}$. Used
/// in convexification of powers with odd exponent

extern "C" {
  CouNumber rootQ (int k);
}

/// class that stores result of previous calls to rootQ into a map for
/// faster access

class Qroot {

 protected:

  /// Maps an integer k with the root of \f$Q^k(x)\f$.

  static std::map <int, CouNumber> Qmap;

 public:

  /// Empty constructor --we only need the method to work on the static
  /// structure.
  
  Qroot  () {}

  /// Empty destructor.
  ~Qroot () {}

  /// Retrieve root of Q with order = k. If no such computation has
  /// been performed yet, do it here

  inline CouNumber operator () (int k) {

    std::map <int, CouNumber>:: iterator pos;
    CouNumber root;

    k /= 2; // becomes true index

    if ((pos = Qmap.find (k)) == Qmap.end()) {

      std::pair <int, CouNumber> newpair;

      newpair.first  = k;
      newpair.second = (root = rootQ (k));
 
      Qmap.insert (newpair);
    }
    else 
      root = pos -> second;

    return root;
  }
};
