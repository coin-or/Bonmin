# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Bonmin],[0.9],[coin-bonmin@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package BONMIN which is distributed
under the Common Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/CbcBonmin/BonCbc.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE

# on blue gene, we need to define -DBLUE_GENE_XLC for the COIN headers
if test x"$enable_bg" = xyes; then
  ADD_CXXFLAGS="$ADD_CXXFLAGS -DBLUE_GENE_XLC"
fi

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                              COIN components                              #
#############################################################################

AC_COIN_HAS_PROJECT(Cbc)
AC_COIN_HAS_PROJECT(Cgl)
AC_COIN_HAS_PROJECT(Clp)
AC_COIN_HAS_PROJECT(CoinUtils)
AC_COIN_HAS_PROJECT(DyLP)
AC_COIN_HAS_PROJECT(Ipopt)
AC_COIN_HAS_PROJECT(Osi)
AC_COIN_HAS_PROJECT(Sym)
AC_COIN_HAS_PROJECT(Vol)
AC_COIN_HAS_PROJECT(Bcp)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_USER_LIBRARY([Glpk],[GLPK],[glpk.h],[glp_lpx_simplex])
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosekdl.h],[MSK_openmosek])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])
AC_COIN_HAS_USER_LIBRARY([filterSqp],[FSQP]) #ToDo

#############################################################################
#                       Other third party software                          #
#############################################################################

# Ampl Solver library
AC_COIN_HAS_ASL

# Some internal variables
if test x"$BUILD_FP" != x; then
  build_fp=yes
fi
AM_CONDITIONAL(BUILD_FP, test x$build_fp = xyes)

##########
# ASTYLE #
##########
ASTYLE=astyle
AC_SUBST(ASTYLE)
ASTYLEFLAGS="--mode=c --indent=spaces=2 --indent-cases --indent-namespaces --min-conditional-indent=1 --brackets=linux --brackets=break-closing-headers --max-instatement-indent=2"
AC_SUBST(ASTYLEFLAGS)    

#############################################################################
#                            Stuff for Examples                             #
#############################################################################

AC_COIN_VPATH_LINK(test/bonmin.opt)
AC_COIN_VPATH_LINK(test/mytoy.nl)
AC_CONFIG_LINKS([test/MyBonmin.cpp:examples/CppExample/MyBonmin.cpp
                 test/MyTMINLP.cpp:examples/CppExample/MyTMINLP.cpp
                 test/MyTMINLP.hpp:examples/CppExample/MyTMINLP.hpp])

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

AC_CONFIG_FILES([Makefile
                 src/Apps/Makefile
                 experimental/Bcp/Makefile
                 src/CbcBonmin/Makefile
                 src/Interfaces/Makefile
                 src/Interfaces/Ipopt/Makefile
                 src/Interfaces/Filter/Makefile
                 src/Interfaces/Ampl/Makefile
                 src/OaAlgorithms/Makefile
                 experimental/FP/Makefile
                 test/Makefile
                 examples/CppExample/Makefile])

AC_CONFIG_HEADER([inc/config_bonmin.h])

AC_COIN_FINALIZE
