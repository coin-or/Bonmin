# Copyright (C) 2006, 2008 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Bonmin],[1.5.1],[bonmin@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006, 2008 International Business Machines and others.
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
AC_COIN_PROJECTDIR_INIT(Bonmin,9:0:5,9:1:5)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(Bonmin)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# AW: I'm adding the Fortran compiler now for FilterSQP - we should
# make this optional

# Get the name of the Fortran compiler and appropriate compiler options
AC_COIN_PROG_F77

# Find out how to call Fortran from C and determine Fortran runtime libraries
AC_COIN_F77_WRAPPERS

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                           COIN-OR components                              #
#############################################################################

AC_COIN_CHECK_PACKAGE(CoinDepend, [cbc osi-clp ipopt], [BonminLib])
if test $coin_has_coindepend != yes ; then
  AC_MSG_ERROR([Required package Cbc or Ipopt not available.])
fi

AC_COIN_CHECK_PACKAGE(ASL, [ipoptamplinterface], [BonminAmplInterfaceLib])

#AC_COIN_CHECK_PACKAGE(Bcp, [bcp])
#
#AC_ARG_WITH([bonminbcp],
#            AC_HELP_STRING([--with-bonminbcp],
#                           [Compile Bonmin with Bcp-based parallel version]),
#            [use_bcp="$withval"], [use_bcp=no])
#if test "$use_bcp" = yes; then
#  if test $coin_has_bcp != yes ; then
#     AC_MSG_ERROR([You specified --with-bonminbcp, but Bcp is not available])
#  fi
#fi
#AM_CONDITIONAL(COMPILE_BONMINBCP, test "$use_bcp" = yes)

#############################################################################
#                                     CPLEX                                 #
#############################################################################

# Check whether OsiCplex is available
AC_COIN_CHECK_PACKAGE(OsiCpx, [osi-cplex], [BonminLib])

# To get also the directory where cplex.h can be found, we also check for cplex directly.
AC_COIN_CHECK_USER_LIBRARY([Cplex], [CPX],[cplex.h],[CPXgetstat])

#############################################################################
#                                FilterSQP                                  #
#############################################################################

# for backward compatibility, map --with-filtersqp to --with-filtersqp-libs
AC_ARG_WITH([filtersqp],
            AC_HELP_STRING([--with-filtersqp],
                           [specify library for FilterSQP (or BUILD for compilation)]),
            [filtersqp_libs="$withval"])

AC_COIN_CHECK_PACKAGE(FilterSQP, [coinfiltersqp], [BonminLib])

# we want to keep the linking check for FilterSQP, but cannot do this if FilterSQP will be build later (in ThirdParty)
#if test $coin_has_filtersqp = yes ; then
#  AC_MSG_CHECKING([whether FilterSQP library $FILTERSQP_LIBS works])
#  SAVE_LIBS="$LIBS"
#  LIBS="$FILTERSQP_LIBS $FLIBS $LIBS"
#  AC_LANG_PUSH([Fortran 77])
#  AC_TRY_LINK([],
#[      call filterSQP
#      end
#      subroutine confun
#      end
#      subroutine objfun
#      end
#      subroutine gradient
#      end
#      subroutine hessian],
#              [AC_MSG_RESULT([yes])
#               FILTERSQP_LIBS="$FILTERSQP_LIBS $FLIBS"],
#              [AC_MSG_RESULT([no])
#               AC_MSG_ERROR([user supplied FilterSQP library $FILTERSQP_LIBS does not work])])
#  AC_LANG_POP([Fortran 77])
#  LIBS="$SAVE_LIBS"
#fi

#############################################################################
#                           Feasibility Pump                                #
#############################################################################

#if test x"$BUILD_FP" != x; then
#  build_fp=yes
#fi
#AM_CONDITIONAL(BUILD_FP, test x$build_fp = xyes)

#############################################################################
#                                 ASTYLE                                    #
#############################################################################

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
                 src/CbcBonmin/Makefile
                 src/Interfaces/Makefile
                 src/Interfaces/Ipopt/Makefile
                 src/Interfaces/Filter/Makefile
                 src/Interfaces/Ampl/Makefile
                 src/Algorithms/Makefile
                 src/Algorithms/Ampl/Makefile
                 src/Algorithms/Branching/Makefile
                 src/Algorithms/QuadCuts/Makefile
                 src/Algorithms/OaGenerators/Makefile
                 src/CbcBonmin/Heuristics/Makefile
                 examples/OptionDocGen/Makefile
                 experimental/Separable/Makefile
                 test/Makefile
                 examples/CppExample/Makefile
                 bonmin.pc bonmin-uninstalled.pc])

if test $coin_has_asl = yes ; then
  AC_CONFIG_FILES([bonminamplinterface.pc:src/CbcBonmin/bonminamplinterface.pc.in
                   bonminamplinterface-uninstalled.pc:src/CbcBonmin/bonminamplinterface-uninstalled.pc.in
                  ])
fi

AC_CONFIG_HEADER([src/Interfaces/config.h src/Interfaces/config_bonmin.h])

AC_COIN_FINALIZE
