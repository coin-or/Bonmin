# Copyright (C) 2006, 2008 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Bonmin],[1.3.3],[coin-bonmin@list.coin-or.org])

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
AC_COIN_PROJECTDIR_INIT

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
AC_COIN_HAS_PROJECT(Couenne)

AC_ARG_WITH([bonminbcp],
            AC_HELP_STRING([--with-bonminbcp],
                           [Compile Bonmin with Bcp-based parallel version]),
            [use_bcp="$withval"], [use_bcp=no])
if test "$use_bcp" = yes; then
  if test $m4_tolower(coin_has_$1) = unavailable ||
     test $m4_tolower(coin_has_$1) = skipping; then
     AC_MSG_ERROR([You specified --with-bonminbcp, but Bcp is not available])
  fi
fi
AM_CONDITIONAL(COMPILE_BONMINBCP, test "$use_bcp" = yes)

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available ToDo!
AC_COIN_HAS_USER_LIBRARY([Cplex],[CPX],[cplex.h],[CPXgetstat])
AC_COIN_HAS_GLPK(Glpk)
AC_COIN_HAS_USER_LIBRARY([FortMP],[FMP])  #ToDo! is Fortran?
AC_COIN_HAS_USER_LIBRARY([Mosek],[MSK],[mosek.h],[MSK_openmosek])
AC_COIN_HAS_USER_LIBRARY([Osl],[OSL],[ekk_c_api.h],[ekk_simplex])
AC_COIN_HAS_USER_LIBRARY([Soplex],[SPX],[spxsolver.h]) # ToDo
AC_COIN_HAS_USER_LIBRARY([Xpress],[XPR],[xprs.h],[XPRSsetintcontrol])

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

#####################
# FILTER SQP SOLVER #
#####################

filtersqpobjdir=../ThirdParty/FilterSQP
filtersqpsrcdir=$srcdir/$filtersqpobjdir

AC_MSG_NOTICE([Checking for Filter SQP library])

AC_ARG_WITH([filtersqp],
            AC_HELP_STRING([--with-filtersqp],
                           [specify library for FilterSQP (or BUILD for compilation)]),
            [use_filtersqp=$withval], [use_filtersqp=detect])

if test "$use_filtersqp" = BUILD || test "$use_filtersqp" = detect; then
# Check which files are available
  AC_COIN_CHECK_FILE([$filtersqpsrcdir/source/src/filter.f],
	             [AC_DEFINE(COIN_HAS_FILTERSQP,1,
                                [Define to 1 if FilterSQP is available])
                      use_filtersqp=BUILD],
                     [if test "$use_filtersqp" = BUILD; then
                       AC_MSG_ERROR([Filter SQP files are not available])
                      fi])
fi
if test "$use_filtersqp" = detect; then
  use_filtersqp=no
fi

if test "$use_filtersqp" = BUILD; then
  FILTERSQPLIBADD=`cd $filtersqpobjdir; pwd`/libfiltersqp.la
  coin_need_flibs=yes
elif test "$use_filtersqp" != no; then
  AC_MSG_CHECKING([whether user supplied FILTERSQPLIB=\"$use_filtersqp\" works])
  SAVE_LIBS="$LIBS"
  LIBS="$use_filtersqp $LIBS $FLIBS"
  ADDLIBS="$use_filtersqp $ADDLIBS"
  AC_LANG_PUSH([Fortran 77])
  AC_TRY_LINK([],
[      call filterSQP
      end
      subroutine confun
      end
      subroutine objfun
      end
      subroutine gradient
      end
      subroutine hessian],
              [AC_MSG_RESULT([yes])],
              [AC_MSG_RESULT([no])
               AC_MSG_ERROR([user supplied FilterSQP library \"$use_filtersqp\" does not work])])
  AC_LANG_POP([Fortran 77])
  LIBS="$SAVE_LIBS"
fi

AM_CONDITIONAL([COIN_HAS_FILTERSQP], [test "$use_filtersqp" != no])
AC_SUBST(FILTERSQPLIBADD)

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
                 test/Makefile
                 examples/CppExample/Makefile])

AC_CONFIG_HEADER([inc/config_bonmin.h])

AC_COIN_FINALIZE
