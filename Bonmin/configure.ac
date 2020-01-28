# Copyright (C) 2006, 2008 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Bonmin],[trunk],[bonmin@list.coin-or.org])

AC_COPYRIGHT([
Copyright 2006, 2008 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package BONMIN which is distributed
under the Eclipse Public License.])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/CbcBonmin/BonCbc.hpp)

# Do some initialization work (version numbers, change prefix default, ...)

AC_COIN_INITIALIZE(Bonmin)

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the name of the C++ compiler and appropriate compiler options

AC_COIN_PROG_CXX

# Initialize libtool

AC_COIN_PROG_LIBTOOL

# set RPATH_FLAGS to the compiler link flags required to hardcode location
# of the shared objects (expanded_libdir is set somewhere in configure before)
# (use in examples Makefile)
AC_COIN_RPATH_FLAGS([$expanded_libdir])

#############################################################################
#                           COIN-OR components                              #
#############################################################################

AC_COIN_CHK_PKG(Cbc,[BonminLib])
if test $coin_has_cbc != yes ; then
  AC_MSG_ERROR([Required package Cbc is not available.])
fi

AC_COIN_CHK_PKG(Ipopt,[BonminLib])
if test $coin_has_ipopt != yes ; then
  AC_MSG_ERROR([Required package Ipopt is not available.])
fi

AC_COIN_CHK_PKG(OsiClp,[BonminLib],[osi-clp])
if test $coin_has_osiclp != yes ; then
  AC_MSG_ERROR([Required package OsiClp is not available.])
fi

AC_COIN_CHK_PKG(ASL,[BonminAmplInterfaceLib],[ipoptamplinterface])

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
AC_COIN_CHK_PKG(OsiCpx,[BonminLib],[osi-cplex])

# To get also the directory where cplex.h can be found, we also check for cplex directly.
AC_COIN_CHK_LIB(Cplex,[BonminLib],[-lcplex -lpthread -lm -ldl],[],[],[CPXgetstat])

#############################################################################
#                                FilterSQP                                  #
#############################################################################

AC_ARG_WITH([filtersqp],
            AC_HELP_STRING([--with-filtersqp],
                           [Linker flags for FilterSQP]),
            [filtersqp_libs="$withval"])

if test -n "$filtersqp_libs" ; then
  AC_LANG_PUSH(C)
  ac_save_LIBS="$LIBS"
  LIBS="$filtersqp_libs $LIBS"
  for ac_extra in "no extra underscore" ; do
    for ac_case in "lower case" "upper case" ; do
      for ac_trail in "underscore" "no underscore" ; do
        case $ac_case in
          "lower case")
            ac_name=filtersqp
            ac_confun_name=confun
            ac_objfun_name=objfun
            ac_gradient_name=gradient
            #ac_objgrad_name=objgrad
            ac_hessian_name=hessian
            ;;
          "upper case")
            ac_name=FILTERSQP
            ac_confun_name=CONFUN
            ac_objfun_name=OBJFUN
            ac_gradient_name=GRADIENT
            #ac_objgrad_name=OBJGRAD
            ac_hessian_name=HESSIAN
            ;;
        esac
        if test "$ac_trail" = underscore ; then
          ac_underscore=_
        else
          ac_underscore=
        fi
        AC_MSG_CHECKING([for function $ac_name in $LIBS])
        AC_LINK_IFELSE(
          [AC_LANG_PROGRAM(
            [void $ac_confun_name$ac_underscore(void) {}
             void $ac_objfun_name$ac_underscore(void) {}
             void $ac_gradient_name$ac_underscore(void) {}
             /* void $ac_objgrad_name$ac_underscore(void) {} */
             void $ac_hessian_name$ac_underscore(void) {}
             void $ac_name$ac_underscore();],
            [$ac_name$ac_underscore()])],
          [filtersqp_namemangling="${ac_case}, ${ac_trail}, ${ac_extra}"
           ac_success=yes],
          [ac_success=no])
        AC_MSG_RESULT([$ac_success])
        if test $ac_success = yes ; then
          break 3
        fi
      done
    done
  done
  LIBS=$ac_save_LIBS
  AC_LANG_POP(C)

  if test $ac_success = yes ; then
    AC_COIN_DEFINENAMEMANGLING([FILTERSQP],[$filtersqp_namemangling])
    BONMINLIB_LFLAGS="$filtersqp_libs $BONMINLIB_LFLAGS"
    AC_DEFINE(COIN_HAS_FILTERSQP,1,[Define to 1 if FilterSQP is available])
  else
    AC_MSG_ERROR([Symbol filterSQP not found with FilterSQP flags $filtersqp_libs.])
  fi
fi

AM_CONDITIONAL(COIN_HAS_FILTERSQP, test -n "$filtersqp_libs")

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
                 test/MyTMINLP.hpp:examples/CppExample/MyTMINLP.hpp
                 test/MyBonminC.c:examples/CExample/MyBonmin.c
                 test/MyBonmin.h:examples/CExample/MyBonmin.h
                 test/main.c:examples/CExample/main.c])

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

#AC_COIN_DOXYGEN(CoinUtils Osi Clp Cgl Cbc Ipopt)
AC_COIN_DOXYGEN

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

AC_COIN_FINALIZE_FLAGS([BonminLib BonminAmplInterfaceLib])

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
                 experimental/RobotBonmin/Makefile
                 test/Makefile
                 examples/CppExample/Makefile
                 examples/CExample/Makefile
                 bonmin.pc
                 doxydoc/doxygen.conf])

if test $coin_has_asl = yes ; then
  AC_CONFIG_FILES([bonminamplinterface.pc:src/CbcBonmin/bonminamplinterface.pc.in])
fi

AC_CONFIG_HEADER([src/Interfaces/config.h src/Interfaces/config_bonmin.h])

AC_COIN_FINALIZE
