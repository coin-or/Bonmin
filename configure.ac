# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Bonmin],[trunk],[bonmin@lists.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Eclipse Public License.])

AC_REVISION(1.3)

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(configure.ac)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#              Create the libtool script for the subdirectories             #
#############################################################################

AC_COIN_CREATE_LIBTOOL

#############################################################################
#                  Check which subprojects are there                        #
#############################################################################

AC_COIN_MAIN_PACKAGEDIR(Blas,     ThirdParty,[daxpy.f])
AC_COIN_MAIN_PACKAGEDIR(Lapack,   ThirdParty,[LAPACK/SRC/dlarf.f])
AC_COIN_MAIN_PACKAGEDIR(ASL,      ThirdParty,[solvers/asl.h])
AC_COIN_MAIN_PACKAGEDIR(FilterSQP,ThirdParty,[source/src/filter.f])
AC_COIN_MAIN_PACKAGEDIR(Metis,    ThirdParty,[metis-4.0/Lib/balance.c])
AC_COIN_MAIN_PACKAGEDIR(Mumps,    ThirdParty,[MUMPS/src/dmumps_part1.F])
AC_COIN_MAIN_PACKAGEDIR(HSL,      ThirdParty,[coinhsl/common/deps.f])
AC_COIN_MAIN_PACKAGEDIR(Sample,   Data)
AC_COIN_MAIN_PACKAGEDIR(CoinUtils)
AC_COIN_MAIN_PACKAGEDIR(Osi)
AC_COIN_MAIN_PACKAGEDIR(Clp)
AC_COIN_MAIN_PACKAGEDIR(Cgl)
AC_COIN_MAIN_PACKAGEDIR(Cbc)
AC_COIN_MAIN_PACKAGEDIR(Bcp)
AC_COIN_MAIN_PACKAGEDIR(Ipopt)
AC_COIN_MAIN_PACKAGEDIR(Bonmin)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile doxydoc/doxygen.conf])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
