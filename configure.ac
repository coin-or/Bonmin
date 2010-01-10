# Copyright (C) 2006 International Business Machines.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter            IBM    2006-04-13

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([Bonmin],[1.2.1],[coin-bonmin@lists.coin-or.org])

AC_COPYRIGHT([
Copyright 2006 International Business Machines and others.
All Rights Reserved.
This file is part of the open source package Coin which is distributed
under the Common Public License.])

AC_REVISION(0.99)

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
#             Check which ThirdParty sources are available                  #
#############################################################################

AC_COIN_THIRDPARTY_SUBDIRS([Blas],[ThirdParty/Blas],[daxpy.f])
AC_COIN_THIRDPARTY_SUBDIRS([Lapack],[ThirdParty/Lapack],[LAPACK/SRC/dlarf.f])
AC_COIN_THIRDPARTY_SUBDIRS([HSL],[ThirdParty/HSL],[configure.ac])
AC_COIN_THIRDPARTY_SUBDIRS([ASL],[ThirdParty/ASL],[solvers/asl.h])
AC_COIN_THIRDPARTY_SUBDIRS([FilterSQP],[ThirdParty/FilterSQP],[source/src/filter.f])
AC_COIN_THIRDPARTY_SUBDIRS([Mumps],[ThirdParty/Mumps],[MUMPS/src/dmumps_part1.F])
AC_COIN_THIRDPARTY_SUBDIRS([Metis],[ThirdParty/Metis],[metis-4.0/Lib/balance.c])

#############################################################################
#                  Check which subprojects are there                        #
#############################################################################

AC_COIN_MAIN_SUBDIRS(CoinUtils)
AC_COIN_MAIN_SUBDIRS(Data/Netlib)
AC_COIN_MAIN_SUBDIRS(Data/Sample)
AC_COIN_MAIN_SUBDIRS(Data/miplib3)
AC_COIN_MAIN_SUBDIRS(Clp)
AC_COIN_MAIN_SUBDIRS(Osi)
AC_COIN_MAIN_SUBDIRS(Cgl)
AC_COIN_MAIN_SUBDIRS(Cbc)
AC_COIN_MAIN_SUBDIRS(Bcp)
AC_COIN_MAIN_SUBDIRS(Ipopt)
AC_COIN_MAIN_SUBDIRS(Couenne)
AC_COIN_MAIN_SUBDIRS(Bonmin)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile doxydoc/doxygen.conf])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
