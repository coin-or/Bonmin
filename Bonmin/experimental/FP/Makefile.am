# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

bin_PROGRAMS = oa_fp

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = $(BONMINLIB_CFLAGS) $(ASL_CFLAGS) \
	-I`$(CYGPATH_W) $(srcdir)/../../src/OaAlgorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin`

# OAFeasibilityPump.cpp needs to include cplex.h
if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
endif

#TODO: get CbcCompareUser.cpp and CbcBranchUser.cpp from Cbc/examples
#TODO: get libOsiCbc via check in configure

########################################################################
#                                 oa_fp                                #
########################################################################

oa_fp_SOURCES = OAFeasibilityPump.cpp FP.cpp FP.hpp

oa_fp_LDADD = \
	../../src/Interfaces/Ampl/libbonamplinterface.la \
	../../src/CbcBonmin/libbonmin.la \
	$(ASL_LIBS) \
	$(BONMINLIB_LIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
oa_fp_LDFLAGS = -rpath $(libdir)

########################################################################
#                         make distclean stuff                         #
########################################################################

CLEANFILES = 

DISTCLEANFILES = 
