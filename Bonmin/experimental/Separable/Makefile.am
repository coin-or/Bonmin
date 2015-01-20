# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 434 2007-03-28 12:58:54Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

bin_PROGRAMS = sepa

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Apps` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin/Heuristics` \
	$(BONMINLIB_CFLAGS) $(ASL_CFLAGS)

########################################################################
#                                sepa                                  #
########################################################################

sepa_SOURCES = Sepa.cpp \
               SepaSetup.cpp SepaSetup.hpp \
               SepaTMINLP2OsiLP.cpp  SepaTMINLP2OsiLP.hpp \
               SepaHeuristicInnerApproximation.hpp \
               SepaHeuristicInnerApproximation.cpp \
               BonOuterDescription.cpp \
               BonOuterDescription.hpp

sepa_LDADD = \
	../../src/Interfaces/Ampl/libbonamplinterface.la \
	../../src/CbcBonmin/libbonmin.la \
	$(ASL_LIBS) \
	$(BONMINLIB_LIBS)

########################################################################
#                         make distclean stuff                         #
########################################################################

CLEANFILES = 

DISTCLEANFILES = 
