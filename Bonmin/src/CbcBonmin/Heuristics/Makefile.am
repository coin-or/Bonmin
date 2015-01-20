# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 433 2007-03-28 05:01:22Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory
noinst_LTLIBRARIES = libbonheuristics.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonheuristics_la_SOURCES = \
       BonLocalSolverBasedHeuristic.cpp  BonLocalSolverBasedHeuristic.hpp \
       BonFixAndSolveHeuristic.cpp  BonFixAndSolveHeuristic.hpp \
       BonDummyPump.cpp BonDummyPump.hpp \
       BonPumpForMinlp.cpp BonPumpForMinlp.hpp \
       BonHeuristicRINS.cpp BonHeuristicRINS.hpp \
       BonHeuristicLocalBranching.cpp BonHeuristicLocalBranching.hpp \
       BonHeuristicFPump.cpp BonHeuristicFPump.hpp \
       BonHeuristicDive.cpp BonHeuristicDive.hpp \
       BonHeuristicDiveFractional.cpp BonHeuristicDiveFractional.hpp \
       BonHeuristicDiveVectorLength.cpp BonHeuristicDiveVectorLength.hpp \
       BonHeuristicDiveMIP.cpp BonHeuristicDiveMIP.hpp \
       BonHeuristicDiveMIPFractional.cpp BonHeuristicDiveMIPFractional.hpp \
       BonMilpRounding.cpp BonMilpRounding.hpp \
       BonHeuristicDiveMIPVectorLength.cpp BonHeuristicDiveMIPVectorLength.hpp

# This is for libtool
libbonheuristics_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../../Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces/Filter` \
	$(BONMINLIB_CFLAGS)

# Currently, the following has to be included for BonChooseVariable.hpp, but
# we should try to unmangle the connection to Cbc

#AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../../CbcBonmin`

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Interfaces


########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
                     BonLocalSolverBasedHeuristic.hpp \
                     BonFixAndSolveHeuristic.hpp \
                     BonDummyPump.hpp \
                     BonPumpForMinlp.hpp \
                     BonMilpRounding.hpp \
                     BonHeuristicRINS.hpp \
                     BonHeuristicLocalBranching.hpp \
                     BonHeuristicFPump.hpp \
                     BonHeuristicDive.hpp \
                     BonHeuristicDiveFractional.hpp \
                     BonHeuristicDiveVectorLength.hpp \
                     BonHeuristicDiveMIP.hpp \
                     BonHeuristicDiveMIPFractional.hpp \
                     BonHeuristicDiveMIPVectorLength.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
   BonLocalSolverBasedHeuristic.cppbak  BonLocalSolverBasedHeuristic.hppbak \
   BonFixAndSolveHeuristic.cppbak  BonFixAndSolveHeuristic.hppbak \
   BonDummyPump.cppbak BonDummyPump.hppbak \
   BonHeuristicRINS.cppbak BonHeuristicRINS.hppbak \
   BonHeuristicLocalBranching.cppbak BonHeuristicLocalBranching.hppbak \
   BonHeuristicFPump.cppbak BonHeuristicFPump.hppbak \
   BonHeuristicDive.cppbak BonHeuristicDive.hppbak \
   BonHeuristicDiveFractional.cppbak BonHeuristicDiveFractional.hppbak \
   BonHeuristicDiveVectorLength.cppbak BonHeuristicDiveVectorLength.hppbak \
   BonHeuristicDiveMIP.cppbak BonHeuristicDiveMIP.hppbak \
   BonHeuristicDiveMIPFractional.cppbak BonHeuristicDiveMIPFractional.hppbak \
   BonHeuristicDiveMIPVectorLength.cppbak BonHeuristicDiveMIPVectorLength.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@
