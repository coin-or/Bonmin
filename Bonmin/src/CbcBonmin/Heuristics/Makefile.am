# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

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
       BonHeuristicDiveMIPVectorLength.cpp BonHeuristicDiveMIPVectorLength.hpp

# This is for libtool (on Windows)
libbonheuristics_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../../Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces/Filter` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/contrib/CGPenalty` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm/LinearSolvers` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglPreProcess` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc`

# Currently, the following has to be included for BonChooseVariable.hpp, but
# we should try to unmangle the connection to Cbc

AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../../CbcBonmin`

if COIN_HAS_CPX
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc


########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
                     BonLocalSolverBasedHeuristic.hpp \
                     BonFixAndSolveHeuristic.hpp \
                     BonDummyPump.hpp \
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
