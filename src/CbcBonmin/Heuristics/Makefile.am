# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

# Name of the library in this directory
noinst_LTLIBRARIES = libbonheuristics.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonheuristics_la_SOURCES = \
       BonLocalSolverBasedHeuristic.cpp \
       BonFixAndSolveHeuristic.cpp \
       BonDummyPump.cpp \
       BonPumpForMinlp.cpp \
       BonHeuristicRINS.cpp \
       BonHeuristicLocalBranching.cpp \
       BonHeuristicFPump.cpp \
       BonHeuristicDive.cpp \
       BonHeuristicDiveFractional.cpp \
       BonHeuristicDiveVectorLength.cpp \
       BonHeuristicDiveMIP.cpp \
       BonHeuristicDiveMIPFractional.cpp \
       BonMilpRounding.cpp \
       BonHeuristicDiveMIPVectorLength.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../../Algorithms \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../Interfaces/Ipopt \
	-I$(srcdir)/../../Interfaces/Filter \
	$(BONMINLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
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
