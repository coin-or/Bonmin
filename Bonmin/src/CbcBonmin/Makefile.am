# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign


SUBDIRS = Heuristics
# Name of the library in this directory
lib_LTLIBRARIES = libbonmin.la
if COIN_HAS_ASL
lib_LTLIBRARIES += libbonminampl.la
endif
########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonmin_la_SOURCES = \
	BonCbc.cpp BonCbc.hpp \
	BonCbcNlpStrategy.cpp BonCbcNlpStrategy.hpp \
	BonCbcNode.cpp BonCbcNode.hpp \
        BonBabInfos.cpp BonBabInfos.hpp \
	BonGuessHeuristic.cpp BonGuessHeuristic.hpp \
	BonDiver.cpp BonDiver.hpp


#nodist_libbonmin_la_SOURCES = CbcCompareUser.cpp CbcBranchUser.cpp

libbonmin_la_LIBADD = \
	../Algorithms/libbonalgorithms.la \
	../Interfaces/libbonmininterfaces.la \
         Heuristics/libbonheuristics.la \
	$(FILTERSQPLIBADD)

libbonmin_la_DEPENDENCIES = $(libbonmin_la_LIBADD)

# This is for libtool (on Windows)
libbonmin_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                         amplbonmin library                           #
########################################################################

# List all source files, including headers
libbonminampl_la_SOURCES = 

libbonminampl_la_LIBADD = \
	../Algorithms/Ampl/libbonamplsetup.la \
	../Interfaces/Ampl/libbonamplinterface.la 

libbonminampl_la_DEPENDENCIES = $(libbonminampl_la_LIBADD)

# This is for libtool (on Windows)
libbonminampl_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/Branching` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/QuadCuts` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Filter` \
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

if COIN_HAS_CPX
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# The following is a trick to get the source code of Cbc/Samples correctly in
# here, and also to get a link to libCbc into this directory (so that it
# is made)
#CbcCompareUser.cpp: $(CBCSRCDIR)/examples/CbcCompareUser.cpp
#rm -f $@
#$(LN_S) $? $@

#CbcBranchUser.cpp: $(CBCSRCDIR)/examples/CbcBranchUser.cpp
#rm -f $@
#$(LN_S) $? $@

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
     BonCbcNlpStrategy.hpp \
     BonCbc.hpp \
     BonBabInfos.hpp \
     BonCbcNode.hpp \
     BonDiver.hpp \
     BonGuessHeuristic.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonCbc.cppbak BonCbc.hppbak \
	BonCbcNlpStrategy.cppbak BonCbcNlpStrategy.hppbak \
	BonCbcNode.cppbak BonCbcNode.hppbak \
        BonBabInfos.cppbak BonBabInfos.hppbak \
	BonGuessHeuristic.cppbak BonGuessHeuristic.hppbak \
	BonDiver.cppbak BonDiver.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_FILES)

DISTCLEANFILES = \
	CbcCompareUser.cpp CbcBranchUser.cpp

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
