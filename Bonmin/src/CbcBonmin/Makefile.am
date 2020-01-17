# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

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
	BonCbc.cpp \
	BonCbcNlpStrategy.cpp \
	BonCbcNode.cpp \
	BonBabInfos.cpp \
	BonGuessHeuristic.cpp \
	BonDiver.cpp

libbonmin_la_LIBADD = \
	../Algorithms/libbonalgorithms.la \
	../Interfaces/libbonmininterfaces.la \
	Heuristics/libbonheuristics.la \
	$(BONMINLIB_LFLAGS)

# This is for libtool
libbonmin_la_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                         amplbonmin library                           #
########################################################################

# List all source files, including headers
libbonminampl_la_SOURCES =

# Dummy C++ source to cause C++ linking
nodist_EXTRA_libbonminampl_la_SOURCES = dummy.cpp 

# List all additionally required libraries
libbonminampl_la_LIBADD = \
	../Algorithms/Ampl/libbonamplsetup.la \
	../Interfaces/Ampl/libbonamplinterface.la \
	libbonmin.la \
	$(BONMINAMPLINTERFACELIB_LFLAGS)

# This is for libtool
libbonminampl_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/../Algorithms \
	-I$(srcdir)/../Algorithms/OaGenerators \
	-I$(srcdir)/../Algorithms/Branching \
	-I$(srcdir)/../Algorithms/QuadCuts \
	-I$(srcdir)/../Interfaces \
	-I$(srcdir)/../Interfaces/Ipopt \
	-I$(srcdir)/../Interfaces/Ampl \
	-I$(srcdir)/../Interfaces/Filter \
	$(BONMINLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
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
