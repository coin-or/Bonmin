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
	BonCbc.cpp BonCbc.hpp \
	BonCbcNlpStrategy.cpp BonCbcNlpStrategy.hpp \
	BonCbcNode.cpp BonCbcNode.hpp \
        BonBabInfos.cpp BonBabInfos.hpp \
	BonGuessHeuristic.cpp BonGuessHeuristic.hpp \
	BonDiver.cpp BonDiver.hpp

libbonmin_la_LIBADD = $(libbonmin_la_DEPENDENCIES)
if DEPENDENCY_LINKING
libbonmin_la_LIBADD += $(BONMINLIB_LIBS)
endif

libbonmin_la_DEPENDENCIES = \
	../Algorithms/libbonalgorithms.la \
	../Interfaces/libbonmininterfaces.la \
         Heuristics/libbonheuristics.la

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
libbonminampl_la_LIBADD = $(libbonminampl_la_DEPENDENCIES)
if DEPENDENCY_LINKING
libbonminampl_la_LIBADD += $(BONMINLIB_LIBS) \
	$(BONMINAMPLINTERFACELIB_LIBS)
endif

# List all dependency libraries
libbonminampl_la_DEPENDENCIES = \
	../Algorithms/Ampl/libbonamplsetup.la \
	../Interfaces/Ampl/libbonamplinterface.la 
if DEPENDENCY_LINKING
libbonminampl_la_DEPENDENCIES += libbonmin.la
endif

# This is for libtool
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
	$(BONMINLIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Interfaces

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
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
