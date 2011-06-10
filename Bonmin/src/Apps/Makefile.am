# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

if COIN_HAS_ASL
  bin_PROGRAMS = bonmin
  noinst_PROGRAMS = BonNodeSolver
endif

########################################################################
#                                bonmin                                #
########################################################################

bonmin_SOURCES = BonMin.cpp

bonmin_LDADD = \
	../CbcBonmin/libbonminampl.la \
	../CbcBonmin/libbonmin.la \
	$(BONMINLIB_LIBS) \
	$(ASL_LIBS)

bonmin_DEPENDENCIES = \
	../CbcBonmin/libbonminampl.la \
	../CbcBonmin/libbonmin.la \
	$(BONMINLIB_DEPENDENCIES) \
	$(ASL_DEPENDENCIES)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = $(BONMINLIB_CFLAGS) \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/Branching` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../CbcBonmin`

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Interfaces

########################################################################
#                              nodeSolver                              #
########################################################################

BonNodeSolver_SOURCES = BonNodeSolver.cpp

BonNodeSolver_LDADD = \
	../CbcBonmin/libbonminampl.la \
	../CbcBonmin/libbonmin.la \
	$(BONMINLIB_LIBS) \
	$(ASL_LIBS)

BonNodeSolver_DEPENDENCIES = \
	../CbcBonmin/libbonminampl.la \
	../CbcBonmin/libbonmin.la \
	$(BONMINLIB_DEPENDENCIES) \
	$(ASL_DEPENDENCIES)

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
             BonMin.cppbak \
             BonNodeSolver.cppbak

CLEANFILES = $(ASTYLE_FILES)

DISTCLEANFILES =

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
