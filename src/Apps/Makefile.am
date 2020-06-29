# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

# What is done in this directory:
if COIN_HAS_ASL
  bin_PROGRAMS = bonmin
  noinst_PROGRAMS = BonNodeSolver
endif

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/../Algorithms \
	-I$(srcdir)/../Algorithms/Branching \
	-I$(srcdir)/../Algorithms/OaGenerators \
	-I$(srcdir)/../Algorithms/Ampl \
	-I$(srcdir)/../Interfaces \
	-I$(srcdir)/../Interfaces/Ipopt \
	-I$(srcdir)/../Interfaces/Ampl \
	-I$(srcdir)/../CbcBonmin \
	$(BONMINAMPLINTERFACELIB_CFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                                bonmin                                #
########################################################################

bonmin_SOURCES = BonMin.cpp

bonmin_LDADD = \
	../CbcBonmin/libbonminampl.la

########################################################################
#                              nodeSolver                              #
########################################################################

BonNodeSolver_SOURCES = BonNodeSolver.cpp

BonNodeSolver_LDADD = \
	../CbcBonmin/libbonminampl.la

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
             BonMin.cppbak \
             BonNodeSolver.cppbak

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
