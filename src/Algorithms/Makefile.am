# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 433 2007-03-28 05:01:22Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

SUBDIRS = \
       OaGenerators \
       Branching \
       QuadCuts

if COIN_HAS_ASL
  SUBDIRS += Ampl
endif

# Name of the library in this directory
noinst_LTLIBRARIES = libbonalgorithms.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonalgorithms_la_SOURCES = \
         BonBabSetupBase.cpp \
         BonBonminSetup.cpp \
         BonSubMipSolver.cpp \
         BonCbcLpStrategy.cpp \
         BonSolverHelp.cpp

libbonalgorithms_la_LIBADD = Branching/libbonbranching.la
libbonalgorithms_la_LIBADD += OaGenerators/libbonoagenerators.la 
libbonalgorithms_la_LIBADD += QuadCuts/libbonquadcuts.la 

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/Branching \
	-I$(srcdir)/QuadCuts \
	-I$(srcdir)/OaGenerators \
	-I$(srcdir)/../CbcBonmin \
	-I$(srcdir)/../CbcBonmin/Heuristics \
	-I$(srcdir)/../Interfaces \
	-I$(srcdir)/../Interfaces/Ipopt \
	-I$(srcdir)/../Interfaces/Filter \
	$(BONMINLIB_CFLAGS) -DBONMIN_BUILD

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
includecoin_HEADERS = \
          BonBabSetupBase.hpp \
          BonSubMipSolver.hpp \
          BonCbcLpStrategy.hpp \
          BonBonminSetup.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
               BonBabSetupBase.cppbak \
               BonBabSetupBase.hppbak \
               BonBonminSetup.cppbak \
               BonCbcLpStrategy.cppbak BonCbcLpStrategy.hppbak \
               BonBonminSetup.hppbak

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
