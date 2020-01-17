# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 433 2007-03-28 05:01:22Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory
noinst_LTLIBRARIES = libbonbranching.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonbranching_la_SOURCES = \
        BonChooseVariable.cpp \
        BonRandomChoice.cpp \
        BonPseudoCosts.cpp \
        BonLpBranchingSolver.cpp \
        BonQpBranchingSolver.cpp

# This is for libtool
libbonbranching_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.
# Currently, CbcBonmin has to be included for BonChooseVariable.hpp, but
# we should try to unmangle the connection to Cbc.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../OaGenerators \
	-I$(srcdir)/../../CbcBonmin \
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
	BonChooseVariable.hpp \
	BonPseudoCosts.hpp \
	BonCurvBranchingSolver.hpp \
	BonLpBranchingSolver.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonChooseVariable.cppbak BonChooseVariable.hppbak \
	BonPseudoCosts.cppbak BonPseudoCosts.hppbak \
	BonQpBranchingSolver.cppbak BonQpBranchingSolver.hppbak \
	BonCurvBranchingSolver.cppbak BonCurvBranchingSolver.hppbak \
	BonLpBranchingSolver.hppbak BonLpBranchingSolver.cppbak

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
