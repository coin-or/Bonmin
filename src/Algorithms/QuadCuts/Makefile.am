# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 433 2007-03-28 05:01:22Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

# Name of the library in this directory
noinst_LTLIBRARIES = libbonquadcuts.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonquadcuts_la_SOURCES = \
                  BonOuterApprox.cpp \
                  BonQuadCut.cpp \
                  BonQuadRow.cpp \
                  BonTMatrix.cpp \
                  BonLinearCutsGenerator.cpp \
                  BonTMINLPLinObj.cpp \
                  BonTMINLP2Quad.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/.. \
	-I$(srcdir)/../Branching \
	-I$(srcdir)/../OaGenerators \
	-I$(srcdir)/../../CbcBonmin \
	-I$(srcdir)/../../CbcBonmin/Heuristics \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../Interfaces/Ipopt \
	-I$(srcdir)/../../Interfaces/Filter \
	$(BONMINLIB_CFLAGS) -DBONMIN_BUILD

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
includecoin_HEADERS = \
            BonArraysHelpers.hpp \
            BonOuterApprox.hpp \
            BonQuadCut.hpp \
            BonQuadRow.hpp \
            BonTMatrix.hpp \
            BonTMINLPLinObj.hpp \
            BonTMINLP2Quad.hpp \
            BonLinearCutsGenerator.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
              BonArraysHelpers.hppbak \
              BonOuterApprox.cppbak \
              BonOuterApprox.hppbak \
              BonQuadCut.cppbak \
              BonQuadCut.hppbak \
              BonQuadRow.cppbak \
              BonQuadRow.hppbak \
              BonTMatrix.cppbak \
              BonTMatrix.hppbak \
              BonTMINLPLinObj.hppbak \
              BonTMINLPLinObj.cppbak \
              BonTMINLP2Quad.cppbak \
              BonTMINLP2Quad.hppbak

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
