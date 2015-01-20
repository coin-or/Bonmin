# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 204 2006-12-28 18:57:46Z pbonami $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                            liboainterface                            #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libbonoagenerators.la

# List all source files, including headers
libbonoagenerators_la_SOURCES = \
	BonDummyHeuristic.cpp BonDummyHeuristic.hpp \
	BonOaNlpOptim.cpp BonOaNlpOptim.hpp \
	BonOACutGenerator2.cpp BonOACutGenerator2.hpp \
        BonOaFeasChecker.cpp BonOaFeasChecker.hpp \
	BonOaDecBase.cpp BonOaDecBase.hpp \
        BonEcpCuts.cpp BonEcpCuts.hpp \
        BonFpForMinlp.cpp BonFpForMinlp.hpp \
	BonOAMessages.cpp BonOAMessages.hpp

# This is for libtool
libbonoagenerators_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../CbcBonmin` \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	$(BONMINLIB_CFLAGS)

# BonOaDecBase.cpp needs to include cplex.h
if COIN_HAS_CPX
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Interfaces

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
        BonDummyHeuristic.hpp \
	BonOaNlpOptim.hpp \
	BonOACutGenerator2.hpp \
        BonOaFeasChecker.hpp \
	BonOaDecBase.hpp \
        BonEcpCuts.hpp \
	BonOAMessages.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonDummyHeuristic.cppbak BonDummyHeuristic.hppbak \
	BonOaNlpOptim.cppbak BonOaNlpOptim.hppbak \
	BonOACutGenerator2.cppbak BonOACutGenerator2.hppbak \
        BonOaFeasChecker.cppbak BonOaFeasChecker.hppbak \
	BonOaDecBase.cppbak BonOaDecBase.hppbak \
        BonEcpCuts.cppbak BonEcpCuts.hppbak \
	BonOAMessages.cppbak BonOAMessages.hppbak

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
