# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 83 2006-10-23 16:53:13Z pbonami $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          libipoptinterface                           #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libipoptinterface.la

# List all source files, including headers
libipoptinterface_la_SOURCES = \
	BonIpoptInteriorWarmStarter.cpp  BonIpoptSolver.cpp \
        BonIpoptWarmStart.cpp \
        BonIpoptInteriorWarmStarter.hpp  BonIpoptSolver.hpp \
        BonIpoptWarmStart.hpp

# This is for libtool
libipoptinterface_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/..` $(BONMINLIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	BonIpoptInteriorWarmStarter.hpp \
	BonIpoptSolver.hpp \
	BonIpoptWarmStart.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
         BonIpoptInteriorWarmStarter.cppbak \
         BonIpoptSolver.cppbak \
         BonIpoptWarmStart.cppbak \
         BonIpoptInteriorWarmStarter.hppbak \
         BonIpoptSolver.hppbak \
         BonIpoptWarmStart.hppbak

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
