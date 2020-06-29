# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                          libipoptinterface                           #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libipoptinterface.la

# List all source files, including headers
libipoptinterface_la_SOURCES = \
	BonIpoptInteriorWarmStarter.cpp \
	BonIpoptSolver.cpp \
	BonIpoptWarmStart.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/.. $(BONMINLIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
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
