# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 83 2006-10-23 16:53:13Z pbonami $

# Author:  Andreas Waechter           IBM    2006-04-13

# Source: Bonmin/src/IpoptInterface/Makefile.am
# Modified for filterInterface
#          P. Bonami                  IBM    2006-10-02

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          libipoptinterface                           #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libfilterinterface.la

# List all source files, including headers
libfilterinterface_la_SOURCES = \
	BonBqpdSolver.cpp \
	BonBqpdWarmStart.cpp \
	BonFilterSolver.cpp \
	BonFilterWarmStart.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = -I$(srcdir)/.. $(BONMINLIB_CFLAGS) -DBONMIN_BUILD

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
includecoin_HEADERS = \
                  BonFilterSolver.hpp \
                  BonBqpdSolver.hpp \
                  BonFilterWarmStart.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
         BonFilterSolver.cppbak \
         BonFilterSolver.hppbak \
         BonFilterWarmStart.cppbak \
         BonFilterWarmStart.hppbak \
         BonBqpdSolver.hppbak \
         BonBqpdSolver.cppbak \
         BonBqpdWarmStart.hppbak \
         BonBqpdWarmStart.cppbak

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
