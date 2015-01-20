# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 433 2007-03-28 05:01:22Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory
noinst_LTLIBRARIES = libbonamplsetup.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonamplsetup_la_SOURCES = BonAmplSetup.cpp BonAmplSetup.hpp 

# This is for libtool
libbonamplsetup_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../Branching` \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces/Filter` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces/Ampl` \
	$(BONMINLIB_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Interfaces

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = BonAmplSetup.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonAmplSetup.cppbak BonAmplSetup.hppbak

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
