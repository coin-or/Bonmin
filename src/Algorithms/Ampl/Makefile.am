# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

# Name of the library in this directory
noinst_LTLIBRARIES = libbonamplsetup.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonamplsetup_la_SOURCES = BonAmplSetup.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/../OaGenerators \
	-I$(srcdir)/../Branching \
	-I$(srcdir)/.. \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../Interfaces/Ipopt \
	-I$(srcdir)/../../Interfaces/Filter \
	-I$(srcdir)/../../Interfaces/Ampl \
	$(BONMINAMPLINTERFACELIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
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
