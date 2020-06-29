# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                          libipoptinterface                           #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libbonamplinterface.la

# List all source files, including headers
libbonamplinterface_la_SOURCES = \
	BonAmplInterface.cpp \
	BonAmplTMINLP.cpp \
	BonSolReader.cpp \
	sos_kludge.cpp

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
        -I$(srcdir)/.. \
        -I$(srcdir)/../Ipopt \
        -I$(srcdir)/../Filter \
        $(BONMINAMPLINTERFACELIB_CFLAGS)

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin-or'
includecoindir = $(includedir)/coin-or
includecoin_HEADERS = BonAmplInterface.hpp BonAmplTMINLP.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonAmplTMINLP.cppbak BonAmplTMINLP.hppbak \
	BonAmplInterface.cppbak BonAmplInterface.hppbak 

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
