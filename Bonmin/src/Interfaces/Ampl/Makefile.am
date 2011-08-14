# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          libipoptinterface                           #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libbonamplinterface.la

# List all source files, including headers
libbonamplinterface_la_SOURCES = \
	BonAmplInterface.cpp BonAmplInterface.hpp \
	BonAmplTMINLP.cpp BonAmplTMINLP.hpp \
	BonSolReader.cpp BonSolReader.hpp \
        sos_kludge.cpp

# This is for libtool
libbonamplinterface_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
        -I`$(CYGPATH_W) $(srcdir)/..` \
        -I`$(CYGPATH_W) $(srcdir)/../Ipopt` \
        $(BONMINLIB_CFLAGS) $(ASL_CFLAGS)

if COIN_HAS_FILTERSQP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Filter` 
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
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
