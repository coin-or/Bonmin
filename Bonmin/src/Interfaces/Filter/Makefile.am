# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

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
	BonBqpdSolver.cpp BonBqpdSolver.hpp \
	BonBqpdWarmStart.cpp BonBqpdWarmStart.hpp \
        BonFilterSolver.cpp BonFilterSolver.hpp \
        BonFilterWarmStart.cpp BonFilterWarmStart.hpp

# This is for libtool (on Windows)
libfilterinterface_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Apps/AmplSolver` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
        -I`$(CYGPATH_W) $(srcdir)/..`

AM_CPPFLAGS += $(ASL_CPPFLAGS)


# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
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
