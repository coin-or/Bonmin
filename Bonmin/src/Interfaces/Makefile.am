# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

#We need to do autothings in the subdirectories for the different interfaces
SUBDIRS = Ipopt Ampl

if COIN_HAS_FSQP
SUBDIRS += Filter
endif

########################################################################
#                          libipoptinterface                           #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libbonmininterfaces.la 
# List all source files, including headers
libbonmininterfaces_la_SOURCES = \
	BonBoundsReader.cpp BonBoundsReader.hpp \
	BonColReader.cpp BonColReader.hpp \
	BonCurvatureEstimator.cpp BonCurvatureEstimator.hpp \
	BonStartPointReader.cpp BonStartPointReader.hpp \
	BonOsiTMINLPInterface.cpp BonOsiTMINLPInterface.hpp \
	BonTMINLP2TNLP.cpp BonTMINLP2TNLP.hpp \
	BonTMINLP.cpp BonTMINLP.hpp \
	BonTNLPSolver.hpp BonTNLPSolver.cpp \
	BonTNLP2FPNLP.cpp BonTNLP2FPNLP.hpp

#sublibraries for various installed solvers
libbonmininterfaces_la_LIBADD = Ipopt/libipoptinterface.la \
                              Ampl/libamplinterface.la

if COIN_HAS_FSQP
libbonmininterfaces_la_LIBADD += Filter/libfilterinterface.la
endif

libbonmininterfaces_la_DEPENDENCIES = $(libbonmininterfaces_la_LIBADD)

# This is for libtool (on Windows)
libbonmininterfaces_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm/LinearSolvers` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Apps/AmplSolver` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` 

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

include_HEADERS = \
	BonOsiTMINLPInterface.hpp \
	BonTMINLP2TNLP.hpp \
	BonTMINLP.hpp \
	BonTNLP2FPNLP.hpp \
        BonTNLPSolver.hpp 

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = BonminConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonBoundsReaderbak.cpp BonBoundsReaderbak.hpp \
	BonColReaderbak.cpp BonColReaderbak.hpp \
	BonStartPointReaderbak.cpp BonStartPointReader.hpp \
	BonOsiTMINLPInterfacebak.cpp BonOsiTMINLPInterface.hpp \
	BonTMINLP2TNLPbak.cpp BonTMINLP2TNLP.hpp \
	BonTMINLPbak.cpp BonTMINLP.hpp \
	BonTNLPSolverbak.hpp BonTNLPSolver.cpp \
	BonTNLP2FPNLPbak.cpp BonTNLP2FPNLP.hpp


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
