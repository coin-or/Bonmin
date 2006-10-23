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
noinst_LTLIBRARIES = libipoptinterface.la 
# List all source files, including headers
libipoptinterface_la_SOURCES = \
	BonBoundsReader.cpp BonBoundsReader.hpp \
	BonColReader.cpp BonColReader.hpp \
	BonCurvatureEstimator.cpp BonCurvatureEstimator.hpp \
	BonStartPointReader.cpp BonStartPointReader.hpp \
	BonOsiTMINLPInterface.cpp BonOsiTMINLPInterface.hpp \
	BonIpoptInteriorWarmStarter.cpp BonIpoptInteriorWarmStarter.hpp \
	BonIpoptWarmStart.cpp BonIpoptWarmStart.hpp \
	BonTMINLP2TNLP.cpp BonTMINLP2TNLP.hpp \
	BonTMINLP.cpp BonTMINLP.hpp \
        BonIpoptSolver.hpp BonIpoptSolver.cpp \
	BonTNLPSolver.hpp BonTNLPSolver.cpp \
	BonTNLP2FPNLP.cpp BonTNLP2FPNLP.hpp

if COIN_HAS_FSQP
libipoptinterface_la_SOURCES += \
	BonFilterSolver.cpp BonFilterSolver.hpp
endif
# This is for libtool (on Windows)
libipoptinterface_la_LDFLAGS = $(LT_LDFLAGS)

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
	BonIpoptInteriorWarmStarter.hpp \
	BonTMINLP2TNLP.hpp \
	BonTMINLP.hpp \
	BonTNLP2FPNLP.hpp \
        BonTNLPSolver.hpp \
        BonIpoptSolver.hpp \
	BonFilterSolver.hpp

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
	BonIpoptInteriorWarmStarterbak.cpp BonIpoptInteriorWarmStarter.hpp \
	BonIpoptWarmStartbak.cpp BonIpoptWarmStart.hpp \
	BonTMINLP2TNLPbak.cpp BonTMINLP2TNLP.hpp \
	BonTMINLPbak.cpp BonTMINLP.hpp \
        BonIpoptSolverbak.hpp BonIpoptSolver.cpp \
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
