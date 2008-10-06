# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

#We need to do autothings in the subdirectories for the different interfaces
SUBDIRS = Ipopt
if COIN_HAS_ASL
SUBDIRS += Ampl
endif
if COIN_HAS_FILTERSQP
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
        BonAuxInfos.cpp BonAuxInfos.hpp \
	BonColReader.cpp BonColReader.hpp \
	BonCurvatureEstimator.cpp BonCurvatureEstimator.hpp \
	BonCutStrengthener.cpp BonCutStrengthener.hpp \
	BonStartPointReader.cpp BonStartPointReader.hpp \
	BonOsiTMINLPInterface.cpp BonOsiTMINLPInterface.hpp \
	BonTMINLP2TNLP.cpp BonTMINLP2TNLP.hpp \
	BonTMINLP.cpp BonTMINLP.hpp \
	BonTNLPSolver.hpp BonTNLPSolver.cpp \
	BonTNLP2FPNLP.cpp BonTNLP2FPNLP.hpp \
	BonBranchingTQP.cpp BonBranchingTQP.hpp \
	BonStrongBranchingSolver.cpp BonStrongBranchingSolver.hpp \
        BonRegisteredOption.hpp BonRegisteredOptions.cpp

#sublibraries for various installed solvers
libbonmininterfaces_la_LIBADD = Ipopt/libipoptinterface.la 
#                              Ampl/libamplinterface.la

if COIN_HAS_FILTERSQP
libbonmininterfaces_la_LIBADD += Filter/libfilterinterface.la
endif

libbonmininterfaces_la_DEPENDENCIES = $(libbonmininterfaces_la_LIBADD)

# This is for libtool (on Windows)
libbonmininterfaces_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/Branching` \
	-I`$(CYGPATH_W) $(srcdir)/../CbcBonmin` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglOddHole` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglClique` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglFlowCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglTwomir` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglPreProcess` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(top_srcdir)/src/Algorithms` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/contrib/CGPenalty` \
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

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
     ../../inc/config_bonmin.h \
     BonOsiTMINLPInterface.hpp \
     BonTMINLP2TNLP.hpp \
     BonAuxInfos.hpp \
     BonTMINLP.hpp \
     BonTNLP2FPNLP.hpp \
     BonCutStrengthener.hpp \
     BonTNLPSolver.hpp \
     BonBranchingTQP.hpp \
     BonTypes.hpp \
     BonRegisteredOptions.hpp


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
	$(install_sh_DATA) bla $(DESTDIR)$(includecoindir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/$(ConfigHeader)

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonAuxInfos.cppbak \
	BonAuxInfos.hppbak \
	BonBoundsReader.cppbak \
	BonBoundsReader.hppbak \
	BonBranchingTQP.cppbak \
	BonBranchingTQP.hppbak \
	BonColReader.cppbak \
	BonColReader.hppbak \
	BonCurvatureEstimator.cppbak \
	BonCurvatureEstimator.hppbak \
	BonCutStrengthener.cppbak \
	BonCutStrengthener.hppbak \
	BonExitCodes.hppbak \
	BonOsiTMINLPInterface.cppbak \
	BonOsiTMINLPInterface.hppbak \
	BonRegisteredOptions.cppbak \
	BonRegisteredOptions.hppbak \
	BonStartPointReader.cppbak \
	BonStartPointReader.hppbak \
	BonStrongBranchingSolver.cppbak \
	BonStrongBranchingSolver.hppbak \
	BonTMINLP2TNLP.cppbak \
	BonTMINLP2TNLP.hppbak \
	BonTMINLP.cppbak \
	BonTMINLP.hppbak \
	BonTNLP2FPNLP.cppbak \
	BonTNLP2FPNLP.hppbak \
	BonTNLPSolver.cppbak \
	BonTNLPSolver.hppbak \
	BonTypes.hppbak


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
