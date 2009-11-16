# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

if COIN_HAS_ASL
  bin_PROGRAMS = bonmin 
endif


COINLIBS = \
	$(CBCOBJDIR)/src/libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(OSIOBJDIR)/src/OsiClp/libOsiClp.la

if COIN_HAS_CPX
  COINLIBS += \
	$(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la \
	$(CPXLIB)
endif


COINLIBS += \
	$(OSIOBJDIR)/src/libOsi.la \
	$(CLPOBJDIR)/src/libClp.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

IPOPTLIBS = \
	$(IPOPTOBJDIR)/src/Interfaces/libipopt.la

LIBS += `cat $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt`
LIBS += `cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`

########################################################################
#                                bonmin                                #
########################################################################

bonmin_SOURCES = BonMin.cpp

bonmin_LDADD = \
	../CbcBonmin/libbonminampl.la \
	../CbcBonmin/libbonmin.la \
	$(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la \
	$(IPOPTLIBS) \
	$(COINLIBS) \
	$(ASLLIB) \
        $(ADDLIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
bonmin_LDFLAGS = -rpath $(libdir)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/Branching` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithms/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/contrib/CGPenalty` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/examples`

if COIN_HAS_CPX
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                              nodeSolver                              #
########################################################################

# BonNodeSolver_SOURCES = BonNodeSolver.cpp
# 
# BonNodeSolver_LDADD = \
# 	../CbcBonmin/libbonminampl.la \
# 	../CbcBonmin/libbonmin.la \
# 	$(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la \
# 	$(IPOPTLIBS) \
# 	$(COINLIBS) \
# 	$(ASLLIB) \
# 	$(ADDLIBS)
# 
# # Finally, the -rpath flag is used by libtool to make sure that the shared
# # library is found (in the lib install directory) when we are using dynamic
# # libraries.
# BonNodeSolver_LDFLAGS = -rpath $(libdir)


########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
             BonMin.cppbak \
             BonNodeSolver.cppbak

CLEANFILES = $(ASTYLE_FILES)

DISTCLEANFILES =

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
