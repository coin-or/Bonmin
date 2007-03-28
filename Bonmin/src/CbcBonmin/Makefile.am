# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# Name of the library in this directory
lib_LTLIBRARIES = libbonmin.la

########################################################################
#                          bonmin library                              #
########################################################################

# List all source files, including headers
libbonmin_la_SOURCES = \
	BonCbc.cpp BonCbc.hpp \
	BonCbcLpStrategy.cpp BonCbcLpStrategy.hpp \
	BonCbcNlpStrategy.cpp BonCbcNlpStrategy.hpp \
	BonCbcNode.cpp BonCbcNode.hpp \
	BonCbcParam.cpp BonCbcParam.hpp \
	BonChooseVariable.cpp BonChooseVariable.hpp \
	BonQPStrongBranching.cpp BonQPStrongBranching.hpp \
	BonCurvBranching.cpp BonCurvBranching.hpp \
        BonLpStrongBranching.cpp BonLpStrongBranching.hpp

nodist_libbonmin_la_SOURCES = CbcCompareUser.cpp CbcBranchUser.cpp

libbonmin_la_LIBADD = \
	../OaAlgorithms/liboainterface.la \
	../Interfaces/libbonmininterfaces.la \
	$(FILTERSQPLIBADD)

libbonmin_la_DEPENDENCIES = $(libbonmin_la_LIBADD)

# This is for libtool (on Windows)
libbonmin_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../OaAlgorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces/Filter` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm/LinearSolvers` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
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
	-I`$(CYGPATH_W) $(CBCSRCDIR)/examples`

if COIN_HAS_CPX
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# The following is a trick to get the source code of Cbc/Samples correctly in
# here, and also to get a link to libCbc into this directory (so that it
# is made)
CbcCompareUser.cpp: $(CBCSRCDIR)/examples/CbcCompareUser.cpp
	rm -f $@
	$(LN_S) $? $@

CbcBranchUser.cpp: $(CBCSRCDIR)/examples//CbcBranchUser.cpp
	rm -f $@
	$(LN_S) $? $@

########################################################################
#                Headers that need to be installed                     #
########################################################################

include_HEADERS = \
	BonCbc.hpp \
	BonCbcParam.hpp

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	BonCbcParam.cppbak BonCbcParam.hppbak \
	BonCbc.cppbak BonCbc.hppbak \
	BonCbcLpStrategy.cppbak BonCbcLpStrategy.hppbak \
	BonCbcNlpStrategy.cppbak BonCbcNlpStrategy.hppbak \
	BonCbcNode.cppbak BonCbcNode.hppbak \
	BonChooseVariable.cppbak BonChooseVariable.hppbak \
	BonQPStrongBranching.cppbak BonQPStrongBranching.hppbak \
	BonCurvBranching.cppbak BonCurvBranching.hppbak \
        BonLpStrongBranching.hppbak BonLpStrongBranching.cppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_FILES)

DISTCLEANFILES = \
	CbcCompareUser.cpp CbcBranchUser.cpp

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
