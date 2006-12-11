# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 105 2006-11-03 19:31:19Z pbonami $

# Author: Pierre Bonami, International Business Machines 12/04/2006

AUTOMAKE_OPTIONS = foreign

bin_PROGRAMS = fun

COINLIBS = \
	$(CBCOBJDIR)/src/libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la \
        $(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
        $(OSIOBJDIR)/src/libOsi.la \
        $(CLPOBJDIR)/src/libClp.la \
        $(COINUTILSOBJDIR)/src/libCoinUtils.la

if COIN_HAS_CPX
COINLIBS += $(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la \
        $(COINUTILSOBJDIR)/src/libCoinUtils.la \
        $(CPXLIB)
endif

if COIN_HAS_COUENNE
COINLIBS += $(COUENNEOBJDIR)/src/libCouenne.la
endif

IPOPTLIBS = \
        $(IPOPTOBJDIR)/src/Interfaces/libipopt.la

LIBS += `cat $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt`

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../../src/OaAlgorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Apps/AmplSolver` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/examples` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCbc`\
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
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` 

if COIN_HAS_CPX
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

AM_CPPFLAGS += $(ASL_CPPFLAGS)

if COIN_HAS_FILTERSQP
  AM_CPPFLAGS += -I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Filter` 
endif

if COIN_HAS_COUENNE
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/include`
endif

fun_SOURCES = main.cpp

fun_LDADD = \
	../../src/Interfaces/Ampl/libamplinterface.la \
        ../../src/Apps/libbonmin.la \
        $(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la \
	$(IPOPTLIBS) \
	$(ASLLIB) \
	$(LIBCBCCOPY) \
	$(COINLIBS) \
	$(ADDLIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
fun_LDFLAGS = -rpath $(libdir)


