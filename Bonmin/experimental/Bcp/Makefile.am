# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

bin_PROGRAMS = bonminbcp

########################################################################

COINLIBS = \
	$(BCPOBJDIR)/src/libBcp.la \
	$(CBCOBJDIR)/src/libCbc.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(OSIOBJDIR)/src/libOsi.la \
	$(CLPOBJDIR)/src/libClp.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

if COIN_HAS_CPX
  COINLIBS += \
	$(OSIOBJDIR)/src/OsiCpx/libOsiCpx.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la \
	$(CPXLIB)
endif

IPOPTLIBS = \
	$(IPOPTOBJDIR)/src/Interfaces/libipopt.la

LIBS += `cat $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt`

########################################################################
#                                bonmin                                #
########################################################################

bonminbcp_SOURCES = BB_cut.hpp BB_cut.cpp BM.hpp BM.cpp BM_pack.cpp bm_var.hpp

bonminbcp_LDADD = \
	../../src/Apps/libbonmin.la \
	../../src/BonminAmplInterface/libamplinterface.la \
	$(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la \
	$(IPOPTLIBS) \
	$(ASLLIB) \
	$(LIBCBCCOPY) \
	$(COINLIBS) \
	$(ADDLIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
bonminbcp_LDFLAGS = -rpath $(libdir)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../../src/OaInterface` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/IpoptInterface` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/BonminAmplInterface` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(BCPSRCDIR)/src/include`

if COIN_HAS_CPX
  AM_CPPFLAGS += \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiCpx` \
	-I`$(CYGPATH_W) $(CPXINCDIR)`
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                         make distclean stuff                         #
########################################################################

CLEANFILES = $(LIBCBCCOPY)

DISTCLEANFILES = 
