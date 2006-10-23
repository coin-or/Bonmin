# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

bin_PROGRAMS = bonmin nodeSolver 


lib_LTLIBRARIES = libbonmin.la

########################################################################
#                          bonmin library                              #
########################################################################

libbonmin_la_SOURCES =

nodist_libbonmin_la_SOURCES = CbcCompareUser.cpp CbcBranchUser.cpp

libbonmin_la_LIBADD = \
	../CbcBonmin/libbonminbb.la \
	../OaInterface/liboainterface.la \
	../IpoptInterface/libipoptinterface.la 

libbonmin_la_LDFLAGS = $(LT_LDFLAGS)

COINLIBS = \
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

if COIN_HAS_FSQP
  COINLIBS += \
	$(FSQPLIB)
endif

IPOPTLIBS = \
	$(IPOPTOBJDIR)/src/Interfaces/libipopt.la

LIBS += `cat $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt`

########################################################################
#                                bonmin                                #
########################################################################

bonmin_SOURCES = BonMin.cpp

bonmin_LDADD = \
	../BonminAmplInterface/libamplinterface.la\
	libbonmin.la \
	$(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la \
	$(IPOPTLIBS) \
	$(ASLLIB) \
	$(LIBCBCCOPY) \
	$(COINLIBS) \
        $(ADDLIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
bonmin_LDFLAGS = -rpath $(libdir)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../OaInterface` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../BonminAmplInterface` \
	-I`$(CYGPATH_W) $(srcdir)/../CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
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
#                              nodeSolver                              #
########################################################################

nodeSolver_SOURCES = nodeSolver.cpp

nodeSolver_LDADD = \
	../BonminAmplInterface/libamplinterface.la \
	libbonmin.la \
	$(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la \
	$(IPOPTLIBS) \
	$(ASLLIB) \
	$(LIBCBCCOPY) \
	$(COINLIBS)

# Finally, the -rpath flag is used by libtool to make sure that the shared
# library is found (in the lib install directory) when we are using dynamic
# libraries.
nodeSolver_LDFLAGS = -rpath $(libdir)


########################################################################
#                         make distclean stuff                         #
########################################################################

CLEANFILES = $(LIBCBCCOPY)

DISTCLEANFILES = \
	CbcCompareUser.cpp CbcBranchUser.cpp

