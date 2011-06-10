# Copyright (C) 2006, 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

if COIN_HAS_ASL
  bin_PROGRAMS = bonminbcp
endif

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
	$(CPXLIB)
endif

IPOPTLIBS = \
	$(IPOPTOBJDIR)/src/Interfaces/libipopt.la

LIBS += `cat $(IPOPTOBJDIR)/src/Interfaces/ipopt_addlibs_cpp.txt` \
	`cat $(BCPOBJDIR)/bcp_addlibs.txt` \
	`cat $(CBCOBJDIR)/cbc_addlibs.txt` \
	`cat $(CGLOBJDIR)/cgl_addlibs.txt` \
	`cat $(OSIOBJDIR)/osi_addlibs.txt` \
	`cat $(CLPOBJDIR)/clp_addlibs.txt` \
	`cat $(COINUTILSOBJDIR)/coinutils_addlibs.txt`

########################################################################
#                                bonminbcp                                #
########################################################################

bonminbcp_SOURCES = BB_cut.hpp BM.hpp \
	BB_cut.cpp BM.cpp BM_tm.cpp BM_lp.cpp BM_lp_branch.cpp BM_pack.cpp

bonminbcp_DEPENDENCIES = \
	../../src/CbcBonmin/libbonminampl.la \
	../../src/CbcBonmin/libbonmin.la \
	$(IPOPTOBJDIR)/src/Apps/AmplSolver/libamplinterface.la \
	$(IPOPTLIBS) \
	$(LIBCBCCOPY) \
	$(COINLIBS) \
	amplsolver.a

AMPL_PATCHED_FILES = \
	asl.h \
	jac0dim.c \
	misc.c \
	pfg_read.c
AMPL_SRC = $(abs_source_dir)/../ThirdParty/ASL/solvers

amplsolver.a:
	rm -rf solvers
	mkdir -p solvers ; \
	cd solvers ; \
	for f in $(AMPL_SRC)/*.c* $(AMPL_SRC)/*.h* $(AMPL_SRC)/makefile.u; do \
		$(LN_S) $$f ; \
	done ; \
	for f in $(AMPL_PATCHED_FILES); do \
		rm $$f ; \
		cp $(AMPL_SRC)/$$f . ; \
	done ; \
	patch -p1 < $(abs_source_dir)/experimental/Bcp/ampl_bcp.patch ; \
	for f in `grep -l 'strtod' *.c *.h`; do \
		mv $$f $$f.orig ; \
		sed -e 's/strtod/ASLstrtod/g' -e 's/ASLASLstrtod/ASLstrtod/g' $$f.orig > $$f ; \
		rm $$f.orig ; \
	done ; \
	for f in `grep -l 'atof(' *.c`; do \
		mv $$f $$f.orig ; \
		sed -e 's/atof/ASLatof/g' -e 's/ASLASLatof/ASLatof/g' $$f.orig > $$f ; \
		rm $$f.orig ; \
	done ; \
	rm -f aslflags.txt ; \
	grep ASLMAKEFLAGS ../../../../ThirdParty/ASL/config.log | tail -1 > aslflags.txt ; \
	cat aslflags.txt ; \
	eval `sed -e "s/ CFLAGS=\".*\"//" aslflags.txt`; \
	eval `sed -e "s/.*CFLAGS=\"\(.*\)\".*/CFLAGS=\"-DNO_FUNCADD \1\"/" aslflags.txt` ; \
	if test "$(EXEEXT)" = ".exe"; then \
		sed -e 's/a.out/a.exe/' makefile.u > blabla ; \
		mv blabla makefile.u ; \
	fi ; \
	$(MAKE) CC="$(CC)" $$ASLMAKEFLAGS CFLAGS="$$CFLAGS" -f makefile.u ; \
	cp amplsolver.a ..
	rm -rf solvers

bonminbcp_LDADD = $(bonminbcp_DEPENDENCIES) \
	$(ADDLIBS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglProbing` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglMixedIntegerRounding` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglKnapsackCover` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src/CglGomory` \
	-I`$(CYGPATH_W) $(BCPSRCDIR)/src/include` \
	-I`$(CYGPATH_W) $(BCPOBJDIR)/inc`

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

CLEANFILES = $(LIBCBCCOPY) amplsolver.a

DISTCLEANFILES = 
