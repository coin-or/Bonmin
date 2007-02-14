# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 405 2006-07-10 18:52:52Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                               libCbc                                 #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the 'lib' directory
lib_LTLIBRARIES = libCouenne.la

# List all source files for this library, including headers
libCouenne_la_SOURCES = \
	branch/CouenneBranchingObject.cpp \
	branch/CouenneObject.cpp \
	branch/CouenneChooseVariable.cpp \
	expression/operators/exprAbs.cpp \
	expression/operators/exprDiv.cpp \
	expression/operators/exprExp.cpp \
	expression/operators/exprInv.cpp \
	expression/operators/exprLog.cpp \
	expression/operators/exprMul.cpp \
	expression/operators/exprOpp.cpp \
	expression/operators/exprPow.cpp \
	expression/operators/exprSin.cpp \
	expression/operators/exprCos.cpp \
	expression/operators/exprSub.cpp \
	expression/operators/exprSum.cpp \
	expression/operators/exprMax.cpp \
	expression/operators/exprMin.cpp \
	expression/exprVar.cpp \
	expression/expression.cpp \
	simplify/simplifiers.cpp \
	simplify/simplify.cpp \
	convex/addEnvelope.cpp \
	convex/createCuts.cpp \
	convex/conv-exprOp.cpp \
	convex/conv-exprUnary.cpp \
	convex/operators/conv-exprAbs.cpp \
	convex/operators/conv-exprDiv.cpp \
	convex/operators/conv-exprMul.cpp \
	convex/operators/conv-exprOpp.cpp \
	convex/operators/conv-exprPow.cpp \
	convex/operators/conv-exprPow-getBounds.cpp \
	convex/operators/conv-exprPow-envelope.cpp \
	convex/operators/conv-exprSub.cpp \
	convex/operators/conv-exprSum.cpp \
	convex/operators/conv-exprInv.cpp \
	convex/operators/conv-exprSinCos.cpp \
	convex/operators/conv-exprExp.cpp \
	convex/operators/conv-exprLog.cpp \
	convex/CouenneCutGenerator.cpp \
	convex/generateCuts.cpp \
	problem/problem.cpp \
	problem/problemIO.cpp \
	problem/constraint.cpp \
	convex/update_conv.cpp \
	util/rootQ.c

COINLIBS = \
	$(OSIOBJDIR)/src/libOsi.la \
	$(OSIOBJDIR)/src/OsiClp/libOsiClp.la \
	$(CGLOBJDIR)/src/libCgl.la \
	$(CLPOBJDIR)/src/libClp.la \
	$(COINUTILSOBJDIR)/src/libCoinUtils.la

# This is for libtool (on Windows)
libCouenne_la_LDFLAGS = $(LT_LDFLAGS)

# If we have the Ampl solver library, we need to add additional things
if COIN_HAS_ASL
libCouenne_la_SOURCES += readnl/readnl.cpp \
			 readnl/nl2e.cpp \
			 readnl/invmap.c \
			 readnl/fix_asl_group.c
endif

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = -I`$(CYGPATH_W) $(COUENNESRCDIR)/src/include` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/include/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc`

# If we have the Ampl solver library, we need to add additional things
if COIN_HAS_ASL
AM_CPPFLAGS += $(ASL_CPPFLAGS)
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'install'
include_HEADERS = \
	../inc/config_couenne.h \
	branch/CouenneObject.hpp \
	branch/CouenneBranchingObject.hpp \
	branch/CouenneChooseVariable.hpp \
	include/operators/exprAbs.h \
	include/operators/exprExp.h \
	include/operators/exprInv.h \
	include/operators/exprMul.h \
	include/operators/exprOpp.h \
	include/operators/exprBMul.h \
	include/operators/exprLog.h \
	include/operators/exprDiv.h \
	include/operators/exprMin.h \
	include/operators/exprBDiv.h \
	include/operators/exprCos.h \
	include/operators/exprSum.h \
	include/operators/exprMax.h \
	include/operators/exprPow.h \
	include/operators/exprSin.h \
	include/operators/exprSub.h \
	include/exprConst.h \
	include/exprIVar.h \
	include/exprVar.h \
	include/exprAux.h \
	include/exprUnary.h \
	include/expression.h \
	include/exprClone.h \
	include/exprBound.h \
	include/exprCopy.h \
	include/exprOp.h \
	include/CouennePrecisions.h \
	include/CouenneProblemElem.h \
	include/CouenneProblem.h \
	include/CouenneTypes.h \
	include/rootQ.h \
	include/CouenneCutGenerator.h

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = CouenneConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)
