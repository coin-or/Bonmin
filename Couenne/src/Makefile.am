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
	branch/CouenneThreeWayBranchObj.cpp \
	branch/CouenneBranchingObject.cpp \
	branch/CouenneObject.cpp \
	branch/CouenneChooseVariable.cpp \
	branch/getFixVarBinFun.cpp \
	branch/projections.c \
	branch/operators/branchExprAbs.cpp \
	branch/operators/branchExprExp.cpp \
	branch/operators/branchExprDiv.cpp \
	branch/operators/branchExprInv.cpp \
	branch/operators/branchExprLog.cpp \
	branch/operators/branchExprMul.cpp \
	branch/operators/branchExprPow.cpp \
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
	expression/operators/exprMinMax.cpp \
	expression/operators/exprGroup.cpp \
	expression/operators/exprQuad.cpp \
	expression/operators/impliedBounds-exprSum.cpp \
	expression/operators/impliedBounds-exprDiv.cpp \
	expression/operators/impliedBounds-exprMul.cpp \
	expression/exprOp.cpp \
	expression/exprUnary.cpp \
	expression/exprVar.cpp \
	expression/exprAux.cpp \
	expression/expression.cpp \
	simplify/simplifiers.cpp \
	simplify/simplify.cpp \
	convex/checkNLP.cpp \
	convex/addEnvelope.cpp \
	convex/createCuts.cpp \
	convex/operators/conv-exprAbs.cpp \
	convex/operators/conv-exprDiv.cpp \
	convex/operators/conv-exprMul.cpp \
	convex/operators/conv-exprMul-genCuts.cpp \
	convex/operators/conv-exprOpp.cpp \
	convex/operators/conv-exprPow.cpp \
	convex/operators/conv-exprPow-getBounds.cpp \
	convex/operators/conv-exprPow-envelope.cpp \
	convex/operators/powNewton.cpp \
	convex/operators/conv-exprSub.cpp \
	convex/operators/conv-exprSum.cpp \
	convex/operators/linearize.cpp \
	convex/operators/conv-exprInv.cpp \
	convex/operators/conv-exprSinCos.cpp \
	convex/operators/conv-exprExp.cpp \
	convex/operators/conv-exprLog.cpp \
	convex/operators/conv-exprGroup.cpp \
	convex/operators/trigNewton.cpp \
	convex/CouenneCutGenerator.cpp \
	convex/generateCuts.cpp \
	convex/boundTightening.cpp \
	convex/genColCuts.cpp \
	convex/genRowCuts.cpp \
	convex/obbt.cpp \
	convex/alphaConvexify.cpp \
	convex/quadCuts.cpp \
	problem/problem.cpp \
	problem/CouenneProblem.cpp \
	problem/problemIO.cpp \
	problem/extended.cpp \
	problem/constraint.cpp \
	problem/tightenBounds.cpp \
	problem/impliedBounds.cpp \
	problem/readASLfg.cpp \
	problem/fillQuadIndices.cpp \
	problem/constrStandardize.cpp \
	problem/elementBreak.cpp \
	problem/splitAux.cpp \
	problem/auxiliarize.cpp \
	problem/standardize.cpp \
	problem/depGraph/depGraph.cpp \
	util/drawCuts.cpp \
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
			 readnl/invmap.c
endif

########################################################################
#                            Additional flags                          #
########################################################################

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = -I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/expression/operators` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/branch` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/util` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/convex/sdpcuts` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem` \
	-I`$(CYGPATH_W) $(COUENNESRCDIR)/src/problem/depGraph` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(OSIOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CBCSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CBCOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CGLSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CGLOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(BONMINSRCDIR)/src/CbcBonmin` \
	-I`$(CYGPATH_W) $(objdir)/inc` \
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
	branch/CouenneThreeWayBranchObj.hpp \
	branch/CouenneBranchingObject.hpp \
	branch/CouenneChooseVariable.hpp \
	branch/projections.h \
	expression/operators/exprAbs.hpp \
	expression/operators/exprExp.hpp \
	expression/operators/exprInv.hpp \
	expression/operators/exprMul.hpp \
	expression/operators/exprOpp.hpp \
	expression/operators/exprBMul.hpp \
	expression/operators/exprLog.hpp \
	expression/operators/exprDiv.hpp \
	expression/operators/exprMin.hpp \
	expression/operators/exprBDiv.hpp \
	expression/operators/exprSum.hpp \
	expression/operators/exprMax.hpp \
	expression/operators/exprPow.hpp \
	expression/operators/exprSin.hpp \
	expression/operators/exprCos.hpp \
	expression/operators/exprBSin.hpp \
	expression/operators/exprBCos.hpp \
	expression/operators/exprSub.hpp \
	expression/operators/exprGroup.hpp \
	expression/operators/exprQuad.hpp \
	expression/exprConst.hpp \
	expression/exprIVar.hpp \
	expression/exprVar.hpp \
	expression/exprAux.hpp \
	expression/exprUnary.hpp \
	expression/expression.hpp \
	expression/exprClone.hpp \
	expression/exprBound.hpp \
	expression/exprCopy.hpp \
	expression/exprOp.hpp \
	expression/CouennePrecisions.h \
	problem/depGraph/depGraph.hpp \
	problem/CouenneProblemElem.hpp \
	problem/CouenneProblem.hpp \
	expression/CouenneTypes.h \
	util/rootQ.hpp \
	convex/CouenneCutGenerator.hpp

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
