# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 434 2007-03-28 12:58:54Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

# What is done in this directory:
bin_PROGRAMS = sepa

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/../../src/Apps \
	-I$(srcdir)/../../src/Algorithms \
	-I$(srcdir)/../../src/Algorithms/Branching \
	-I$(srcdir)/../../src/Algorithms/OaGenerators \
	-I$(srcdir)/../../src/Algorithms/Ampl \
	-I$(srcdir)/../../src/Interfaces \
	-I$(srcdir)/../../src/Interfaces/Ipopt \
	-I$(srcdir)/../../src/Interfaces/Ampl \
	-I$(srcdir)/../../src/CbcBonmin \
	-I$(srcdir)/../../src/CbcBonmin/Heuristics \
	$(BONMINLIB_CFLAGS) $(BONMINAMPLINTERFACELIB_CFLAGS)

# add libtool flags
AM_LDFLAGS = $(LT_LDFLAGS)

########################################################################
#                                sepa                                  #
########################################################################

sepa_SOURCES = Sepa.cpp \
               SepaSetup.cpp \
               SepaTMINLP2OsiLP.cpp \
               SepaHeuristicInnerApproximation.cpp \
               BonOuterDescription.cpp

sepa_LDADD = \
	../../src/Interfaces/Ampl/libbonamplinterface.la \
	../../src/CbcBonmin/libbonmin.la
