# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id: Makefile.am 434 2007-03-28 12:58:54Z andreasw $

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

# What is done in this directory:

bin_PROGRAMS = robotbonmin

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Apps` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/Branching` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/OaGenerators` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Algorithms/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ipopt` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/Interfaces/Ampl` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin` \
	-I`$(CYGPATH_W) $(srcdir)/../../src/CbcBonmin/Heuristics` \
	$(BONMINLIB_CFLAGS) $(ASL_CFLAGS)

########################################################################
#                                robotbonmin                                  #
########################################################################

robotbonmin_SOURCES = Robot-bonmin.cpp \
               RobotSetup.cpp RobotSetup.hpp \
               BonNWayChoose.cpp  BonNWayChoose.hpp  BonNWayObject.cpp  BonNWayObject.hpp
robotbonmin_LDADD = \
	../../src/Interfaces/Ampl/libbonamplinterface.la \
	../../src/CbcBonmin/libbonmin.la \
	$(ASL_LIBS) \
	$(BONMINLIB_LIBS)

########################################################################
#                         make distclean stuff                         #
########################################################################

CLEANFILES = 

DISTCLEANFILES = 
