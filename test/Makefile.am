# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.


# Author:  Andreas Waechter           IBM    2006-04-13

########################################################################
#                       unitTest for Bonmin                            #
########################################################################

noinst_PROGRAMS =

if COIN_HAS_ASL
  noinst_PROGRAMS += unitTest
endif

noinst_PROGRAMS += CppExample
#CExample

unitTest_SOURCES = InterfaceTest.cpp 

# list libraries that need to be linked in
unitTest_LDADD = ../src/CbcBonmin/libbonmin.la
if COIN_HAS_ASL
  unitTest_LDADD += ../src/CbcBonmin/libbonminampl.la
endif

#########################################################################
##                      Example C++ program                             #
#########################################################################

CppExample_SOURCES = MyBonmin.cpp  MyTMINLP.cpp  MyTMINLP.hpp

# list libraries that need to be linked in
CppExample_LDADD = ../src/CbcBonmin/libbonmin.la

#########################################################################
##                      Example C program                             #
#########################################################################

#CExample_SOURCES = MyBonminC.c  MyBonmin.h main.c

# list libraries that need to be linked in
#CExample_LDADD        = ../src/CbcBonmin/libbonmin.la

# Here list all include flags, relative to this "srcdir" directory.
AM_CPPFLAGS = \
	-I$(srcdir)/../src/Apps \
	-I$(srcdir)/../src/Interfaces \
	-I$(srcdir)/../src/Interfaces/Ipopt \
	-I$(srcdir)/../src/Interfaces/Filter \
	-I$(srcdir)/../src/Interfaces/Ampl \
	-I$(srcdir)/../src/CbcBonmin \
	-I$(srcdir)/../src/Algorithms \
	-I$(srcdir)/../src/Algorithms/Branching \
	-I$(srcdir)/../src/Algorithms/OaGenerators \
	-I$(srcdir)/../src/Algorithms/Ampl \
	$(BONMINLIB_CFLAGS) -UBONMINLIB_BUILD

#test: unitTest$(EXEEXT) CppExample$(EXEEXT) CExample$(EXEEXT)
test: unitTest$(EXEEXT) CppExample$(EXEEXT)
	./unitTest$(EXEEXT)
	./CppExample$(EXEEXT)
#	./CExample$(EXEEXT)

.PHONY: test
