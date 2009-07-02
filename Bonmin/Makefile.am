# Copyright (C) 2006, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          Subdirectories                              #
########################################################################

# Here list all subdirectories into which to recurse
SUBDIRS = src/Interfaces \
	  src/Algorithms \
	  src/CbcBonmin \
	  src/Apps 

if BUILD_FP 
SUBDIRS += experimental/FP
endif

if COMPILE_BONMINBCP
SUBDIRS += experimental/Bcp
endif

# We don't want to compile the test subdirectory, unless the test target is
# specified.  But we need to list it as subdirectory to make sure that it is
# included in the tarball

if ALWAYS_FALSE
  SUBDIRS += test
endif

########################################################################
#             Additional files to be included in tarball               #
########################################################################

# Here we need include all files that are not mentioned in other Makefiles

EXTRA_DIST = \
	doc/Bonmin_ReferenceManual.hpp \
	doc/Bonm-Min_UsersManual.pdf \
	doc/Bonm-Min_UsersManual.tex \
	doc/BOUM_Intro.tex \
	doc/options.tex \
	test/InterfaceTest.cpp \
	test/CppExample/MyTMINLP.hpp \
	test/CppExample/MyTMINLP.cpp \
	test/CppExample/MyBonmin.cpp \
	test/mytoy.nl

########################################################################
#                           Extra Targets                              #
########################################################################

.PHONY: doc test unitTest userman astyle nodeSolver

test: unitTest

unitTest: all
	cd test; $(MAKE) test

nodeSolver: all
	cd src/Apps; $(MAKE) nodeSolver

doc:	userman

userman: $(top_srcdir)/doc/BONMIN_UsersManual.tex
	cd $(top_srcdir)/doc && pdflatex BONMIN_UsersManual.tex

install-exec-local: install-doc

uninstall-local: uninstall-doc

astyle:
	cd src/Algorithms; $(MAKE) astyle
	cd src/Algorithms/Ampl; $(MAKE) astyle
	cd src/Algorithms/Branching; $(MAKE) astyle
	cd src/Algorithms/OaGenerators; $(MAKE) astyle
	cd src/Apps; $(MAKE) astyle
	cd src/CbcBonmin; $(MAKE) astyle
	cd src/Interfaces; $(MAKE) astyle
	cd src/Interfaces/Ampl; $(MAKE) astyle
	cd src/Interfaces/Filter; $(MAKE) astyle
	cd src/Interfaces/Ipopt; $(MAKE) astyle

########################################################################
#                         Maintainer Stuff                             #
########################################################################

# Files that are generated and should be cleaned with make distclean
DISTCLEANFILES = \
	Test/*.$(OBJEXT) \
	Test/unitTest$(EXEEXT) \
	Test/CppExample/*.$(OBJEXT) \
	Test/CppExample/CppExample$(EXEEXT) \
	Test/testtoy.nl

include ../BuildTools/Makemain.inc
