# Copyright (C) 2006, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author:  Andreas Waechter           IBM    2006-04-13

include BuildTools/Makemain.inc

########################################################################
#                          Subdirectories                              #
########################################################################

# Here list all subdirectories into which to recurse
SUBDIRS = src/Interfaces \
	src/Algorithms \
	src/CbcBonmin \
	src/Apps

########################################################################
#                           Extra Targets                              #
########################################################################

.PHONY: doc test unitTest userman astyle nodeSolver doxydoc

test: unitTest

unitTest: all
	cd test; $(MAKE) test

nodeSolver: all
	cd src/Apps; $(MAKE) nodeSolver

doc:	userman

userman: $(top_srcdir)/doc/BONMIN_UsersManual.tex
	cd $(top_srcdir)/doc && pdflatex BONMIN_UsersManual.tex

# Doxygen documentation

doxydoc:
	doxygen doxydoc/doxygen.conf

clean-doxydoc:
	( cd doxydoc ; rm -rf html *.log *.tag )

clean-local: clean-doxydoc
	if test -r test/Makefile; then cd test; $(MAKE) clean; fi

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
#                  Installation of the .pc file                        #
########################################################################

pkgconfiglib_DATA = bonmin.pc
if COIN_HAS_ASL
  pkgconfiglib_DATA += bonminamplinterface.pc
endif
