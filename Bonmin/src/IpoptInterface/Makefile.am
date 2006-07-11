# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Author:  Andreas Waechter           IBM    2006-04-13

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          libipoptinterface                           #
########################################################################

# Name of the library in this directory
noinst_LTLIBRARIES = libipoptinterface.la 
# List all source files, including headers
libipoptinterface_la_SOURCES = \
	IpCbcBoundsReader.cpp IpCbcBoundsReader.hpp \
	IpCbcColReader.cpp IpCbcColReader.hpp \
	IpCbcStartPointReader.cpp IpCbcStartPointReader.hpp \
	IpoptIntegerBias.hpp \
	IpoptInterface.cpp IpoptInterface.hpp \
	IpoptInteriorWarmStarter.cpp IpoptInteriorWarmStarter.hpp \
	IpoptWarmStart.cpp IpoptWarmStart.hpp \
	TMINLP2TNLP.cpp TMINLP2TNLP.hpp \
	TMINLP.hpp \
	TNLP2FPNLP.cpp TNLP2FPNLP.hpp

# This is for libtool (on Windows)
libipoptinterface_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows
AM_CPPFLAGS = \
	-I`$(CYGPATH_W) $(IPOPTOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Common` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/LinAlg` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Algorithm` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Interfaces` \
	-I`$(CYGPATH_W) $(IPOPTSRCDIR)/src/Apps/AmplSolver` \
	-I`$(CYGPATH_W) $(COINUTILSSRCDIR)/src` \
	-I`$(CYGPATH_W) $(COINUTILSOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src` \
	-I`$(CYGPATH_W) $(OSISRCDIR)/src/OsiClp` \
	-I`$(CYGPATH_W) $(CLPOBJDIR)/inc` \
	-I`$(CYGPATH_W) $(CLPSRCDIR)/src` 

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

########################################################################
#                Headers that need to be installed                     #
########################################################################

include_HEADERS = \
	IpoptInterface.hpp \
	IpoptInteriorWarmStarter.hpp \
	TMINLP2TNLP.hpp \
	TMINLP.hpp \
	TNLP2FPNLP.hpp

#############################################################################
# Create the Config.h file that always defines HAVE_CONFIG_H and install it #
#############################################################################

# You only need to adapt the following line
ConfigHeader = BonminConfig.h

install-exec-local:
	echo "#ifndef HAVE_CONFIG_H" >bla
	echo "#define HAVE_CONFIG_H" >>bla
	echo "#endif" >> bla
	cat $(srcdir)/$(ConfigHeader) >> bla
	$(install_sh_DATA) bla $(DESTDIR)$(includedir)/$(ConfigHeader)
	rm -f bla

uninstall-local:
	rm -f $(DESTDIR)$(includedir)/$(ConfigHeader)

########################################################################
#                            Astyle stuff                              #
########################################################################

# Here repeat all source files, with "bak" appended
ASTYLE_FILES = \
	IpCbcBoundsReader.cppbak IpCbcBoundsReader.hppbak \
	IpCbcColReader.cppbak IpCbcColReader.hppbak \
	IpCbcStartPointReader.cppbak IpCbcStartPointReader.hppbak \
	IpoptIntegerBias.hppbak \
	IpoptInterface.cppbak IpoptInterface.hppbak \
	IpoptInteriorWarmStarter.cppbak IpoptInteriorWarmStarter.hppbak \
	IpoptWarmStart.cppbak IpoptWarmStart.hppbak \
	TMINLP2TNLP.cppbak TMINLP2TNLP.hppbak \
	TMINLP.hppbak \
	TNLP2FPNLP.cppbak TNLP2FPNLP.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@
