#ifndef BonVersions_H
#define BonVersions_H


#ifdef HAVE_CONFIG_H
#include "config_ipopt.h"
#include <string>

/* undefine macros that could conflict with those in other config.h
   files */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include "config_bonmin.h"

/* undefine macros that could conflict with those in other config.h
   files */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION
#else

#include "IpoptConfig.h"
#define BONMIN_VERSION="1.4trunk"
#endif

#endif

