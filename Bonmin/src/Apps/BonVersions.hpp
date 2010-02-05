#ifndef BonVersions_H
#define BonVersions_H


#ifdef HAVE_CONFIG_H
#include "config_ipopt.h"
#include <string>

static std::string IPOPT_VERSION=PACKAGE_VERSION;
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

//static std::string BONMIN_VERSION=PACKAGE_VERSION;
#define BONMIN_VERSION "1.3"
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
#define BONMIN_VERSION="1.3"
#endif

#endif

