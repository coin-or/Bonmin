/*
 * Include file for the configuration of Bonmin.
 *
 * On systems where the code is configured with the configure script
 * (i.e., compilation is always done with HAVE_CONFIG_H defined), this
 * header file includes the automatically generated header file, and
 * undefines macros that might configure with other Config.h files.
 *
 * On systems that are compiled in other ways (e.g., with the
 * Developer Studio), a header files is included to define those
 * macros that depend on the operating system and the compiler.  The
 * macros that define the configuration of the particular user setting
 * (e.g., presence of other COIN packages or third party code) are set
 * here.  The project maintainer needs to remember to update this file
 * and choose reasonable defines.  A user can modify the default
 * setting by editing this file here.
 *
 */

#ifndef __BONMINCONFIG_H__
#define __BONMINCONFIG_H__

#ifdef HAVE_CONFIG_H
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

#else /* HAVE_CONFIG_H */

/* include the COIN-wide system specific configure header */
#include "configall_system.h"

/***************************************************************************/
/*             HERE DEFINE THE CONFIGURATION SPECIFIC MACROS               */
/***************************************************************************/

/* If defined, debug sanity checks are performed during runtime */
/* #define COIN_DEBUG 1 */

/* If defined, the Ampl Solver Library is available. */
#define COIN_HAS_ASL 1

/* Define to 1 if the Ipopt package is used */
#define COIN_HAS_BONMIN 1

/* Define to 1 if MA27 is available */
#define HAVE_MA27 1

/* Define to 1 if MA57 is available */
/* #undef HAVE_MA57 */

/* Define to 1 if MC19 is available */
#define HAVE_MC19 1

/* Define to 1 if MUMPS is available */
/* #undef COIN_HAS_MUMPS */

/* Define to 1 if Pardiso is available */
/* #undef HAVE_PARDISO */

/* Define to 1 if you are using the parallel version of Pardiso */
/* #undef HAVE_PARDISO_PARALLEL */

/* Define to 1 if TAUCS is available */
/* #undef HAVE_TAUCS */

/* Define to 1 if WSMP is available */
/* #undef HAVE_WSMP */

#endif /* HAVE_CONFIG_H */

#endif /*__BONMINCONFIG_H__*/
