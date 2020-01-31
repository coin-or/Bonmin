
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_ipopt_default.h */
#ifndef BONMINLIB_EXPORT
#ifdef _WIN32
/* assuming we build an Ipopt DLL */
#define BONMINLIB_EXPORT __declspec(dllexport)
#else
#define BONMINLIB_EXPORT
#endif
#endif

/* include the public project specific macros */
#include "config_bonmin_default.h"

/***************************************************************************/
/*        HERE DEFINE THE PROJECT SPECIFIC PRIVATE MACROS                  */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* If defined, the Ampl Solver Library is available. */
#define COIN_HAS_ASL 1
