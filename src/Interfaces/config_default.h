
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* this needs to come before the include of config_ipopt_default.h */
#ifndef BONMINLIB_EXPORT
#if defined(_WIN32) && defined(DLL_EXPORT)
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
