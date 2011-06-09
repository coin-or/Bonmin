
/* include the COIN-OR-wide system specific configure header */
#include "configall_system.h"

/* include the public project specific macros */
#include "config_bonmin_default.h"

/***************************************************************************/
/*             HERE DEFINE THE PROJECT SPECIFIC MACROS                     */
/*    These are only in effect in a setting that doesn't use configure     */
/***************************************************************************/

/* Define to the debug sanity check level (0 is no test) */
#define COIN_BONMIN_CHECKLEVEL 0

/* Define to the debug verbosity level (0 is no output) */
#define COIN_BONMIN_VERBOSITY 0

/* If defined, the Ampl Solver Library is available. */
#define COIN_HAS_ASL 1
