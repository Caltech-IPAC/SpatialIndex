/** \file
    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#include "tinyhtm/common.h"

#ifdef __cplusplus
extern "C" {
#endif


static const char * const _htm_inv_errcode =
    "Invalid error code passed to htm_errmsg()";

static const char * const _htm_errmsg[HTM_NUM_CODES] = {
    /* HTM_OK */        "Success",
    /* HTM_ENOMEM */    "Memory allocation or reallocation failed",
    /* HTM_ENULLPTR */  "NULL pointer",
    /* HTM_ENANINF */   "Value is NaN or +/-Inf",
    /* HTM_EZERONORM */ "Vector has zero norm",
    /* HTM_ELAT */      "Latitude angle not in range [-90, 90] degrees",
    /* HTM_EANG */      "Radius, width, height or angle is negative or too large",
    /* HTM_EHEMIS */    "Vectors (vertices/points) are not hemispherical",
    /* HTM_ELEN */      "Too many/too few array elements (vertices/points)",
    /* HTM_EDEGEN */    "Vectors (vertices/points) are degenerate",
    /* HTM_EID */       "Invalid HTM ID",
    /* HTM_ELEVEL */    "Invalid HTM subdivision level",
    /* HTM_EIO */       "IO operation failed",
    /* HTM_EMMAN */     "Failed to mmap(), madvise(), or mlock()",
    /* HTM_EINV */      "Invalid argument",
    /* HTM_ETREE */     "Invalid HTM tree file, or tree/data file mismatch"
};


const char * htm_errmsg(enum htm_errcode err)
{
    if ((int) err < 0 || err >= HTM_NUM_CODES) {
        return _htm_inv_errcode;
    }
    return _htm_errmsg[err];
}


#ifdef __cplusplus
}
#endif

