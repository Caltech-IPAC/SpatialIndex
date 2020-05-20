/** \file
    \brief      Common macros/functions.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#ifndef HTM_COMMON_H
#define HTM_COMMON_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

/* JCG
#include "config.h"
JCG */


/** \def HTM_UNUSED
    Marks a function argument as unused.
 */
#if HAVE_ATTRIBUTE_UNUSED
#   define HTM_UNUSED __attribute__ ((unused))
#else
#   define HTM_UNUSED
#endif

/** \def HTM_INLINE
    Marks a function as inline.
 */
#if __STDC_VERSION__ >= 199901L || __GNUC__
/* "inline" is a keyword */
#   define HTM_INLINE static inline
#else
#   define HTM_INLINE static HTM_UNUSED
#endif

/** \def HTM_ALIGNED(x)
    Marks a structure as requiring alignment to x bytes.
 */
#if HAVE_ATTRIBUTE_ALIGNED
#   define HTM_ALIGNED(x) __attribute__ ((aligned(x)))
#else
#   define HTM_ALIGNED(x)
#endif

/** \def HTM_ISNAN(x)
    Returns 1 if \p x is a floating point NaN, 0 otherwise.
    \def HTM_ISSPECIAL(x)
    Returns 1 if \p x is NaN or +/- Inf, 0 otherwise.
 */
#if __STDC_VERSION__ >= 199901L
#   define HTM_ISNAN(x) isnan(x)
#   define HTM_ISSPECIAL(x) (!isfinite(x))
#else
#   define HTM_ISNAN(x) ((x) != (x))
#   define HTM_ISSPECIAL(x) ((x) != (x) || ((x) != 0.0 && (x) == 2*(x)))
#endif

/* M_PI may not be available in <math.h> */
#ifndef M_PI
#   define M_PI 3.141592653589793238462643
#endif


#ifdef __cplusplus
extern "C" {
#endif

/** \defgroup utils Utilities
    @{
  */

/** Represents a range of integers.
  */
struct htm_range {
    int64_t min;  /**< Lower bound (inclusive) */
    int64_t max;  /**< Upper bound (inclusive) */
};

/** Counts the number of 1 bits in a 64 bit integer.
  */
HTM_INLINE int htm_popcount(uint64_t x)
{
    const uint64_t c1 = UINT64_C(0x5555555555555555); /* -1/3 */
    const uint64_t c2 = UINT64_C(0x3333333333333333); /* -1/5 */
    const uint64_t c4 = UINT64_C(0x0f0f0f0f0f0f0f0f); /* -1/17 */
    const uint64_t f  = UINT64_C(0x0101010101010101); /* -1/255 */
    x = x - ((x >> 1)  & c1);
    x = (x & c2) + ((x >> 2)  & c2);
    x = (x + (x >> 4)) & c4;
    return (int) ((x * f) >> 56);
}

/** Library error codes.
  */
enum htm_errcode {
    HTM_OK = 0,     /**< Success */
    HTM_ENOMEM,     /**< Memory (re)allocation failed */
    HTM_ENULLPTR,   /**< NULL pointer argument */
    HTM_ENANINF,    /**< NaN or +/-Inf argument */
    HTM_EZERONORM,  /**< Input vector has zero norm */
    HTM_ELAT,       /**< Latitude angle out-of-bounds */
    HTM_EANG,       /**< Invalid radius, width, or height (angle) */
    HTM_EHEMIS,     /**< Input vertices are non-hemispherical */
    HTM_ELEN,       /**< Too many/too few array elements */
    HTM_EDEGEN,     /**< Degenerate vertices */
    HTM_EID,        /**< Invalid HTM ID */
    HTM_ELEVEL,     /**< Invalid HTM subdivision level */
    HTM_EIO,        /**< IO operation failed */
    HTM_EMMAN,      /**< Failed to mmap(), madvise() or mlock() */
    HTM_EINV,       /**< Invalid argument */
    HTM_ETREE,      /**< Invalid HTM tree file, or tree/data file mismatch */
    HTM_NUM_CODES
};


/** Returns an error message for the given error code. The memory
    for this error message must never be freed.
  */
const char * htm_errmsg(enum htm_errcode err);


/** Degrees per radian. */
#define HTM_DEG_PER_RAD    57.2957795130823208767981548141

/** Radians per degree. */
#define HTM_RAD_PER_DEG    0.0174532925199432957692369076849

/** Arcseconds per degree. */
#define HTM_ARCSEC_PER_DEG 3600.0


/** Returns the given angle, range-reduced to lie in <tt>[0, 360)</tt> degrees.
  */
HTM_INLINE double htm_angred(double angle_deg)
{
    double angle = fmod(angle_deg, 360.0);
    if (angle < 0.0) {
        angle += 360.0;
        if (angle == 360.0) {
            angle = 0.0;
        }
    }
    return angle;
}

/** Returns the value of \p x clamped to <tt>[min, max]</tt>.
  */
HTM_INLINE double htm_clamp(double x, double min, double max)
{
    return (x < min) ? min : ((x > max) ? max : x);
}

/** @}
  */

#ifdef __cplusplus
}
#endif

#endif /* HTM_COMMON_H */

