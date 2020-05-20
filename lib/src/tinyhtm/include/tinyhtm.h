/** \file
    \brief      Convenience header that includes all tinyhtm library headers. 

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech

    \mainpage   Overview

    Tiny HTM is a minimal C99 library providing spherical geometry
    primitives and methods for HTM indexing them. The implementation
    also exports selection algorithms in its API; these are required
    by the library implementation but may be more generally useful.

    Please refer to the modules tab to see a list of functions and types
    grouped by area of functionality. In general:

    <table style="width:75%; border:none; border-spacing:1em">
    <tr><td> tinyhtm/common.h </td>
        <td> contains utilities, including a defintion of the
             \link htm_errcode error codes \endlink used by the library
             and the means to map those to error messages.
        </td>
    </tr>
    <tr><td> tinyhtm/geometry.h </td>
        <td> contains types and functions related to spherical geometry,
             including \link htm_v3 vectors \endlink,
             \link htm_sc spherical coordinates \endlink, 
             \link htm_s2ellipse ellipses \endlink and
             \link htm_s2cpoly convex polygons \endlink. Spherical circles
             are easily represented as a position and a radius.
        </td>
    </tr>
    <tr><td> tinyhtm/htm.h </td>
        <td> contains types and functions for HTM indexing of points,
             circles, ellipses and polygons.
        </td>
    </tr>
    <tr><td> tinyhtm/select.h </td>
        <td> contains linear time selection algorithms over arrays of
             doubles, e.g. for finding medians.
        </td>
    </tr>
    <tr><td> tinyhtm/tree.h </td>
        <td> contains algorithms for determining how many points are
             inside a region. Points and HTM trees over them are external,
             allowing these algorithms to operate quickly on very large
             data sets (billions of points).
        </td>
    </tr>
    </table>

    The top-level header tinyhtm.h will include all of the above.

    To link against the library (which is static), pass
    \c -ltinyhtm \c -lm to the linker.
  */
#ifndef HTM_TINYHTM_H
#define HTM_TINYHTM_H

#include "tinyhtm/common.h"
#include "tinyhtm/select.h"
#include "tinyhtm/geometry.h"
#include "tinyhtm/htm.h"
#include "tinyhtm/tree.h"

#endif /* HTM_TINYHTM_H */

