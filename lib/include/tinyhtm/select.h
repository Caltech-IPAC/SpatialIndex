/** \file
    \brief      Selection algorithms.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#ifndef HTM_SELECT_H
#define HTM_SELECT_H

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \defgroup select Selection algorithms
    @{
  */

/** Returns the k-th smallest value in an array of doubles (where k = 0 is
    the smallest element). The implementation guarantees O(n) runtime
    even when faced with an array of identical elements. After this function
    returns, the k-th largest element is stored in array[k], and the
    following invariants hold:

    - <tt>array[i] <= array[k]</tt> for <tt>i < k</tt>
    - <tt>array[i] >= array[k]</tt> for <tt>i > k</tt>

    Assumes that:
    - <tt>array != NULL</tt>
    - <tt>n > 0</tt>
    - <tt>k < n</tt>

    \param[in] array    Array to select from.
    \param[in] n        Number of elements in array.
    \param[in] k        The element to select.
  */
double htm_select(double *array, size_t n, size_t k);

/** Finds the k-th smallest value in an array of doubles (where k = 0 is the
    smallest element) using the linear time median-of-medians algorithm.
    Except for pathological inputs (e.g. arrays of mostly/entirely identical
    values or median-of-3 killer sequences), htm_select() will be
    faster and should be preferred. Note that htm_select() detects
    such pathological sequences and switches to the median-of-medians
    algorithm if necessary, so its runtime is also linear.

    Assumes that:
    - <tt>array != NULL</tt>
    - <tt>n > 0</tt>
    - <tt>k < n</tt>

    This function has the same inputs ond enforces the same invariants as
    htm_select().
  */
double htm_selectmm(double *array, size_t n, size_t k);

/** Returns the smallest value in an array of doubles.

    Assumes that:
    - <tt>array != NULL</tt>
    - <tt>n > 0</tt>

    \param[in] array    Array to select from.
    \param[in] n        Number of elements in array.
  */
double htm_min(const double *array, size_t n);

/** @}
  */

#ifdef __cplusplus
}
#endif

#endif /* HTM_SELECT_H */

