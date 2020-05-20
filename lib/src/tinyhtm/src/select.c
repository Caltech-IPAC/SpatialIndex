/** \file
    \brief      Selection algorithm implementations.

    For API documentation, see select.h.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#include "tinyhtm/select.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  Lookup tables for the 4 and 5 element median finding algorithms.
    Computed with the following python 2.6+ script (with n = 4, 5):

    @code
    import itertools

    def computeLut(n):
        nbits = (n * (n - 1)) / 2
        lut = [-1] * 2**nbits
        array = range(n)
        median = array[len(array) >> 1]
        for p in itertools.permutations(array):
            res = []
            for i in xrange(n - 1):
                for j in xrange(i + 1, n):
                    res.append(1 if p[i] < p[j] else 0)
            index = 0
            for i in xrange(len(res)):
                index += res[i] << (len(res) - i - 1)
            lut[index] = p.index(median)
        return lut
    @endcode
 */
static const signed char _htm_lut4[64] = {
     1, 1,-1, 3, 2,-1, 2, 3,-1,-1,-1, 0,-1,-1,-1, 0,
    -1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1, 3, 2,
     0, 0,-1,-1,-1,-1,-1,-1,-1, 3,-1, 1,-1,-1,-1,-1,
     2,-1,-1,-1, 1,-1,-1,-1, 2, 3,-1, 1, 1,-1, 3, 2
};

static const signed char _htm_lut5[1024] = {
     2, 2,-1, 4, 3,-1, 3, 4,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1,-1,-1,-1, 4, 3,
     1, 1,-1,-1,-1,-1,-1,-1,-1, 4,-1, 2,-1,-1,-1,-1, 3,-1,-1,-1, 2,-1,-1,-1, 3, 4,-1, 2, 2,-1, 4, 3,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 3,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1, 3,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1,-1,-1,-1, 4,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1, 2,-1, 4,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     1, 1,-1,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1,-1, 3, 4,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 4, 3,-1, 3, 4,-1, 2, 2,
     2, 2,-1, 4, 3,-1, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0,-1, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1, 4, 3,-1,-1,-1,-1,-1,-1,-1, 3,-1,-1,-1,-1,-1,-1, 4,-1,-1,-1,-1,-1,-1,-1, 1, 1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 4,-1, 2,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1, 4,-1,-1,-1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     3,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     3, 4,-1, 2, 2,-1, 4, 3,-1,-1,-1, 2,-1,-1,-1, 3,-1,-1,-1,-1, 2,-1, 4,-1,-1,-1,-1,-1,-1,-1, 1, 1,
     3, 4,-1,-1,-1,-1,-1,-1,-1, 1,-1, 1,-1,-1,-1,-1, 1,-1,-1,-1, 1,-1,-1,-1, 4, 3,-1, 3, 4,-1, 2, 2
};

/*  Returns the index of the median of 2 doubles.
 */
HTM_INLINE size_t _htm_median2(const double *array)
{
    return (array[0] < array[1]) ? 0 : 1;
}

/*  Returns the index of the median of 3 doubles.
 */
HTM_INLINE size_t _htm_median3(const double *array)
{
    double x0 = array[0];
    double x1 = array[1];
    double x2 = array[2];
    if (x0 < x1) {
        return x1 < x2 ? 1 : (x0 < x2 ? 2 : 0);
    } else {
        return x1 < x2 ? (x0 < x2 ? 0 : 2) : 1;
    }
}

/*  Returns the index of the median of 4 doubles using a branchless algorithm.
 */
HTM_INLINE size_t _htm_median4(const double *array)
{
    double a = array[0];
    double b = array[1];
    double c = array[2];
    double d = array[3];
    /* avoids branches, but always requires 6 comparisons. */
    int i = (((int) (a < b)) << 5) |
            (((int) (a < c)) << 4) |
            (((int) (a < d)) << 3) |
            (((int) (b < c)) << 2) |
            (((int) (b < d)) << 1) |
             ((int) (c < d));
    return (size_t) _htm_lut4[i];
}

/*  Returns the index of the median of 5 doubles using a branchless algorithm.
 */
HTM_INLINE size_t _htm_median5(const double *array)
{
    double a = array[0];
    double b = array[1];
    double c = array[2];
    double d = array[3];
    double e = array[4];
    /* avoids branches, but always performs 10 comparisons */
    int i = (((int) (a < b)) << 9) |
            (((int) (a < c)) << 8) |
            (((int) (a < d)) << 7) |
            (((int) (a < e)) << 6) |
            (((int) (b < c)) << 5) |
            (((int) (b < d)) << 4) |
            (((int) (b < e)) << 3) |
            (((int) (c < d)) << 2) |
            (((int) (c < e)) << 1) |
             ((int) (d < e));
    return (size_t) _htm_lut5[i];
}

/*  Returns the index of the median of medians for an array.

    The following pre-conditions are assumed:
        -   array != 0
        -   n > 0
 */
static size_t _htm_mm(double *array, size_t n)
{
    size_t i, j, m = 0;
    while (1) {
        if (n <= 5) {
            switch (n) {
                case 1: m = 0; break;
                case 2: m = _htm_median2(array); break;
                case 3: m = _htm_median3(array); break;
                case 4: m = _htm_median4(array); break;
                case 5: m = _htm_median5(array); break;
            }
            break;
        }
        for (i = 0, j = 0; i < n - 4; i += 5, ++j) {
            size_t m5 = _htm_median5(array + i) + i;
            double tmp = array[j];
            array[j] = array[m5];
            array[m5] = tmp;
        }
        n = j;
    }
    return m;
}

/*  Partitions the given array around the value of the i-th element, and
    returns the index of the pivot value after partitioning.

    If the array contains many values identical to the pivot value,
    lop-sided partitions can be generated by a naive algorithm. In the
    worst case (e.g. all elements are identical), an array of n elements will
    be partitioned into two sub-arrays of size 1 and n - 1. This leads to
    quadratic run-time for selection/sorting algorithms that use the
    partitioning primitive.

    This implementation therefore counts values identical to the pivot
    while partitioning. If the partitions are overly lop-sided, a second
    pass over the larger partition is performed. This pass assigns pivot
    values to the smaller partition until the sizes are balanced or there
    are no more pivot values left.

    The run-time of this function is O(n).

    The following pre-conditions are assumed:
        -   array != 0
        -   n > 0
        -   i < n
 */
static size_t _htm_wcpart(double *array, size_t n, size_t i)
{
    size_t u, v, neq;
    double pivot = array[i];
    array[i] = array[n - 1];
    /* partition around pivot */
    for (u = 0, v = 0, neq = 0; v < n - 1; ++v) {
        if (array[v] < pivot) {
            double tmp = array[u];
            array[u] = array[v];
            array[v] = tmp;
            ++u;
        } else if (array[v] == pivot) {
            ++neq;
        }
    }
    array[n - 1] = array[u];
    array[u] = pivot;
    if (neq > 0 && u < (n >> 2)) {
        /* lop-sided partition - use values identical
           to the pivot value to increase u */
        if (u + neq > (n >> 1)) {
            neq = (n >> 1) - u;
        }
        for (v = u + 1; neq != 0; ++v) {
            if (array[v] == pivot) {
                ++u;
                array[v] = array[u];
                array[u] = pivot;
                --neq;
            }
        }
    }
    return u;
}

/*  Chooses a pivot value using the median-of-3 strategy.

    The following pre-conditions are assumed:
        -   array != 0
        -   n > 0
 */
static size_t _htm_med3pivot(double *array, size_t n) {
    double a, b, c;
    size_t m;
    if (n <= 5) {
        switch (n) {
            case 1: return 0;
            case 2: return _htm_median2(array);
            case 3: return _htm_median3(array);
            case 4: return _htm_median4(array);
            case 5: return _htm_median5(array);
        }
    }
    m = n >> 1;
    a = array[0];
    b = array[m];
    c = array[n - 1];
    if (a < b) {
        return b < c ? m : (a < c ? n - 1 : 0);
    } else {
        return b < c ? (a < c ? 0 : n - 1) : m;
    }
}

/*  Partitions the given array around the value of the i-th element, and
    returns the index of the pivot value after partitioning.

    The following pre-conditions are assumed:
        -   array != 0
        -   n > 0
        -   i < n
 */
static size_t _htm_part(double *array, size_t n, size_t i)
{
    size_t u, v;
    const double pivot = array[i];
    array[i] = array[n - 1];
    for (u = 0, v = 0; v < n - 1; ++v) {
        if (array[v] < pivot) {
            double tmp = array[u];
            array[u] = array[v];
            array[v] = tmp;
            ++u;
        }
    }
    array[n - 1] = array[u];
    array[u] = pivot;
    return u;
}


double htm_select(double *array, size_t n, size_t k)
{
    size_t tot = 0;
    const size_t thresh = n*3;

    while (1) {
        size_t i = _htm_med3pivot(array, n);
        i = _htm_part(array, n, i);
        if (k == i) {
            break;
        } else if (k < i) {
            n = i;
        } else {
            array += i + 1;
            n -= i + 1;
            k -= i + 1;
        }
        tot += n;
        if (tot > thresh) {
            return htm_selectmm(array, n, k);
        }
    }
    return array[k];
}


double htm_selectmm(double *array, size_t n, size_t k)
{
    while (1) {
        size_t i = _htm_mm(array, n);
        i = _htm_wcpart(array, n, i);
        if (k == i) {
            break;
        } else if (k < i) {
            n = i;
        } else {
            array += i + 1;
            n -= i + 1;
            k -= i + 1;
        }
    }
    return array[k];
}


double htm_min(const double *array, size_t n)
{
    double m;
    size_t i;

    m = array[0];
    for (i = 1; i < n; ++i) {
        if (array[i] < m) {
            m = array[i];
        }
    }
    return m;
}


#ifdef __cplusplus
}
#endif

