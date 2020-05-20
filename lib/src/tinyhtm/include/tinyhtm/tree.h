/** \file
    \brief      HTM tree indexes and queries.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#ifndef HTM_TREE_H
#define HTM_TREE_H

#include "geometry.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ================================================================ */
/** \defgroup tree HTM tree indexes
    @{
  */
/* ================================================================ */

/** An entry (a row ID and a position) in an HTM tree.
  */
struct htm_tree_entry {
    struct htm_v3 v;   /**< Unit vector position. */
    int64_t rowid;     /**< Row ID. */
} HTM_ALIGNED(16);


/** An HTM tree containing a list of points sorted on HTM ID (tree
    entries), and optionally an index over the points that allows for
    fast spatial searches/counts.
  */
struct htm_tree {
    uint64_t leafthresh;    /**< Min # of points in an internal node. */
    uint64_t count;         /**< Total # of points in tree. */
    const unsigned char *root[8];           /**< Pointers to HTM root nodes. */
    const struct htm_tree_entry *entries;   /**< Data file memory map. */
    const void *index;      /**< Tree file memory map. */
    size_t indexsz;         /**< Size of tree file memory-map (bytes). */
    size_t datasz;          /**< Size of data file memory-map (bytes). */
    int indexfd;            /**< File descriptor for tree file. */
    int datafd;             /**< File descriptor for data file. */
} HTM_ALIGNED(16);

/** Initializes an HTM tree from the given tree and data files,
    returning HTM_OK on success. The \p treefile argument may be NULL,
    in which case queries result in scans of the points in \p datafile.
  */
enum htm_errcode htm_tree_init(struct htm_tree *tree,
                               const char * const treefile,
                               const char * const datafile);

/** Releases the resources for the given HTM tree.
  */
void htm_tree_destroy(struct htm_tree *tree);

/** Locks the HTM tree index in memory. If the associated data file
    is of size datathresh or less, it is locked in memory as well.
  */
enum htm_errcode htm_tree_lock(struct htm_tree *tree, size_t datathresh);


/* ================================================================ */
/** @}
    \defgroup tree_query HTM tree index queries
    @{
  */
/* ================================================================ */

/** Returns the number of points in \p tree that are inside the
    spherical circle with the given center and radius.

    The implementation scans over all the points in the tree rather than
    taking advantage of the index. Therefore, this function is primarily
    useful for testing - application code should use htm_tree_s2circle_count()
    instead.

    If an error occurs, the return value is negative, and \p *err
    is set to an error code describing the reason for the failure.
  */
int64_t htm_tree_s2circle_scan(const struct htm_tree *tree,
                               const struct htm_v3 *center,
                               double radius,
                               enum htm_errcode *err);

/** Returns the number of points in \p tree that are inside
    the given spherical ellipse.

    The implementation scans over all the points in the tree rather than
    taking advantage of the index. Therefore, this function is primarily
    useful for testing - application code should use htm_tree_s2ellipse_count()
    instead.

    If an error occurs, the return value is negative, and \p *err
    is set to an error code describing the reason for the failure.
  */
int64_t htm_tree_s2ellipse_scan(const struct htm_tree *tree,
                                const struct htm_s2ellipse *ellipse,
                                enum htm_errcode *err);

/** Returns the number of points in \p tree that are inside
    the given spherical convex polygon.

    The implementation scans over all the points in the tree rather than
    taking advantage of the index. Therefore, this function is primarily
    useful for testing - application code should use htm_tree_s2cpoly_count()
    instead.

    If an error occurs, the return value is negative, and \p *err
    is set to an error code describing the reason for the failure.
  */
int64_t htm_tree_s2cpoly_scan(const struct htm_tree *tree,
                              const struct htm_s2cpoly *poly,
                              enum htm_errcode *err);


/** Returns the number of points in \p tree that are inside the
    spherical circle with the given center and radius.

    If an error occurs, the return value is negative, and \p *err
    is set to an error code describing the reason for the failure.
  */
int64_t htm_tree_s2circle_count(const struct htm_tree *tree,
                                const struct htm_v3 *center,
                                double radius,
                                enum htm_errcode *err);

/** Returns the number of points in \p tree that are inside
    the given spherical ellipse.

    If an error occurs, the return value is negative, and \p *err
    is set to an error code describing the reason for the failure.
  */
int64_t htm_tree_s2ellipse_count(const struct htm_tree *tree,
                                 const struct htm_s2ellipse *ellipse,
                                 enum htm_errcode *err);

/** Returns the number of points in \p tree that are inside
    the given spherical convex polygon.

    If an error occurs, the return value is negative, and \p *err
    is set to an error code describing the reason for the failure.
  */
int64_t htm_tree_s2cpoly_count(const struct htm_tree *tree,
                               const struct htm_s2cpoly *poly,
                               enum htm_errcode *err);

/** Returns a lower and upper bound on the number of points in \p tree
    that are inside the spherical circle with the given center and radius.

    If an error occurs, the returned range will contain an upper
    bound below the lower bound, and \p *err is set to an error code
    describing the reason for the failure.
  */
struct htm_range htm_tree_s2circle_range(const struct htm_tree *tree,
                                         const struct htm_v3 *center,
                                         double radius,
                                         enum htm_errcode *err);

/** Returns a lower and upper bound on the number of points in \p tree
    that are inside the given spherical ellipse.

    If an error occurs, the returned range will contain an upper
    bound below the lower bound, and \p *err is set to an error code
    describing the reason for the failure.
  */
struct htm_range htm_tree_s2ellipse_range(const struct htm_tree *tree,
                                          const struct htm_s2ellipse *ellipse,
                                          enum htm_errcode *err);

/** Returns a lower and upper bound on the number of points in \p tree
    that are inside the given spherical convex polygon.

    If an error occurs, the returned range will contain an upper
    bound below the lower bound, and \p *err is set to an error code
    describing the reason for the failure.
  */
struct htm_range htm_tree_s2cpoly_range(const struct htm_tree *tree,
                                        const struct htm_s2cpoly *poly,
                                        enum htm_errcode *err);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* HTM_TREE_H */

