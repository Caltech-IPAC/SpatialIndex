/** \file
    \brief      Minimalistic HTM indexing implementation.

    This software is based on work by A. Szalay, T. Budavari,
    G. Fekete at The Johns Hopkins University, and Jim Gray,
    Microsoft Research. See the following links for more information:

    http://voservices.net/spherical/
    http://adsabs.harvard.edu/abs/2010PASP..122.1375B

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#include "tinyhtm/htm.h"

#include <stdlib.h>

#include "tinyhtm/tree.h"
#include "tinyhtm/varint.h"


#ifdef __cplusplus
extern "C" {
#endif


/*
    HTM triangles are subdivided into 4 sub-triangles as follows :

            v2
             *
            / \
           /   \
      sv1 *-----* sv0
         / \   / \
        /   \ /   \
    v0 *-----*-----* v1
            sv2

     -  vertices are unit magnitude 3-vectors
     -  edges are great circles on the unit sphere
     -  vertices are stored in counter-clockwise order
       (when viewed from outside the unit sphere in a
       right handed coordinate system)
     -  sv0 = (v1 + v2) / ||v1 + v2||, and likewise for sv1, sv2

    Note that if the HTM triangle given by (v0,v1,v2) has index I, then:
     -  sub triangle T0 = (v0,sv2,sv1) has index I*4
     -  sub triangle T1 = (v1,sv0,sv2) has index I*4 + 1
     -  sub triangle T2 = (v2,sv1,sv0) has index I*4 + 2
     -  sub triangle T3 = (sv0,sv1,sv2) has index I*4 + 3

    All HTM triangles are obtained via subdivision of 8 initial
    triangles, defined from the following set of 6 vertices :
     -  V0 = ( 0,  0,  1) north pole
     -  V1 = ( 1,  0,  0)
     -  V2 = ( 0,  1,  0)
     -  V3 = (-1,  0,  0)
     -  V4 = ( 0, -1,  0)
     -  V5 = ( 0,  0, -1) south pole

    The root triangles (corresponding to subdivision level 0) are :
     -  S0 = (V1, V5, V2), HTM index = 8
     -  S1 = (V2, V5, V3), HTM index = 9
     -  S2 = (V3, V5, V4), HTM index = 10
     -  S3 = (V4, V5, V1), HTM index = 11
     -  N0 = (V1, V0, V4), HTM index = 12
     -  N1 = (V4, V0, V3), HTM index = 13
     -  N2 = (V3, V0, V2), HTM index = 14
     -  N3 = (V2, V0, V1), HTM index = 15

    'S' denotes a triangle in the southern hemisphere,
    'N' denotes a triangle in the northern hemisphere.
 */


/* ---- Types ---- */

/** HTM triangle vs. region classification codes.
  */
enum _htm_cov {
    HTM_DISJOINT = 0,  /**< HTM triangle disjoint from region. */
    HTM_INTERSECT = 1, /**< HTM triangle intersects region. */
    HTM_CONTAINS = 2,  /**< HTM triangle completely contains region. */
    HTM_INSIDE = 3     /**< HTM triangle completely inside region. */
};

/** A node (triangle/trixel) in an HTM tree.
  */
struct _htm_node {
    struct htm_v3 mid_vert[3];      /**< Triangle edge mid-points. */
    struct htm_v3 mid_edge[3];      /**< Subdivision plane normals. */
    const struct htm_v3 *vert[3];   /**< Triangle vertex pointers. */
    const struct htm_v3 *edge[3];   /**< Triangle edge normal pointers. */
    struct htm_v3p *end;            /**< Temporary used for HTM indexing. */
    const unsigned char *s;         /**< Temporary used for tree searches. */
    uint64_t index;                 /**< Temporary used for tree searches. */
    int64_t id;                     /**< HTM ID of the node. */
    int child;                      /**< Index of next child (0-3). */
};

/** A root to leaf path in a depth-first traversal of an HTM tree.
  */
struct _htm_path {
    enum htm_root root; /**< ordinal of root triangle (0-7) */
    struct _htm_node node[HTM_MAX_LEVEL + 1];
};


/* ---- Data ---- */

/*  HTM root triangle vertices/edge plane normals.
 */
static const struct htm_v3 _htm_root_v3[6] = {
  { 0.0,  0.0,  1.0 },
  { 1.0,  0.0,  0.0 },
  { 0.0,  1.0,  0.0 },
  {-1.0,  0.0,  0.0 },
  { 0.0, -1.0,  0.0 },
  { 0.0,  0.0, -1.0 }
};

#define HTM_Z  &_htm_root_v3[0]
#define HTM_X  &_htm_root_v3[1]
#define HTM_Y  &_htm_root_v3[2]
#define HTM_NX &_htm_root_v3[3]
#define HTM_NY &_htm_root_v3[4]
#define HTM_NZ &_htm_root_v3[5]

/*  Vertex pointers for the 3 vertices of each HTM root triangle.
 */
static const struct htm_v3 * const _htm_root_vert[24] = {
  HTM_X,  HTM_NZ, HTM_Y,  /* S0 */
  HTM_Y,  HTM_NZ, HTM_NX, /* S1 */
  HTM_NX, HTM_NZ, HTM_NY, /* S2 */
  HTM_NY, HTM_NZ, HTM_X,  /* S3 */
  HTM_X,  HTM_Z,  HTM_NY, /* N0 */
  HTM_NY, HTM_Z,  HTM_NX, /* N1 */
  HTM_NX, HTM_Z,  HTM_Y,  /* N2 */
  HTM_Y,  HTM_Z,  HTM_X   /* N3 */
};

/*  Edge normal pointers for the 3 edge normals of each HTM root triangle.
 */
static const struct htm_v3 * const _htm_root_edge[24] = {
  HTM_Y,  HTM_X,  HTM_NZ, /* S0 */
  HTM_NX, HTM_Y,  HTM_NZ, /* S1 */
  HTM_NY, HTM_NX, HTM_NZ, /* S2 */
  HTM_X,  HTM_NY, HTM_NZ, /* S3 */
  HTM_NY, HTM_X,  HTM_Z,  /* N0 */
  HTM_NX, HTM_NY, HTM_Z,  /* N1 */
  HTM_Y,  HTM_NX, HTM_Z,  /* N2 */
  HTM_X,  HTM_Y,  HTM_Z   /* N3 */
};


/* ---- Implementation details ---- */

/*  Sets path to the i-th HTM root triangle.
 */
HTM_INLINE void _htm_path_root(struct _htm_path *path, enum htm_root r)
{
    path->node[0].vert[0] = _htm_root_vert[r*3];
    path->node[0].vert[1] = _htm_root_vert[r*3 + 1];
    path->node[0].vert[2] = _htm_root_vert[r*3 + 2];
    path->node[0].edge[0] = _htm_root_edge[r*3];
    path->node[0].edge[1] = _htm_root_edge[r*3 + 1];
    path->node[0].edge[2] = _htm_root_edge[r*3 + 2];
    path->node[0].id = r + 8;
    path->node[0].child = 0;
    path->root = r;
}

/*  Computes the normalized average of two input vertices.
 */
HTM_INLINE void _htm_vertex(struct htm_v3 *out,
                            const struct htm_v3 *v1,
                            const struct htm_v3 *v2)
{
    htm_v3_add(out, v1, v2);
    htm_v3_normalize(out, out);
}

/*  Computes quantities needed by _htm_node_make0(node).
 */
HTM_INLINE void _htm_node_prep0(struct _htm_node *node)
{
    _htm_vertex(&node->mid_vert[1], node->vert[2], node->vert[0]);
    _htm_vertex(&node->mid_vert[2], node->vert[0], node->vert[1]);
    htm_v3_rcross(&node->mid_edge[1], &node->mid_vert[2], &node->mid_vert[1]);
}

/*  Sets node[1] to child 0 of node[0]. Assumes _htm_node_prep0(node)
    has been called.
 */
HTM_INLINE void _htm_node_make0(struct _htm_node *node)
{
    node[1].vert[0] = node[0].vert[0];
    node[1].vert[1] = &node[0].mid_vert[2];
    node[1].vert[2] = &node[0].mid_vert[1];
    node[1].edge[0] = node[0].edge[0];
    node[1].edge[1] = &node[0].mid_edge[1];
    node[1].edge[2] = node[0].edge[2];
    node[0].child = 1;
    node[1].id = node[0].id << 2;
    node[1].child = 0;
}

/*  Computes quantities needed by _htm_node_make1(node). Assumes
    _htm_node_prep0(node) has been called.
 */
HTM_INLINE void _htm_node_prep1(struct _htm_node *node)
{
    _htm_vertex(&node->mid_vert[0], node->vert[1], node->vert[2]);
    htm_v3_rcross(&node->mid_edge[2], &node->mid_vert[0], &node->mid_vert[2]);
}

/*  Sets node[1] to child 1 of node[0]. Assumes _htm_node_prep1(node)
    has been called.
 */
HTM_INLINE void _htm_node_make1(struct _htm_node *node)
{
    node[1].vert[0] = node[0].vert[1];
    node[1].vert[1] = &node[0].mid_vert[0];
    node[1].vert[2] = &node[0].mid_vert[2];
    node[1].edge[0] = node[0].edge[1];
    node[1].edge[1] = &node[0].mid_edge[2];
    node[1].edge[2] = node[0].edge[0];
    node[0].child = 2;
    node[1].id = (node[0].id << 2) + 1;
    node[1].child = 0;
}

/*  Computes quantities needed by _htm_node_make2(node). Assumes
    _htm_node_prep1 has been called.
 */
HTM_INLINE void _htm_node_prep2(struct _htm_node *node)
{
    htm_v3_rcross(&node->mid_edge[0], &node->mid_vert[1], &node->mid_vert[0]);
}

/*  Sets node[1] to child 2 of node[0]. Assumes _htm_node_prep2(node)
    has been called.
 */
HTM_INLINE void _htm_node_make2(struct _htm_node *node)
{
    node[1].vert[0] = node[0].vert[2];
    node[1].vert[1] = &node[0].mid_vert[1];
    node[1].vert[2] = &node[0].mid_vert[0];
    node[1].edge[0] = node[0].edge[2];
    node[1].edge[1] = &node[0].mid_edge[0];
    node[1].edge[2] = node[0].edge[1];
    node[0].child = 3;
    node[1].id = (node[0].id << 2) + 2;
    node[1].child = 0;
}

/*  Sets node[1] to child 3 of node[0]. Assumes _htm_node_prep2(node)
    has been called.
 */
HTM_INLINE void _htm_node_make3(struct _htm_node *node)
{
    htm_v3_neg(&node[0].mid_edge[0], &node[0].mid_edge[0]);
    htm_v3_neg(&node[0].mid_edge[1], &node[0].mid_edge[1]);
    htm_v3_neg(&node[0].mid_edge[2], &node[0].mid_edge[2]);
    node[1].vert[0] = &node[0].mid_vert[0];
    node[1].vert[1] = &node[0].mid_vert[1];
    node[1].vert[2] = &node[0].mid_vert[2];
    node[1].edge[0] = &node[0].mid_edge[0];
    node[1].edge[1] = &node[0].mid_edge[1];
    node[1].edge[2] = &node[0].mid_edge[2];
    node[0].child = 4;
    node[1].id = (node[0].id << 2) + 3;
    node[1].child = 0;
}

/*  Reorders the vectors in [begin,end) such that the resulting array can be
    partitioned into [begin,m) and [m,end), where all vectors in [begin,m) are
    inside the partitioning plane, and all vectors in [m, end) are outside.

    A pointer to the partitioning element m is returned.

    Assumes that plane != 0, begin != 0, end != 0, and begin <= end.
 */
static struct htm_v3p * _htm_partition(const struct htm_v3 *plane,
                                       struct htm_v3p *beg,
                                       struct htm_v3p *end)
{
    struct htm_v3p tmp;
    for (; beg < end; ++beg) {
        if (htm_v3_dot(plane, &beg->v) < 0.0) {
            /* beg is outside plane, find end which is inside,
               swap contents of beg and end. */
            for (--end; end > beg; --end) {
                if (htm_v3_dot(plane, &end->v) >= 0.0) {
                    break;
                }
            }
            if (end <= beg) {
                break;
            }
            tmp = *beg;
            *beg = *end;
            *end = tmp;
        }
    }
    return beg;
}

/*  Perform a depth-first traversal of an HTM tree. At each step of the
    traversal, the input position list is partitioned into the list of
    points inside the current HTM node and those outside - the full
    depth-first traverals of the HTM tree therefore yields a list of
    positions sorted on their HTM indexes. This can be faster than computing
    HTM indexes one at a time when inputs are clustered spatially, since the
    boundary representation of a given HTM triangle is computed at most once.
 */
static void _htm_path_sort(struct _htm_path *path,
                           struct htm_v3p *begin,
                           struct htm_v3p *end,
                           int64_t *ids,
                           int level)
{
    struct _htm_node * const root = path->node;
    struct _htm_node * const leaf = root + level;
    struct _htm_node *curnode = path->node;
    struct htm_v3p *beg = begin;

    curnode->end = end;

    while (1) {
        if (curnode != leaf) {
            /* Not a leaf node, so continue descending the tree.
               Mid-points and edge normals are computed on-demand. */
            int child = curnode->child;
            if (child == 0) {
                _htm_node_prep0(curnode);
                end = _htm_partition(&curnode->mid_edge[1], beg, end);
                if (beg < end) {
                    _htm_node_make0(curnode);
                    ++curnode;
                    curnode->end = end;
                    continue;
                }
                end = curnode->end;
            }
            if (child <= 1) {
                _htm_node_prep1(curnode);
                end = _htm_partition(&curnode->mid_edge[2], beg, end);
                if (beg < end) {
                    _htm_node_make1(curnode);
                    ++curnode;
                    curnode->end = end;
                    continue;
                }
                end = curnode->end;
            }
            if (child <= 2) {
                _htm_node_prep2(curnode);
                end = _htm_partition(&curnode->mid_edge[0], beg, end);
                if (beg < end) {
                    _htm_node_make2(curnode);
                    ++curnode;
                    curnode->end = end;
                    continue;
                }
                end = curnode->end;
            }
            if (beg < end) {
                _htm_node_make3(curnode);
                ++curnode;
                curnode->end = end;
                continue;
            }
        } else {
            /* reached a leaf triangle - all points between beg and end are
               inside and have the same HTM index. */
            int64_t id = curnode->id;
            size_t i = beg - begin;
            size_t j = end - begin;
            for (; i < j; ++i) {
                ids[i] = id;
            }
            beg = end;
        }
        /* walk back up the path until we find a node which
           still contains unsorted points. */
        do {
            if (curnode == root) {
                return;
            }
            --curnode;
            end = curnode->end;
        } while (beg == end);
    }
}

#define HTM_IDS_INIT_CAP 16

static struct htm_ids * _htm_ids_init()
{
    struct htm_ids *ids = (struct htm_ids *) malloc(
        sizeof(struct htm_ids) + HTM_IDS_INIT_CAP * sizeof(struct htm_range));
    if (ids != NULL) {
        ids->n = 0;
        ids->cap = HTM_IDS_INIT_CAP;
    }
    return ids;
}

static struct htm_ids * _htm_ids_grow(struct htm_ids *ids)
{
    size_t cap = ids->cap;
    size_t nbytes = sizeof(struct htm_ids) + 2 * cap * sizeof(struct htm_range);
    struct htm_ids *out = (struct htm_ids *) realloc(ids, nbytes);
    if (out != NULL) {
        out->cap = 2 * cap;
    } else {
        free(ids);
    }
    return out;
}

HTM_INLINE struct htm_ids * _htm_ids_add(struct htm_ids *ids,
                                         int64_t min,
                                         int64_t max)
{
    size_t n = ids->n;
    if (n == 0) {
        ids->n = 1;
        ids->range[0].min = min;
        ids->range[0].max = max;
    } else if (min == ids->range[n - 1].max + 1) {
        ids->range[n - 1].max = max;
    } else {
        if (n == ids->cap) {
            ids = _htm_ids_grow(ids);
            if (ids == NULL) {
                return NULL;
            }
        }
        ids->n = n + 1;
        ids->range[n].min = min;
        ids->range[n].max = max;
    }
    return ids;
}


/*  Returns the coverage code describing the spatial relationship between the
    given HTM triangle and spherical circle.
 */
static enum _htm_cov _htm_s2circle_htmcov(const struct _htm_node *n,
                                          const struct htm_v3 *c,
                                          double dist2)
{
    int nin = htm_v3_dist2(c, n->vert[0]) <= dist2;
    nin += htm_v3_dist2(c, n->vert[1]) <= dist2;
    nin += htm_v3_dist2(c, n->vert[2]) <= dist2;
    if (nin == 3) {
        /* every vertex inside circle */
        return HTM_INSIDE;
    } else if (nin != 0) {
        return HTM_INTERSECT;
    }
    /* no triangle vertices inside circle */
    if (htm_v3_edgedist2(c, n->vert[0], n->vert[1], n->edge[0]) <= dist2 ||
        htm_v3_edgedist2(c, n->vert[1], n->vert[2], n->edge[1]) <= dist2 ||
        htm_v3_edgedist2(c, n->vert[2], n->vert[0], n->edge[2]) <= dist2) {
        /* min distance to at least one edge is <= circle radius */
        return HTM_INTERSECT;
    }
    /* min distance to every edge is > circle radius - circle
       is either inside triangle or disjoint from it */
    if (htm_v3_dot(c, n->edge[0]) >= 0.0 &&
        htm_v3_dot(c, n->edge[1]) >= 0.0 &&
        htm_v3_dot(c, n->edge[2]) >= 0.0) {
        return HTM_CONTAINS;
    }
    return HTM_DISJOINT;
}


/*  Returns 1 if the edge between v1 and v2 intersects the given spherical
    ellipse.

    Let M be the 3x3 symmetric matrix corresponding to the ellipse; the
    ellipse boundary is given by:

    v' M v = 0

    (where ' denotes transpose, and v = [x y z]'). The lines of intersection,
    between the plane defined by v1, v2 and M are given by:

    a*(v1 + v2) + b*(v2 - v1)

    where a,b are scalars to be determined. Note that we can use any basis
    for the plane; (v1 + v2, v2 - v1) is chosen because of numerical stability
    when v1 and v2 are nearly identical, and because it yields a simple test
    for whether a particular solution lies on the edge in question.

    Note that if a,b is a solution, so is k*a,k*b for any scalar k. So
    wlog, set a = 1 and plug into the ellipse boundary equation:

    0 = (v1 + v2)' M (v1 + v2) + b*(v2 - v1)' M (v1 + v2) +
        b*(v1 + v2)' M (v2 - v1) + b^2*(v2 - v1)' M (v2 - v1)

    Since M is symmetric:

    0 = b^2 * (v2 - v1)' M (v2 - v1) +
         2b * (v2 - v1)' M (v1 + v2) +
              (v1 + v2)' M (v1 + v2)

    0 = c22 * b^2 + 2*c21 * b + c11

    Solving this quadratic equation yields 0, 1, or 2 lines of intersection.
    Note that for an intersection to lie on the edge between v1 and v2,
    b must be in the range [-1, 1].
 */
static int _htm_s2ellipse_isect(const struct htm_v3 *v1,
                                const struct htm_v3 *v2,
                                const struct htm_s2ellipse *ellipse)
{
    struct htm_v3 e1, e2, v;
    double c11, c21, c22, delta;

    htm_v3_add(&e1, v1, v2);
    htm_v3_sub(&e2, v2, v1);
    /* compute coeffs of quadratic eqn. */
    c11 = e1.x * e1.x * ellipse->xx +
          e1.y * e1.y * ellipse->yy +
          e1.z * e1.z * ellipse->zz +
          e1.x * e1.y * ellipse->xy * 2.0 +
          e1.x * e1.z * ellipse->xz * 2.0 +
          e1.y * e1.z * ellipse->yz * 2.0;
    c22 = e2.x * e2.x * ellipse->xx +
          e2.y * e2.y * ellipse->yy +
          e2.z * e2.z * ellipse->zz +
          e2.x * e2.y * ellipse->xy * 2.0 +
          e2.x * e2.z * ellipse->xz * 2.0 +
          e2.y * e2.z * ellipse->yz * 2.0;
    c21 = e2.x * e1.x * ellipse->xx +
          e2.y * e1.y * ellipse->yy +
          e2.z * e1.z * ellipse->zz +
          (e2.x * e1.y + e2.y * e1.x) * ellipse->xy +
          (e2.x * e1.z + e2.z * e1.x) * ellipse->xz +
          (e2.y * e1.z + e2.z * e1.y) * ellipse->yz;
    if (c11 == 0.0) {
        /* v1 + v2 is a solution, and lies on the edge */
        if (ellipse->a >= 90.0 || htm_v3_dot(&e1, &ellipse->cen) >= 0.0) {
            return 1;
        }
        /* other solution is given by a linear equation */
        if (c22 == 0.0 || fabs(c22) < fabs(2.0*c21)) {
            return 0;
        }
        /* check whether solution lies in correct hemisphere */
        htm_v3_mul(&v, &e2, -2.0*c21/c22);
        htm_v3_add(&v, &v, &e1);
        return htm_v3_dot(&v, &ellipse->cen) >= 0.0;
    }
    if (c22 == 0.0) {
        /* v2 - v1 is a solution, the other is given by b = -c11/(2*c21). */
        if (c21 == 0.0) {
            return 0;
        }
        if (fabs(c11) <= fabs(2.0*c21)) {
            if (ellipse->a >= 90.0) {
                return 1;
            }
            /* check whether solution lies in correct hemisphere */
            htm_v3_mul(&v, &e2, -0.5*c11/c21);
            htm_v3_add(&v, &v, &e1);
            return htm_v3_dot(&v, &ellipse->cen) >= 0.0;
        }
        return 0;
    }
    delta = c21*c21 - c11*c22;
    if (delta < 0.0) {
        /* no solutions */
        return 0;
    }
    /* 1 or 2 solutions */
    delta = sqrt(delta);
    if (fabs(c22) >= fabs(delta - c21)) {
        if (ellipse->a >= 90.0) {
            return 1;
        }
        /* check whether solution lies in correct hemisphere */
        htm_v3_mul(&v, &e2, (delta - c21)/c22);
        htm_v3_add(&v, &v, &e1);
        return htm_v3_dot(&v, &ellipse->cen) >= 0.0;
    }
    if (fabs(c22) >= fabs(delta + c21)) {
        if (ellipse->a >= 90.0) {
            return 1;
        }
        /* check whether solution lies in correct hemisphere */
        htm_v3_mul(&v, &e2, -(delta + c21)/c22);
        htm_v3_add(&v, &v, &e1);
        return htm_v3_dot(&v, &ellipse->cen) >= 0.0;
    }
    return 0;
}


/*  Returns the coverage code describing the spatial relationship between the
    given HTM triangle and spherical ellipse.
 */
static enum _htm_cov _htm_s2ellipse_htmcov(const struct _htm_node *n,
                                           const struct htm_s2ellipse *e)
{
    int nin = htm_s2ellipse_cv3(e, n->vert[0]);
    nin += htm_s2ellipse_cv3(e, n->vert[1]);
    nin += htm_s2ellipse_cv3(e, n->vert[2]);
    if (nin == 3) {
        return HTM_INSIDE;
    } else if (nin != 0) {
        return HTM_INTERSECT;
    }
    /* no triangle vertices inside ellipse - check for edge/ellipse
       intersections */
    if (_htm_s2ellipse_isect(n->vert[0], n->vert[1], e) ||
        _htm_s2ellipse_isect(n->vert[1], n->vert[2], e) ||
        _htm_s2ellipse_isect(n->vert[2], n->vert[0], e)) {
        return HTM_INTERSECT;
    }
    /* no triangle/ellipse intersections */
    if (htm_v3_dot(&e->cen, n->edge[0]) >= 0.0 &&
        htm_v3_dot(&e->cen, n->edge[1]) >= 0.0 &&
        htm_v3_dot(&e->cen, n->edge[2]) >= 0.0) {
        /* ellipse center inside triangle */
        return HTM_CONTAINS;
    }
    return HTM_DISJOINT;
}


static const double HTM_INF = 1.0 / 0.0;
static const double HTM_NEG_INF = -1.0 / 0.0;

/*  Tests whether poly intersects the edge (v1, v2) with plane normal n.

    The idea is that a solution v = (x,y,z) must satisfy:

        v . n = 0, v != 0
        v . (n ^ v1) >= 0
        v . (v2 ^ n) >= 0
        v . e_i >= 0

    where e_i are the edge plane normals for the polygon, and (n ^ v1),
    (v2 ^ n) are plane normals that bound the lune defined by n, v1, and v2.
    Write this as:

        v . n = 0
        v . c_i >= 0

    Now assume nz > 0 (for a non zero solution, at least one of nx, ny, nz
    must be non-zero, and treatment of negative values is analagous). Use
    the equality to obtain

        z = - (x * nx + y * ny) / nz

    and substitute into the inequalities to obtain:

        x * (c_ix * nz - c_iz * nx) + y * (c_iy * nz - c_iz * ny) >= 0

    which we write

        x * a_i + y * b_i >= 0

    If a solution v exists, then Kv is also a solution (for positive scalar K),
    so we can fix y = 1 and look for solutions to

        x * a_i + b_i >= 0

    If there are none, fix y = -1 and look for solutions to:

        x * a_i - b_i >= 0

    If again there are none, then y = 0, and the problem reduces to checking
    whether

        x * a_i >= 0

    has any solutions. This is the case when the non-zero a_i have the same
    sign.
 */
static int _htm_isect_test(const struct htm_v3 *v1,
                           const struct htm_v3 *v2,
                           const struct htm_v3 *n,
                           const struct htm_s2cpoly *poly,
                           double *ab)
{
    struct htm_v3 c0;
    struct htm_v3 c1;
    double min_1, max_1, min_m1, max_m1;
    size_t nv, i, neg, pos;

    htm_v3_cross(&c0, n, v1);
    htm_v3_cross(&c1, v2, n);
    nv = poly->n;
    if (n->z != 0.0) {
        double s = (n->z > 0.0) ? 1.0 : -1.0;
        ab[0] = s * (c0.x * n->z - c0.z * n->x);
        ab[1] = s * (c0.y * n->z - c0.z * n->y);
        ab[2] = s * (c1.x * n->z - c1.z * n->x);
        ab[3] = s * (c1.y * n->z - c1.z * n->y);
        for (i = 0; i < nv; ++i) {
            ab[2*i + 4] = s * (poly->ve[nv + i].x * n->z - poly->ve[nv + i].z * n->x);
            ab[2*i + 5] = s * (poly->ve[nv + i].y * n->z - poly->ve[nv + i].z * n->y);
        }
    } else if (n->y != 0.0) {
        double s = (n->y > 0.0) ? 1.0 : -1.0;
        ab[0] = s * (c0.x * n->y - c0.y * n->x);
        ab[1] = s * (c0.z * n->y);
        ab[2] = s * (c1.x * n->y - c1.y * n->x);
        ab[3] = s * (c1.z * n->y);
        for (i = 0; i < nv; ++i) {
            ab[2*i + 4] = s * (poly->ve[nv + i].x * n->y - poly->ve[nv + i].y * n->x);
            ab[2*i + 5] = s * (poly->ve[nv + i].z * n->y);
        }
    } else if (n->x != 0.0) {
        double s = (n->x > 0.0) ? 1.0 : -1.0;
        ab[0] = s * (c0.y * n->x);
        ab[1] = s * (c0.z * n->x);
        ab[2] = s * (c1.y * n->x);
        ab[3] = s * (c1.z * n->x);
        for (i = 0; i < nv; ++i) {
            ab[2*i + 4] = s * (poly->ve[nv + i].y * n->x);
            ab[2*i + 5] = s * (poly->ve[nv + i].z * n->x);
        }
    } else {
        return 0;
    }
    /* search for solutions to a*x +/- b >= 0, with constraint coeffs stored in
       ab */
    min_1 = min_m1 = HTM_NEG_INF;
    max_1 = max_m1 = HTM_INF;
    for (i = 0, neg = 0, pos = 0; i < nv + 2; ++i) {
        double a = ab[2*i];
        double b = ab[2*i + 1];
        if (a == 0.0) {
            if (b < 0.0) {
                min_1 = HTM_INF;
                max_1 = HTM_NEG_INF;
            } else if (b > 0.0) {
                min_m1 = HTM_INF;
                max_m1 = HTM_NEG_INF;
            }
        } else if (a < 0.0) {
            ++neg;
            double d = -b / a;
            if (d < max_1) {
                max_1 = d;
            }
            if (-d < max_m1) {
                max_m1 = -d;
            }
        } else {
            ++pos;
            double d = -b / a;
            if (d > min_1) {
                min_1 = d;
            }
            if (-d > min_m1) {
                min_m1 = -d;
            }
        }
    }
    if (min_1 <= max_1 || min_m1 <= max_m1) {
        return 1;
    }
    return (neg == 0 || pos == 0);
}

/*  Returns the coverage code describing the spatial relationship between the
    given HTM triangle and spherical convex polygon.
 */
static enum _htm_cov _htm_s2cpoly_htmcov(const struct _htm_node *n,
                                         const struct htm_s2cpoly *poly,
                                         double *ab)
{
    int nin = htm_s2cpoly_cv3(poly, n->vert[0]);
    nin += htm_s2cpoly_cv3(poly, n->vert[1]);
    nin += htm_s2cpoly_cv3(poly, n->vert[2]);

    if (nin == 3) {
        /* all triangle vertices are inside poly,
           so triangle is inside by convexity. */
        return HTM_INSIDE;
    } else if (nin != 0) {
        return HTM_INTERSECT;
    }
    /* all triangle vertices are outside poly - check for edge intersections */
    if (_htm_isect_test(n->vert[0], n->vert[1], n->edge[0], poly, ab) != 0 ||
        _htm_isect_test(n->vert[1], n->vert[2], n->edge[1], poly, ab) != 0 ||
        _htm_isect_test(n->vert[2], n->vert[0], n->edge[2], poly, ab) != 0) {
        return HTM_INTERSECT;
    }
    /* All triangle vertices are outside poly and there are no edge/edge
       intersections. Polygon is either inside triangle or disjoint from
       it */
    if (htm_v3_dot(&poly->vsum, n->edge[0]) >= 0.0 &&
        htm_v3_dot(&poly->vsum, n->edge[1]) >= 0.0 &&
        htm_v3_dot(&poly->vsum, n->edge[2]) >= 0.0) {
        return HTM_CONTAINS;
    }
    return HTM_DISJOINT;
}

/*  Returns the HTM root triangle for a point.
 */
HTM_INLINE enum htm_root _htm_v3_htmroot(const struct htm_v3 *v)
{
    if (v->z < 0.0) {
        /* S0, S1, S2, S3 */
        if (v->y > 0.0) {
            return (v->x > 0.0) ? HTM_S0 : HTM_S1;
        } else if (v->y == 0.0) {
            return (v->x >= 0.0) ? HTM_S0 : HTM_S2;
        } else {
            return (v->x < 0.0) ? HTM_S2 : HTM_S3;
        }
    } else {
        /* N0, N1, N2, N3 */
        if (v->y > 0.0) {
            return (v->x > 0.0) ? HTM_N3 : HTM_N2;
        } else if (v->y == 0.0) {
            return (v->x >= 0.0) ? HTM_N3 : HTM_N1;
        } else {
            return (v->x < 0.0) ? HTM_N1 : HTM_N0;
        }
    }
}

/*  Partitions an array of points according to their root triangle numbers.
 */
static size_t _htm_rootpart(struct htm_v3p *points,
                            unsigned char *ids,
                            size_t n,
                            enum htm_root root)
{
    struct htm_v3p tmp;
    size_t beg, end;
    unsigned char c;
    for (beg = 0, end = n; beg < end; ++beg) {
        if (ids[beg] >= root) {
            for (; end > beg; --end) {
                if (ids[end - 1] < root) {
                    break;
                }
            }
            if (end == beg) {
                break;
            }
            tmp = points[beg];
            points[beg] = points[end - 1];
            points[end - 1] = tmp;
            c = ids[beg];
            ids[beg] = ids[end - 1];
            ids[end - 1] = c;
        }
    }
    return beg;
}

/*  Sorts the given array of positions by root triangle number.
 */
static void _htm_rootsort(size_t roots[HTM_NROOTS + 1],
                          struct htm_v3p *points,
                          unsigned char *ids,
                          size_t n)
{
    size_t i, n0, n2, s2;

    /* compute root ids for all points */
    for (i = 0; i < n; ++i) {
        ids[i] = (unsigned char) _htm_v3_htmroot(&points[i].v);
    }
    n0 = _htm_rootpart(points, ids, n, HTM_N0);
    s2 = _htm_rootpart(points, ids, n0, HTM_S2);
    roots[HTM_S0] = 0;
    roots[HTM_S1] = _htm_rootpart(points, ids, s2, HTM_S1);
    roots[HTM_S2] = s2;
    roots[HTM_S3] = _htm_rootpart(points + s2, ids + s2, n0 - s2, HTM_S3) + s2;
    n2 = _htm_rootpart(points + n0, ids + n0, n - n0, HTM_N2) + n0;
    roots[HTM_N0] = n0;
    roots[HTM_N1] = _htm_rootpart(points + n0, ids + n0, n2 - n0, HTM_N1) + n0;
    roots[HTM_N2] = n2;
    roots[HTM_N3] = _htm_rootpart(points + n2, ids + n2, n - n2, HTM_N3) + n2;
    roots[HTM_NROOTS] = n;
}

/*  Reduces the effective subdivision level of an HTM id range list by n levels
    and merges adjacent ranges. This typically reduces the number of ranges in
    the list, but also makes it a poorer approximation of the underlying
    geometry. Note that with sufficiently large n, any range list can be shrunk
    down to at most 4 ranges.

    In detail: a range [I1, I2] is mapped to [I1 & ~mask, I2 | mask], where
    mask = (1 << 2*n) - 1.
 */
static void _htm_simplify_ids(struct htm_ids *ids, int n)
{
    size_t i, j, nr;
    int64_t mask;
    if (n <= 0 || ids == 0 || ids->n == 0) {
        return;
    }
    mask = (((int64_t) 1) << 2*n) - 1;
    for (i = 0, j = 0, nr = ids->n; i < nr; ++i, ++j) {
        int64_t min = ids->range[i].min & ~mask;
        int64_t max = ids->range[i].max | mask;
        for (; i < nr - 1; ++i) {
            int64_t next = ids->range[i + 1].min & ~mask;
            if (next > max + 1) {
                break;
            }
            max = ids->range[i + 1].max | mask;
        }
        ids->range[j].min = min;
        ids->range[j].max = max;
    }
    ids->n = j;
}


/*  Constructs the next non-empty child of \p node.
 */
static const unsigned char * _htm_subdivide(struct _htm_node *node,
                                            const unsigned char *s)
{
    uint64_t off;
    switch (node->child) {
        case 0:
            off = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            _htm_node_prep0(node);
            _htm_node_make0(node);
            if (off != 0) {
                break;
            }
            /* fall-through */
        case 1:
            off = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            _htm_node_prep1(node);
            _htm_node_make1(node);
            if (off != 0) {
                break;
            }
            /* fall-through */
        case 2:
            off = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            _htm_node_prep2(node);
            _htm_node_make2(node);
            if (off != 0) {
                break;
            }
            /* fall-through */
        case 3:
            off = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            if (off != 0) {
                _htm_node_make3(node);
                break;
            }
        default:
            return NULL;
    }
    node->s = s;
    return s + (off - 1);
}


/* ---- API ---- */

int64_t htm_v3_id(const struct htm_v3 *point, int level)
{
    struct htm_v3 v0, v1, v2;
    struct htm_v3 sv0, sv1, sv2;
    struct htm_v3 e;
    int64_t id;
    int curlevel;
    enum htm_root r;

    if (point == NULL) {
        return 0;
    }
    if (level < 0 || level > HTM_MAX_LEVEL) {
        return 0;
    }
    r = _htm_v3_htmroot(point);
    v0 = *_htm_root_vert[r*3];
    v1 = *_htm_root_vert[r*3 + 1];
    v2 = *_htm_root_vert[r*3 + 2];
    id = r + 8;
    for (curlevel = 0; curlevel < level; ++curlevel) {
        _htm_vertex(&sv1, &v2, &v0);
        _htm_vertex(&sv2, &v0, &v1);
        htm_v3_rcross(&e, &sv2, &sv1);
        if (htm_v3_dot(&e, point) >= 0) {
            v1 = sv2;
            v2 = sv1;
            id = id << 2;
            continue;
        }
        _htm_vertex(&sv0, &v1, &v2);
        htm_v3_rcross(&e, &sv0, &sv2);
        if (htm_v3_dot(&e, point) >= 0) {
            v0 = v1;
            v1 = sv0;
            v2 = sv2;
            id = (id << 2) + 1;
            continue;
        }
        htm_v3_rcross(&e, &sv1, &sv0);
        if (htm_v3_dot(&e, point) >= 0) {
            v0 = v2;
            v1 = sv1;
            v2 = sv0;
            id = (id << 2) + 2;
        } else {
            v0 = sv0;
            v1 = sv1;
            v2 = sv2;
            id = (id << 2) + 3;
        }
    }
    return id;
}


enum htm_errcode htm_v3p_idsort(struct htm_v3p *points,
                                int64_t *ids,
                                size_t n,
                                int level)
{
    struct _htm_path path;
    size_t roots[HTM_NROOTS + 1];
    enum htm_root r;

    if (n == 0) {
        return HTM_ELEN;
    } else if (points == 0 || ids == 0) {
        return HTM_ENULLPTR;
    } else if (level < 0 || level > HTM_MAX_LEVEL) {
        return HTM_ELEVEL;
    }
    _htm_rootsort(roots, points, (unsigned char *) ids, n);
    for (r = HTM_S0; r <= HTM_N3; ++r) {
        if (roots[r] < roots[r + 1]) {
            _htm_path_root(&path, r);
            _htm_path_sort(&path, points + roots[r],
                                 points + roots[r + 1], ids + roots[r], level);
        }
    }
    return HTM_OK;
}


int htm_level(int64_t id)
{
    uint64_t x = (uint64_t) id;
    int l;
    if (id < 8) {
        return -1;
    }
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (x >> 32);
    l = htm_popcount(x) - 4;
    /* check that l is even, in range, and that the 4 MSBs of id
       give a valid root ID (8-15) */
    if ((l & 1) != 0 || ((id >> l) & 0x8) == 0 || l > HTM_MAX_LEVEL*2) {
        return -1;
    }
    return l / 2;
}


enum htm_errcode htm_tri_init(struct htm_tri *tri, int64_t id)
{
    struct htm_v3 v0, v1, v2;
    struct htm_v3 sv0, sv1, sv2;
    int shift, level;
    enum htm_root r;

    if (tri == NULL) {
        return HTM_ENULLPTR;
    }
    level = htm_level(id);
    if (level < 0) {
        return HTM_EID;
    }
    tri->id = id;
    tri->level = level;
    shift = 2*level;
    r = (id >> shift) & 0x7;
    v0 = *_htm_root_vert[r*3];
    v1 = *_htm_root_vert[r*3 + 1];
    v2 = *_htm_root_vert[r*3 + 2];
    for (shift -= 2; shift >= 0; shift -= 2) {
        int child = (id >> shift) & 0x3;
        _htm_vertex(&sv1, &v2, &v0);
        _htm_vertex(&sv2, &v0, &v1);
        _htm_vertex(&sv0, &v1, &v2);
        switch (child) {
            case 0:
                v1 = sv2;
                v2 = sv1;
                break;
            case 1:
                v0 = v1;
                v1 = sv0;
                v2 = sv2;
                break;
            case 2:

                v0 = v2;
                v1 = sv1;
                v2 = sv0;
                break;
            case 3:
                v0 = sv0;
                v1 = sv1;
                v2 = sv2;
                break;
        }
    }
    tri->verts[0] = v0;
    tri->verts[1] = v1;
    tri->verts[2] = v2;
    htm_v3_add(&sv0, &v0, &v1);
    htm_v3_add(&sv0, &sv0, &v2);
    htm_v3_normalize(&tri->center, &sv0);
    tri->radius = htm_v3_angsep(&sv0, &v0);
    return HTM_OK;
}


struct htm_ids * htm_s2circle_ids(struct htm_ids *ids,
                                  const struct htm_v3 *center,
                                  double radius,
                                  int level,
                                  size_t maxranges,
                                  enum htm_errcode *err)
{
    struct _htm_path path;
    double dist2;
    enum htm_root root;
    int efflevel;

    if (center == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        free(ids);
        return NULL;
    } else if (level < 0 || level > HTM_MAX_LEVEL) {
        if (err != NULL) {
            *err = HTM_ELEVEL;
        }
        free(ids);
        return NULL;
    }
    if (ids == NULL) {
        ids = _htm_ids_init();
        if (ids == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            return NULL;
        }
    } else {
        ids->n = 0;
    }
    /* Deal with degenerate cases */
    if (radius < 0.0) {
        /* empty ID list */
        if (err != NULL) {
            *err = HTM_OK;
        }
        return ids;
    } else if (radius >= 180.0) {
        /* the entire sky */
        int64_t min_id = (8 + HTM_S0) << level * 2;
        int64_t max_id = ((8 + HTM_NROOTS) << level * 2) - 1;
        if (err != NULL) {
            *err = HTM_OK;
        }
        ids = _htm_ids_add(ids, min_id, max_id);
        if (ids == NULL && err != NULL) {
            *err = HTM_ENOMEM;
        }
        return ids;
    }

    efflevel = level;
    /* compute square of secant distance corresponding to radius */
    dist2 = sin(radius * 0.5 * HTM_RAD_PER_DEG);
    dist2 = 4.0 * dist2 * dist2;

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        int curlevel = 0;
        _htm_path_root(&path, root);

        while (1) {
            switch (_htm_s2circle_htmcov(curnode, center, dist2)) {
                case HTM_CONTAINS:
                    if (curlevel == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (curlevel < efflevel) {
                        /* continue subdividing */
                        _htm_node_prep0(curnode);
                        _htm_node_make0(curnode);
                        ++curnode;
                        ++curlevel;
                        continue;
                    }
                    /* fall-through */
                case HTM_INSIDE:
                    /* reached a leaf or fully covered HTM triangle,
                       append HTM ID range to results */
                    {
                        int64_t id = curnode->id << (level - curlevel) * 2;
                        int64_t n = ((int64_t) 1) << (level - curlevel) * 2;
                        ids = _htm_ids_add(ids, id, id + n - 1);
                    }
                    if (ids == NULL) {
                        if (err != NULL) {
                            *err = HTM_ENOMEM;
                        }
                        return NULL;
                    }
                    while (ids->n > maxranges && efflevel != 0) {
                        /* too many ranges:
                           reduce effective subdivision level */
                        --efflevel;
                        if (curlevel > efflevel) {
                           curnode = curnode - (curlevel - efflevel);
                           curlevel = efflevel;
                        }
                        _htm_simplify_ids(ids, level - efflevel);
                    }
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
            --curlevel;
            --curnode;
            while (curlevel >= 0 && curnode->child == 4) {
                --curnode;
                --curlevel;
            }
            if (curlevel < 0) {
                /* finished with this root */
                break;
            }
            if (curnode->child == 1) {
                _htm_node_prep1(curnode);
                _htm_node_make1(curnode);
            } else if (curnode->child == 2) {
                _htm_node_prep2(curnode);
                _htm_node_make2(curnode);
            } else {
                _htm_node_make3(curnode);
            }
            ++curnode;
            ++curlevel;
        }
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return ids;
}


struct htm_ids * htm_s2ellipse_ids(struct htm_ids *ids,
                                   const struct htm_s2ellipse *ellipse,
                                   int level,
                                   size_t maxranges,
                                   enum htm_errcode *err)
{
    struct _htm_path path;
    enum htm_root root;
    int efflevel;

    if (ellipse == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        free(ids);
        return NULL;
    } else if (level < 0 || level > HTM_MAX_LEVEL) {
        if (err != NULL) {
            *err = HTM_ELEVEL;
        }
        free(ids);
        return NULL;
    }
    if (ids == NULL) {
        ids = _htm_ids_init();
        if (ids == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            return NULL;
        }
    } else {
        ids->n = 0;
    }

    efflevel = level;
    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        int curlevel = 0;
        _htm_path_root(&path, root);

        while (1) {
            switch (_htm_s2ellipse_htmcov(curnode, ellipse)) {
                case HTM_CONTAINS:
                    if (curlevel == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (curlevel < efflevel) {
                        /* continue subdividing */
                        _htm_node_prep0(curnode);
                        _htm_node_make0(curnode);
                        ++curnode;
                        ++curlevel;
                        continue;
                    }
                    /* fall-through */
                case HTM_INSIDE:
                    /* reached a leaf or fully covered HTM triangle,
                       append HTM ID range to results */
                    {
                        int64_t id = curnode->id << (level - curlevel) * 2;
                        int64_t n = ((int64_t) 1) << (level - curlevel) * 2;
                        ids = _htm_ids_add(ids, id, id + n - 1);
                    }
                    if (ids == NULL) {
                        if (err != NULL) {
                            *err = HTM_ENOMEM;
                        }
                        return NULL;
                    }
                    while (ids->n > maxranges && efflevel != 0) {
                        /* too many ranges:
                           reduce effective subdivision level */
                        --efflevel;
                        if (curlevel > efflevel) {
                           curnode = curnode - (curlevel - efflevel);
                           curlevel = efflevel;
                        }
                        _htm_simplify_ids(ids, level - efflevel);
                    }
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
            --curlevel;
            --curnode;
            while (curlevel >= 0 && curnode->child == 4) {
                --curnode;
                --curlevel;
            }
            if (curlevel < 0) {
                /* finished with this root */
                break;
            }
            if (curnode->child == 1) {
                _htm_node_prep1(curnode);
                _htm_node_make1(curnode);
            } else if (curnode->child == 2) {
                _htm_node_prep2(curnode);
                _htm_node_make2(curnode);
            } else {
                _htm_node_make3(curnode);
            }
            ++curnode;
            ++curlevel;
        }
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return ids;
}


struct htm_ids * htm_s2cpoly_ids(struct htm_ids *ids,
                                 const struct htm_s2cpoly *poly,
                                 int level,
                                 size_t maxranges,
                                 enum htm_errcode *err)
{
    double stackab[2*256 + 4];
    struct _htm_path path;
    enum htm_root root;
    double *ab;
    size_t nb;
    int efflevel;

    if (poly == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        free(ids);
        return NULL;
    } else if (level < 0 || level > HTM_MAX_LEVEL) {
        if (err != NULL) {
            *err = HTM_ELEVEL;
        }
        free(ids);
        return NULL;
    }
    if (ids == NULL) {
        ids = _htm_ids_init();
        if (ids == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            return NULL;
        }
    } else {
        ids->n = 0;
    }
    nb = (2 * poly->n + 4) * sizeof(double);
    if (nb > sizeof(stackab)) {
        ab = (double *) malloc(nb);
        if (ab == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            free(ids);
            return NULL;
        }
    } else {
        ab = stackab;
    }

    efflevel = level;

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        int curlevel = 0;
        _htm_path_root(&path, root);

        while (1) {
            switch (_htm_s2cpoly_htmcov(curnode, poly, ab)) {
                case HTM_CONTAINS:
                    if (curlevel == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (curlevel < efflevel) {
                        /* continue subdividing */
                        _htm_node_prep0(curnode);
                        _htm_node_make0(curnode);
                        ++curnode;
                        ++curlevel;
                        continue;
                    }
                    /* fall-through */
                case HTM_INSIDE:
                    /* reached a leaf or fully covered HTM triangle,
                       append HTM ID range to results */
                    {
                        int64_t id = curnode->id << (level - curlevel) * 2;
                        int64_t n = ((int64_t) 1) << (level - curlevel) * 2;
                        ids = _htm_ids_add(ids, id, id + n - 1);
                    }
                    if (ids == NULL) {
                        if (ab != stackab) {
                            free(ab);
                        }
                        if (err != NULL) {
                            *err = HTM_ENOMEM;
                        }
                        return ids;
                    }
                    while (ids->n > maxranges && efflevel != 0) {
                        /* too many ranges:
                           reduce effetive subdivision level */
                        --efflevel;
                        if (curlevel > efflevel) {
                           curnode = curnode - (curlevel - efflevel);
                           curlevel = efflevel;
                        }
                        _htm_simplify_ids(ids, level - efflevel);
                    }
                    break;
                default:
                    /* HTM triangle does not intersect polygon */
                    break;
            }
            /* ascend towards the root */
            --curlevel;
            --curnode;
            while (curlevel >= 0 && curnode->child == 4) {
                --curnode;
                --curlevel;
            }
            if (curlevel < 0) {
                /* finished with this root */
                break;
            }
            if (curnode->child == 1) {
                _htm_node_prep1(curnode);
                _htm_node_make1(curnode);
            } else if (curnode->child == 2) {
                _htm_node_prep2(curnode);
                _htm_node_make2(curnode);
            } else {
                _htm_node_make3(curnode);
            }
            ++curnode;
            ++curlevel;
        }
    }
    if (ab != stackab) {
        free(ab);
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return ids;
}


int64_t htm_idtodec(int64_t id)
{
    int64_t dec = 0;
    int64_t factor = 1;
    int level = htm_level(id);
    if (level < 0 || level > HTM_DEC_MAX_LEVEL) {
        return 0;
    }
    for (++level; level > 0; --level, id >>= 2, factor *= 10) {
        dec += factor * (id & 3);
    }
    if ((id & 1) == 1) {
        dec += 2*factor;
    } else {
        dec += factor;
    }
    return dec;
}


int64_t htm_tree_s2circle_count(const struct htm_tree *tree,
                                const struct htm_v3 *center,
                                double radius,
                                enum htm_errcode *err)
{
    struct _htm_path path;
    enum htm_root root;
    double d2;
    int64_t count;

    if (tree == NULL || center == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return -1;
    }
    if (tree->indexfd == -1) {
        return htm_tree_s2circle_scan(tree, center, radius, err);
    } else if (radius < 0.0) {
        /* circle is empty */
        return 0;
    } else if (radius >= 180.0) {
        /* entire sky */
        return (int64_t) tree->count;
    }

    count = 0;
    /* compute square of secant distance corresponding to radius */
    d2 = sin(radius * 0.5 * HTM_RAD_PER_DEG);
    d2 = 4.0 * d2 * d2;

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        const unsigned char *s = tree->root[root];
        uint64_t index = 0;
        int level = 0;

        if (s == NULL) {
            /* root contains no points */
            continue;
        }
        _htm_path_root(&path, root);

        while (1) {
            uint64_t curcount = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            index += htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            curnode->index = index;

            switch (_htm_s2circle_htmcov(curnode, center, d2)) {
                case HTM_CONTAINS:
                    if (level == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (level < 20 && curcount >= tree->leafthresh) {
                        s = _htm_subdivide(curnode, s);
                        if (s == NULL) {
                            /* tree is invalid */
                            if (err != NULL) {
                                *err = HTM_EINV;
                            }
                            return -1;
                        }
                        ++level;
                        ++curnode;
                        continue;
                    }
                    /* scan points in leaf */
                    {
                        uint64_t i;
                        for (i = index; i < index + curcount; ++i) {
                            if (htm_v3_dist2(center, &tree->entries[i].v) <= d2) {
                                ++count;
                            }
                        }
                    }
                    break;
                case HTM_INSIDE:
                    /* fully covered HTM triangle */
                    count += (int64_t) curcount;
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
ascend:
            --level;
            --curnode;
            while (level >= 0 && curnode->child == 4) {
                --curnode;
                --level;
            }
            if (level < 0) {
                /* finished with this root */
                break;
            }
            index = curnode->index;
            s = _htm_subdivide(curnode, curnode->s);
            if (s == NULL) {
                /* no non-empty children remain */
                goto ascend;
            }
            ++level;
            ++curnode;
        }
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return count;
}


int64_t htm_tree_s2ellipse_count(const struct htm_tree *tree,
                                 const struct htm_s2ellipse *ellipse,
                                 enum htm_errcode *err)
{
    struct _htm_path path;
    enum htm_root root;
    int64_t count;

    if (tree == NULL || ellipse == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return -1;
    }
    if (tree->indexfd == -1) {
        return htm_tree_s2ellipse_scan(tree, ellipse, err);
    }

    count = 0;

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        const unsigned char *s = tree->root[root];
        uint64_t index = 0;
        int level = 0;

        if (s == NULL) {
            /* root contains no points */
            continue;
        }
        _htm_path_root(&path, root);

        while (1) {
            uint64_t curcount = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            index += htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            curnode->index = index;

            switch (_htm_s2ellipse_htmcov(curnode, ellipse)) {
                case HTM_CONTAINS:
                    if (level == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (level < 20 && curcount >= tree->leafthresh) {
                        s = _htm_subdivide(curnode, s);
                        if (s == NULL) {
                            /* tree is invalid */
                            if (err != NULL) {
                                *err = HTM_EINV;
                            }
                            return -1;
                        }
                        ++level;
                        ++curnode;
                        continue;
                    }
                    /* scan points in leaf */
                    {
                        uint64_t i;
                        for (i = index; i < index + curcount; ++i) {
                            if (htm_s2ellipse_cv3(ellipse, &tree->entries[i].v)) {
                                ++count;
                            }
                        }
                    }
                    break;
                case HTM_INSIDE:
                    /* fully covered HTM triangle */
                    count += (int64_t) curcount;
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
ascend:
            --level;
            --curnode;
            while (level >= 0 && curnode->child == 4) {
                --curnode;
                --level;
            }
            if (level < 0) {
                /* finished with this root */
                break;
            }
            index = curnode->index;
            s = _htm_subdivide(curnode, curnode->s);
            if (s == NULL) {
                /* no non-empty children remain */
                goto ascend;
            }
            ++level;
            ++curnode;
        }
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return count;
}


int64_t htm_tree_s2cpoly_count(const struct htm_tree *tree,
                               const struct htm_s2cpoly *poly,
                               enum htm_errcode *err)
{
    double stackab[2*256 + 4];
    struct _htm_path path;
    enum htm_root root;
    double *ab;
    size_t nb;
    int64_t count;

    if (tree == NULL || poly == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return -1;
    }
    if (tree->indexfd == -1) {
        return htm_tree_s2cpoly_scan(tree, poly, err);
    }
    nb = (2 * poly->n + 4) * sizeof(double);
    if (nb > sizeof(stackab)) {
        ab = (double *) malloc(nb);
        if (ab == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            return -1;
        }
    } else {
        ab = stackab;
    }
    count = 0;

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        const unsigned char *s = tree->root[root];
        uint64_t index = 0;
        int level = 0;

        if (s == NULL) {
            /* root contains no points */
            continue;
        }
        _htm_path_root(&path, root);

        while (1) {
            uint64_t curcount = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            index += htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            curnode->index = index;

            switch (_htm_s2cpoly_htmcov(curnode, poly, ab)) {
                case HTM_CONTAINS:
                    if (level == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (level < 20 && curcount >= tree->leafthresh) {
                        s = _htm_subdivide(curnode, s);
                        if (s == NULL) {
                            /* tree is invalid */
                            if (ab != stackab) {
                                free(ab);
                            }
                            if (err != NULL) {
                                *err = HTM_EINV;
                            }
                            return -1;
                        }
                        ++level;
                        ++curnode;
                        continue;
                    }
                    /* scan points in leaf */
                    {
                        uint64_t i;
                        for (i = index; i < index + curcount; ++i) {
                            if (htm_s2cpoly_cv3(poly, &tree->entries[i].v)) {
                                ++count;
                            }
                        }
                    }
                    break;
                case HTM_INSIDE:
                    /* fully covered HTM triangle */
                    count += (int64_t) curcount;
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
ascend:
            --level;
            --curnode;
            while (level >= 0 && curnode->child == 4) {
                --curnode;
                --level;
            }
            if (level < 0) {
                /* finished with this root */
                break;
            }
            index = curnode->index;
            s = _htm_subdivide(curnode, curnode->s);
            if (s == NULL) {
                /* no non-empty children remain */
                goto ascend;
            }
            ++level;
            ++curnode;
        }
    }
    if (ab != stackab) {
        free(ab);
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return count;
}


struct htm_range htm_tree_s2circle_range(const struct htm_tree *tree,
                                         const struct htm_v3 *center,
                                         double radius,
                                         enum htm_errcode *err)
{
    struct _htm_path path;
    enum htm_root root;
    double d2;
    struct htm_range range;

    range.min = 0;
    range.max = 0;
    if (tree == NULL || center == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        range.max = -1;
        return range;
    }
    if (tree->indexfd == -1) {
        range.max = htm_tree_s2circle_scan(tree, center, radius, err);
        if (range.max >= 0) {
            range.min = range.max;
        }
        return range;
    } else if (radius < 0.0) {
        /* circle is empty */
        return range;
    } else if (radius >= 180.0) {
        /* entire sky */
        range.min = (int64_t) tree->count;
        range.max = (int64_t) tree->count;
        return range;
    }
    /* compute square of secant distance corresponding to radius */
    d2 = sin(radius * 0.5 * HTM_RAD_PER_DEG);
    d2 = 4.0 * d2 * d2;

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        const unsigned char *s = tree->root[root];
        int level = 0;

        if (s == NULL) {
            /* root contains no points */
            continue;
        }
        _htm_path_root(&path, root);

        while (1) {
            uint64_t curcount = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            s += 1 + htm_varint_nfollow(*s);
            switch (_htm_s2circle_htmcov(curnode, center, d2)) {
                case HTM_CONTAINS:
                    if (level == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (level < 20 && curcount >= tree->leafthresh) {
                        s = _htm_subdivide(curnode, s);
                        if (s == NULL) {
                            /* tree is invalid */
                            if (err != NULL) {
                                *err = HTM_EINV;
                            }
                            range.max = -1;
                            return range;
                        }
                        ++level;
                        ++curnode;
                        continue;
                    }
                    range.max += (int64_t) curcount;
                    break;
                case HTM_INSIDE:
                    /* fully covered HTM triangle */
                    range.min += (int64_t) curcount;
                    range.max += (int64_t) curcount;
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
ascend:
            --level;
            --curnode;
            while (level >= 0 && curnode->child == 4) {
                --curnode;
                --level;
            }
            if (level < 0) {
                /* finished with this root */
                break;
            }
            s = _htm_subdivide(curnode, curnode->s);
            if (s == NULL) {
                /* no non-empty children remain */
                goto ascend;
            }
            ++level;
            ++curnode;
        }
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return range;
}


struct htm_range htm_tree_s2ellipse_range(const struct htm_tree *tree,
                                          const struct htm_s2ellipse *ellipse,
                                          enum htm_errcode *err)
{
    struct _htm_path path;
    enum htm_root root;
    struct htm_range range;

    range.min = 0;
    range.max = 0;
    if (tree == NULL || ellipse== NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        range.max = -1;
        return range;
    }
    if (tree->indexfd == -1) {
        range.max = htm_tree_s2ellipse_scan(tree, ellipse, err);
        if (range.max >= 0) {
            range.min = range.max;
        }
        return range;
    }

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        const unsigned char *s = tree->root[root];
        int level = 0;

        if (s == NULL) {
            /* root contains no points */
            continue;
        }
        _htm_path_root(&path, root);

        while (1) {
            uint64_t curcount = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            s += 1 + htm_varint_nfollow(*s);
            switch (_htm_s2ellipse_htmcov(curnode, ellipse)) {
                case HTM_CONTAINS:
                    if (level == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (level < 20 && curcount >= tree->leafthresh) {
                        s = _htm_subdivide(curnode, s);
                        if (s == NULL) {
                            /* tree is invalid */
                            if (err != NULL) {
                                *err = HTM_EINV;
                            }
                            range.max = -1;
                            return range;
                        }
                        ++level;
                        ++curnode;
                        continue;
                    }
                    range.max += (int64_t) curcount;
                    break;
                case HTM_INSIDE:
                    /* fully covered HTM triangle */
                    range.min += (int64_t) curcount;
                    range.max += (int64_t) curcount;
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
ascend:
            --level;
            --curnode;
            while (level >= 0 && curnode->child == 4) {
                --curnode;
                --level;
            }
            if (level < 0) {
                /* finished with this root */
                break;
            }
            s = _htm_subdivide(curnode, curnode->s);
            if (s == NULL) {
                /* no non-empty children remain */
                goto ascend;
            }
            ++level;
            ++curnode;
        }
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return range;
}


struct htm_range htm_tree_s2cpoly_range(const struct htm_tree *tree,
                                        const struct htm_s2cpoly *poly,
                                        enum htm_errcode *err)
{
    double stackab[2*256 + 4];
    struct _htm_path path;
    enum htm_root root;
    double *ab;
    size_t nb;
    struct htm_range range;

    range.min = 0;
    range.max = 0;
    if (tree == NULL || poly == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        range.max = -1;
        return range;
    }
    if (tree->indexfd == -1) {
        range.max = htm_tree_s2cpoly_scan(tree, poly, err);
        if (range.max >= 0) {
            range.min = range.max;
        }
        return range;
    }
    nb = (2 * poly->n + 4) * sizeof(double);
    if (nb > sizeof(stackab)) {
        ab = (double *) malloc(nb);
        if (ab == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            range.max = -1;
            return range;
        }
    } else {
        ab = stackab;
    }

    for (root = HTM_S0; root <= HTM_N3; ++root) {
        struct _htm_node *curnode = path.node;
        const unsigned char *s = tree->root[root];
        int level = 0;

        if (s == NULL) {
            /* root contains no points */
            continue;
        }
        _htm_path_root(&path, root);

        while (1) {
            uint64_t curcount = htm_varint_decode(s);
            s += 1 + htm_varint_nfollow(*s);
            s += 1 + htm_varint_nfollow(*s);
            switch (_htm_s2cpoly_htmcov(curnode, poly, ab)) {
                case HTM_CONTAINS:
                    if (level == 0) {
                        /* no need to consider other roots */
                        root = HTM_N3;
                    } else {
                        /* no need to consider other children of parent */
                        curnode[-1].child = 4;
                    }
                    /* fall-through */
                case HTM_INTERSECT:
                    if (level < 20 && curcount >= tree->leafthresh) {
                        s = _htm_subdivide(curnode, s);
                        if (s == NULL) {
                            /* tree is invalid */
                            if (ab != stackab) {
                                free(ab);
                            }
                            if (err != NULL) {
                                *err = HTM_EINV;
                            }
                            range.max = -1;
                            return range;
                        }
                        ++level;
                        ++curnode;
                        continue;
                    }
                    range.max += (int64_t) curcount;
                    break;
                case HTM_INSIDE:
                    /* fully covered HTM triangle */
                    range.min += (int64_t) curcount;
                    range.max += (int64_t) curcount;
                    break;
                default:
                    /* HTM triangle does not intersect circle */
                    break;
            }
            /* ascend towards the root */
ascend:
            --level;
            --curnode;
            while (level >= 0 && curnode->child == 4) {
                --curnode;
                --level;
            }
            if (level < 0) {
                /* finished with this root */
                break;
            }
            s = _htm_subdivide(curnode, curnode->s);
            if (s == NULL) {
                /* no non-empty children remain */
                goto ascend;
            }
            ++level;
            ++curnode;
        }
    }
    if (ab != stackab) {
        free(ab);
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return range;
}


#ifdef __cplusplus
}
#endif
