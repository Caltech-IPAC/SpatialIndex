/** \file
    \brief      Minimalistic functions and types for spherical geometry.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#ifndef HTM_GEOMETRY_H
#define HTM_GEOMETRY_H

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ================================================================ */
/** \defgroup geom_point Spherical coordinates and 3-vectors
    @{
  */
/* ================================================================ */

/** Cartesian coordinates for a point in R3.
  */
struct htm_v3 {
    double x; /**< x coordinate value. */
    double y; /**< y coordinate value. */
    double z; /**< z coordinate value. */
};

/** Spherical coordinates (in degrees) for a point in S2.
  */
struct htm_sc {
    double lon; /**< Longitude angle (right ascension), degrees. */
    double lat; /**< Latitude angle (declination), degrees. */
};

/** Stores the given cartesian coordinates in \p out.

    \return
            - HTM_ENANINF    if \p x, \p y, or \p z is non-finite.
            - HTM_ENULLPTR   if \p out is NULL.
            - HTM_OK         on success.
  */
HTM_INLINE enum htm_errcode htm_v3_init(struct htm_v3 *out,
                                        double x,
                                        double y,
                                        double z)
{
    if (HTM_ISSPECIAL(x) || HTM_ISSPECIAL(y) || HTM_ISSPECIAL(z)) {
        return HTM_ENANINF;
    } else if (out == NULL) {
        return HTM_ENULLPTR;
    }
    out->x = x;
    out->y = y;
    out->z = z;
    return HTM_OK;
}

/** Stores the given spherical coordinates in \p out; coordinates are assumed
    to be in degrees.

    \return
            - HTM_ENANINF    if any input coordinate is non-finite.
            - HTM_ELAT       if \p lat is not in the <tt>[-90, 90]</tt>
                             degree range.
            - HTM_ENULLPTR   if \p out is NULL.
            - HTM_OK         on success.
  */
HTM_INLINE enum htm_errcode htm_sc_init(struct htm_sc *out,
                                        double lon_deg,
                                        double lat_deg)
{
    if (HTM_ISSPECIAL(lon_deg) || HTM_ISSPECIAL(lat_deg)) {
        return HTM_ENANINF;
    } else if (lat_deg < -90.0 || lat_deg > 90.0) {
        return HTM_ELAT;
    } else if (out == NULL) {
        return HTM_ENULLPTR;
    }
    out->lon = lon_deg;
    out->lat = lat_deg;
    return HTM_OK;
}

/** Stores the vector sum <tt>v1 + v2</tt> in \p out.
    Arguments must not be \c NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_add(struct htm_v3 *out,
                           const struct htm_v3 *v1,
                           const struct htm_v3 *v2)
{
    out->x = v1->x + v2->x;
    out->y = v1->y + v2->y;
    out->z = v1->z + v2->z;
}

/** Stores the vector difference <tt>v1 - v2</tt> in \p out.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_sub(struct htm_v3 *out,
                           const struct htm_v3 *v1,
                           const struct htm_v3 *v2)
{
    out->x = v1->x - v2->x;
    out->y = v1->y - v2->y;
    out->z = v1->z - v2->z;
}

/** Stores the vector <tt>v * -1</tt> in \p out.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_neg(struct htm_v3 *out, const struct htm_v3 *v)
{
    out->x = - v->x;
    out->y = - v->y;
    out->z = - v->z;
}

/** Stores the vector-scalar product <tt>v * s</tt> in \p out.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_mul(struct htm_v3 *out,
                           const struct htm_v3 *v,
                           double s)
{
    out->x = v->x * s;
    out->y = v->y * s;
    out->z = v->z * s;
}

/** Stores the vector-scalar quotient <tt>v / s</tt> in \p out.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_div(struct htm_v3 *out,
                           const struct htm_v3 *v,
                           double s)
{
    out->x = v->x / s;
    out->y = v->y / s;
    out->z = v->z / s;
}

/** Returns the dot product of the 3-vectors \p v1 and \p v2.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE double htm_v3_dot(const struct htm_v3 *v1,
                             const struct htm_v3 *v2)
{
    return v1->x * v2->x + v1->y * v2->y + v1->z * v2->z;
}

/** Stores the component-wise product of the 3-vectors \p v1 and \p v2 in
    \p out.  Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_cwise_mul(struct htm_v3 *out,
                                 const struct htm_v3 *v1,
                                 const struct htm_v3 *v2)
{
    out->x = v1->x * v2->x;
    out->y = v1->y * v2->y;
    out->z = v1->z * v2->z;
}

/** Returns the squared norm of the 3-vector \p v, which must not be NULL.
    Equivalent to htm_v3_dot(v, v).
  */
HTM_INLINE double htm_v3_norm2(const struct htm_v3 *v)
{
    return v->x * v->x + v->y * v->y + v->z * v->z;
}

/** Returns the norm of the 3-vector \p v, which must not be NULL.
  */
HTM_INLINE double htm_v3_norm(const struct htm_v3 *v)
{
    return sqrt(htm_v3_norm2(v));
}

/** Stores a normalized copy of \p v in \p out.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_normalize(struct htm_v3 *out, const struct htm_v3 *v)
{
    double norm = htm_v3_norm(v);
    out->x = v->x / norm;
    out->y = v->y / norm;
    out->z = v->z / norm;
}

/** Stores twice the vector cross product of \p v1 and \p v2 in \p out.

    Arguments must not be NULL pointers, but may alias.

    This is equivalent to the cross product of <tt>v2 + v1</tt> and
    <tt>v2 - v1</tt>, and is more robust than htm_v3_cross() for nearly
    identical input unit vectors. This is because when \p v1 and \p v2
    are nearly identical, the evaluation of <tt>v2 - v1</tt> is exact and,
    for unit \p v1 and \p v2, <tt>v1 + v2</tt> is perpendicular to
    <tt>v2 - v1</tt>.
  */
HTM_INLINE void htm_v3_rcross(struct htm_v3 *out,
                              const struct htm_v3 *v1,
                              const struct htm_v3 *v2)
{
    double x1 = v2->x + v1->x;
    double x2 = v2->x - v1->x;
    double y1 = v2->y + v1->y;
    double y2 = v2->y - v1->y;
    double z1 = v2->z + v1->z;
    double z2 = v2->z - v1->z;
    out->x = y1 * z2 - z1 * y2;
    out->y = z1 * x2 - x1 * z2;
    out->z = x1 * y2 - y1 * x2;
}

/** Stores the vector cross product of \p v1 and \p v2 in \p out.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE void htm_v3_cross(struct htm_v3 *out,
                             const struct htm_v3 *v1,
                             const struct htm_v3 *v2)
{
    double x = v1->y * v2->z - v1->z * v2->y;
    double y = v1->z * v2->x - v1->x * v2->z;
    double z = v1->x * v2->y - v1->y * v2->x;
    out->x = x;
    out->y = y;
    out->z = z;
}

/** Computes the N,E basis unit vectors at a position \p v on the unit sphere.

    \return
            - HTM_ENULLPTR      if any argument is NULL.
            - HTM_EZERONORM     if \p v has a norm of 0.
            - HTM_OK            on success.
  */
enum htm_errcode htm_v3_ne(struct htm_v3 *north,
                           struct htm_v3 *east,
                           const struct htm_v3 *v);

/** Computes the angle by which the normal of the plane defined by the origin,
    \p v1 and \p v2 should be rotated about the axis <tt>v1 - v2</tt>, such
    that the resulting plane is tangent to two small circles of radius \p r
    centered at \p v1 and \p v2.

    Note that \p v1 and \p v2 must be unit vectors!

    \return
            - HTM_ENULLPTR  if \p angle, \p v1 or \p v2 is NULL.
            - HTM_EDEGEN    if \p v1 is too close to \p v2.
            - HTM_EANG      if \p r is negative or too large.
            - HTM_OK        on success.
  */
enum htm_errcode htm_v3_tanrot(double *angle,
                               const struct htm_v3 *v1,
                               const struct htm_v3 *v2,
                               double r);

/** Rotates vector \p v around the axis \p k by
    the specified number of degrees.

    \return
            - HTM_ENULLPTR      if \p out, \p v or \p k is NULL.
            - HTM_EZERONORM     if \p k has a norm of 0.
            - HTM_OK            on success.
  */
enum htm_errcode htm_v3_rot(struct htm_v3 *out,
                            const struct htm_v3 *v,
                            const struct htm_v3 *k,
                            double angle_deg);

/** Computes the normalized sum of the given list of vectors.

    \return
            - HTM_ENULLPTR      if \p cen or \p points is NULL.
            - HTM_ELEN          if n == 0.
            - HTM_OK            on success.
  */
enum htm_errcode htm_v3_centroid(struct htm_v3 *cen,
                                 const struct htm_v3 *points,
                                 size_t n);

/** Converts the spherical coordinate pair \p p to a unit 3-vector and stores
    the results in \p out.

    \return
            - HTM_ENULLPTR  if \p out or \p p is NULL.
            - HTM_OK        on success.
  */
enum htm_errcode htm_sc_tov3(struct htm_v3 *out, const struct htm_sc *p);

/** Converts the 3-vector \p v to spherical coordinates and stores the results
    in \p out.  The vector \p v is not required to have unit norm.

    \return
            - HTM_ENULLPTR  if \p out or \p v is NULL.
            - HTM_OK        on success.
  */
enum htm_errcode htm_v3_tosc(struct htm_sc *out, const struct htm_v3 *v);


/* ================================================================ */
/** @}
    \defgroup geom_dist Angular separation and distance
    @{
  */
/* ================================================================ */

/** Returns the square of the distance between the unit vectors corresponding
    to points \p p1 and \p p2.  Arguments must not be NULL pointers, but may
    alias.
  */
double htm_sc_dist2(const struct htm_sc *p1, const struct htm_sc *p2);

/** Returns the angular separation (in degrees) between the points \p p1 and
    \p p2.  Arguments must not not be NULL pointers, but may alias.
  */
double htm_sc_angsep(const struct htm_sc *p1, const struct htm_sc *p2);

/** Returns the square of the distance betwen vectors \p v1 and \p v2.
    Arguments must not be NULL pointers, but may alias.
  */
HTM_INLINE double htm_v3_dist2(const struct htm_v3 *v1,
                               const struct htm_v3 *v2)
{
    struct htm_v3 delta;
    htm_v3_sub(&delta, v1, v2);
    return htm_v3_norm2(&delta);
}

/** Returns the angular separation (in degrees) between unit vectors
    \p v1 and \p v2.  Arguments must not be NULL pointers, but may alias.
 */
double htm_v3_angsepu(const struct htm_v3 *v1,
                      const struct htm_v3 *v2);

/** Returns the angular separation (in degrees) between vectors \p v1 and
    \p v2, which need not have unit norm.  Arguments must not be NULL
    pointers, but may alias.
 */
double htm_v3_angsep(const struct htm_v3 *v1, const struct htm_v3 *v2);

/** Returns the minimum square distance between \p v, and points on the edge
    from \p v1 to \p v2 (where \p e is a vector parallel to the cross product
    of \p v1 and \p v2). The vectors \p v, \p v1, and \p v2 are assumed to be
    normalized, but \p e need not have unit norm.
 */
double htm_v3_edgedist2(const struct htm_v3 *v,
                        const struct htm_v3 *v1,
                        const struct htm_v3 *v2,
                        const struct htm_v3 *e);


/* ================================================================ */
/** @}
    \defgroup geom_ellipse Spherical Ellipses
    @{
 */
/* ================================================================ */

/** An ellipse on the sphere, defined as the points p such that
    angsep(p, f1) + angsep(p, f2) <= 2a, where a is the semi-major axis
    angle of the ellipse and f1/f2 are its focii. The boundary is
    given by the intersection of an elliptical cone passing through the
    origin with the unit sphere; the cone is represented via the
    symmetric 3 by 3 matrix M of the cones quadratic form.
  */
struct htm_s2ellipse {
    struct htm_v3 cen; /**< Ellipse center (unit vector). */
    double xx;  /**< M[1,1] */
    double yy;  /**< M[2,2] */
    double zz;  /**< M[3,3] */
    double xy;  /**< M[1,2] (== M[2,1]) */
    double xz;  /**< M[1,3] (== M[3,1]) */
    double yz;  /**< M[2,3] (== M[3,2]) */
    double a;   /**< Semi-major axis angle (degrees). */
};

/** Initializes a spherical ellipse from the given focii (which must be
    unit vectors) and semi-major axis angle \p a (in degrees).
  */
enum htm_errcode htm_s2ellipse_init(struct htm_s2ellipse *e,
                                    const struct htm_v3 *f1,
                                    const struct htm_v3 *f2,
                                    double a);

/** Initializes a spherical ellipse from the given center (which must be
    a unit vector), axis angles \p a and \p b, and orientation (in degrees).
    The orientation is the position angle (east of north) of the first axis
    with respect to the north pole.
  */
enum htm_errcode htm_s2ellipse_init2(struct htm_s2ellipse *e,
                                     const struct htm_v3 *cen,
                                     double a,
                                     double b,
                                     double angle);

/** Returns 1 if the spherical ellipse \p e contains vector \p v,
    and 0 otherwise.  Arguments must not be NULL pointers.
  */
HTM_INLINE int htm_s2ellipse_cv3(const struct htm_s2ellipse *e,
                                 const struct htm_v3 *v)
{
    double qf = e->xx * v->x * v->x +
                e->yy * v->y * v->y +
                e->zz * v->z * v->z +
          2.0 * e->xy * v->x * v->y +
          2.0 * e->xz * v->x * v->z +
          2.0 * e->yz * v->y * v->z;
    double dp = htm_v3_dot(&e->cen, v);
    if (e->a <= 90.0) {
        return dp >= 0.0 && qf <= 0.0;
    }
    return dp >= 0.0 || qf >= 0.0;
}


/* ================================================================ */
/** @}
    \defgroup geom_poly Convex Spherical Polygons
    @{
 */
/* ================================================================ */

/** A convex polygon on the sphere.

    Note that vertices are stored as unit vectors, but that
    edge plane normals are not necessarily normalized.
  */
struct htm_s2cpoly {
    size_t n;               /**< number of vertices/edges. */
    struct htm_v3 vsum;     /**< sum of all vertices in polygon. */
    struct htm_v3 ve[];     /**< polygon vertices, followed by edges. */
};

/** Returns a pointer to the vertices of \p poly.
  */
HTM_INLINE struct htm_v3 * htm_s2cpoly_verts(struct htm_s2cpoly *poly)
{
    return poly->ve;
}

/** Returns a pointer to the edges of \p poly.
  */
HTM_INLINE struct htm_v3 * htm_s2cpoly_edges(struct htm_s2cpoly *poly)
{
    return poly->ve + poly->n;
}

/** Creates a polygon from a list of at least 3 vertices.
    Vertices can be in clockwise or counter-clockwise order, but:

    - vertices must be unit vectors
    - vertices must to be hemispherical
    - vertices must define edges that do not intersect except at vertices
    - vertices must define edges forming a convex polygon

    To release resources for the polygon, call free() on the returned pointer.

    \return
            - a newly allocated polygon on success.
            - NULL if an error occurs or the inputs are invalid. In this
              case, \p *err is additionally set to indicate the reason
              for the failure.
  */
struct htm_s2cpoly * htm_s2cpoly_init(const struct htm_v3 *verts,
                                      size_t n,
                                      enum htm_errcode *err);

/** Creates a polygon corresponding to the box with the given
    width, height, center and orientation. The orientation is specified
    as a rotation angle in degrees around the box center. An angle
    of zero degrees will result in a box aligned to the north (height)
    and east (width) axes at the box center.

    To release resources for the polygon, call free() on the returned pointer.

    \return
            - a newly allocated polygon on success.
            - NULL if an error occurs or the inputs are invalid. In this
              case, \p *err is additionally set to indicate the reason
              for the failure.
  */
struct htm_s2cpoly * htm_s2cpoly_box(const struct htm_v3 *cen,
                                     double width,
                                     double height,
                                     double angle,
                                     enum htm_errcode *err);

/** Creates an <tt>n</tt>-gon inscribed in the circle of the given
    radius and center.

    To release resources for the polygon, call free() on the returned pointer.

    \return
            - a newly allocated polygon on success.
            - NULL if an error occurs or the inputs are invalid. In this
              case, \p *err is additionally set to indicate the reason
              for the failure.
  */
struct htm_s2cpoly * htm_s2cpoly_ngon(const struct htm_v3 *cen,
                                      double r,
                                      size_t n,
                                      enum htm_errcode *err);

/** Creates a polygon corresponding to the quadrilateral completely enclosing
    two small circles with the given centers and a radius of r degrees.  The
    quadrilateral is oriented in the direction <tt>v2 - v1</tt>.

    To release resources for the polygon, call free() on the returned pointer.

    \return
            - a newly allocated polygon on success.
            - NULL if an error occurs or the inputs are invalid. In this
              case, \p *err is additionally set to indicate the reason
              for the failure.
  */
struct htm_s2cpoly * htm_s2cpoly_line(const struct htm_v3 *v1,
                                      const struct htm_v3 *v2,
                                      double r,
                                      enum htm_errcode *err);

/** Returns 1 if the spherical convex polygon \p poly contains
    vector \p v, and 0 otherwise.  Arguments must not be NULL pointers.
  */
int htm_s2cpoly_cv3(const struct htm_s2cpoly *poly, const struct htm_v3 *v);

/** Returns the area (in steradians) enclosed by the given spherical
    convex polygon. If \p poly is NULL, 0.0 is returned.
  */
double htm_s2cpoly_area(const struct htm_s2cpoly *poly);

/** Returns a copy of the given polygon, or \c NULL if \p poly is \c NULL
    or the required memory allocation failed.

    To release resources for the polygon, call free() on the returned pointer.
  */
struct htm_s2cpoly * htm_s2cpoly_clone(const struct htm_s2cpoly *poly);

/** Pads \p poly by the angle \p r (degrees).
  */
enum htm_errcode htm_s2cpoly_pad(struct htm_s2cpoly *poly, double r);


/* ================================================================ */
/** @}
    \defgroup geom_convex Convex Hull
    @{
 */
/* ================================================================ */

/** Returns 1 if the given set of points (which need not consist of unit
    vectors) is hemispherical and 0 otherwise. If an error occurs, 0 is
    returned and \p *err is set to indicate the reason for the failure.
  */
int htm_v3_hemispherical(const struct htm_v3 *points,
                         size_t n,
                         enum htm_errcode *err);

/** Tests whether an ordered list of points form a spherical convex polygon.

    \return
            -  1: \p points form a spherical convex polygon and
                  are in counter-clockwise order.
            -  0: \p points do not form a spherical convex polygon,
                  or an error occurred. In case of error, \p *err is
                  set to indicate the reason for the failure.
            - -1: \p points form a spherical convex polygon and are
                  in clockwise order.
  */
int htm_v3_convex(const struct htm_v3 *points,
                  size_t n,
                  enum htm_errcode *err);

/** Creates a polygon corresponding to the convex hull of the given point set.
    Points must be specified as unit vectors.

    To release resources for the polygon, call free() on the returned pointer.

    \return
            - a newly allocated polygon on success.
            - NULL if an error occurs or the inputs are invalid. In this
              case, \p *err is additionally set to indicate the reason
              for the failure.
  */
struct htm_s2cpoly * htm_s2cpoly_hull(const struct htm_v3 *points,
                                      size_t n,
                                      enum htm_errcode *err);

/** @} */


#ifdef __cplusplus
}
#endif

#endif /* HTM_GEOMETRY_H */

