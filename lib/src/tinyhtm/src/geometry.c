/** \file
    \brief      Spherical geometry implementation

    For API documentation, see geometry.h.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#include "tinyhtm/geometry.h"

#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "tinyhtm/select.h"

#ifdef __cplusplus
extern "C" {
#endif


/* ---- Spherical coordinates and 3-vectors ---- */

enum htm_errcode htm_v3_ne(struct htm_v3 *north,
                           struct htm_v3 *east,
                           const struct htm_v3 *v)
{
    if (north == NULL || east == NULL || v == NULL) {
        return HTM_ENULLPTR;
    } else if (htm_v3_norm2(v) == 0.0) {
        return HTM_EZERONORM;
    }
    north->x = - v->x * v->z;
    north->y = - v->y * v->z;
    north->z = v->x * v->x + v->y * v->y;
    if (north->x == 0.0 && north->y == 0.0 && north->z == 0.0) {
        /* pick an arbitrary orthogonal basis with z = 0 */
        north->x = -1.0;
        east->x = 0.0;
        east->y = 1.0;
        east->z = 0.0;
    } else {
        htm_v3_normalize(north, north);
        htm_v3_rcross(east, north, v);
        htm_v3_normalize(east, east);
    }
    return HTM_OK;
}


static const double HTM_NAN = 0.0 / 0.0;
static const double HTM_RMAX = 90.0 - 0.001/3600.0;

enum htm_errcode htm_v3_tanrot(double *angle,
                               const struct htm_v3 *v1,
                               const struct htm_v3 *v2,
                               double r)
{
    double a, s;
    if (angle == NULL || v1 == NULL || v2 == NULL) {
        return HTM_ENULLPTR;
    }
    if (r <= 0.0) {
        return HTM_EANG;
    }
    a = htm_v3_angsep(v1, v2);
    if (a == 0.0) {
        return HTM_EDEGEN;
    }
    if (a + 2.0*r > 2.0*HTM_RMAX) {
        return HTM_EANG;
    }
    r *= HTM_RAD_PER_DEG;
    a *= HTM_RAD_PER_DEG;
    s = 2.0 * sin(r) * sin(0.5*a) / sin(a);
    if (s >= 1.0) {
        *angle = 90.0;
    } else {
        *angle = asin(s) * HTM_DEG_PER_RAD;
    }
    return HTM_OK;
}


enum htm_errcode htm_v3_rot(struct htm_v3 *out,
                            const struct htm_v3 *v,
                            const struct htm_v3 *k,
                            double angle)
{
    struct htm_v3 kxv, tmp;
    double sina, cosa, nk, kdotv;
    if (out == NULL || v == NULL || k == NULL) {
        return HTM_ENULLPTR;
    }
    /* Rodrigues' rotation formula:
       v_rot = v cos(a) + (k ^ v) sin(a) + k (k . v) (1 - cos(a))
     */
    nk = htm_v3_norm(k);
    if (nk == 0.0) {
        return HTM_EZERONORM;
    }
    sina = sin(angle * HTM_RAD_PER_DEG);
    cosa = cos(angle * HTM_RAD_PER_DEG);
    kdotv = htm_v3_dot(k, v) / nk;
    htm_v3_rcross(&kxv, k, v);
    htm_v3_mul(&kxv, &kxv, 0.5*sina / nk);
    htm_v3_mul(&tmp, v, cosa);
    htm_v3_add(out, &kxv, &tmp);
    htm_v3_mul(&tmp, k, kdotv*(1.0 - cosa));
    htm_v3_add(out, out, &tmp);
    return HTM_OK;
}


enum htm_errcode htm_v3_centroid(struct htm_v3 *cen,
                                 const struct htm_v3 *points,
                                 size_t n)
{
    size_t i;
    if (cen == NULL || points == NULL) {
        return HTM_ENULLPTR;
    } else if (n == 0) {
        return HTM_ELEN;
    }
    cen->x = cen->y = cen->z = 0.0;
    for (i = 0; i < n; ++i) {
        cen->x += points[i].x;
        cen->y += points[i].y;
        cen->z += points[i].z;
    }
    htm_v3_normalize(cen, cen);
    return HTM_OK;
}


enum htm_errcode htm_sc_tov3(struct htm_v3 *out, const struct htm_sc *p)
{
    double lon, lat, cos_lat;
    if (out == NULL || p == NULL) {
        return HTM_ENULLPTR;
    }
    lon = p->lon * HTM_RAD_PER_DEG;
    lat = p->lat * HTM_RAD_PER_DEG;
    cos_lat = cos(lat);
    out->x = cos(lon) * cos_lat;
    out->y = sin(lon) * cos_lat;
    out->z = sin(lat);
    return HTM_OK;
}


enum htm_errcode htm_v3_tosc(struct htm_sc *out, const struct htm_v3 *v)
{
    double d2;
    if (out == NULL || v == NULL) {
        return HTM_ENULLPTR;
    }
    d2  = v->x*v->x + v->y*v->y;
    if (d2 == 0.0) {
        out->lon = 0.0;
    } else {
        double lon = atan2(v->y, v->x) * HTM_DEG_PER_RAD;
        if (lon < 0.0) {
            lon += 360.0;
            if (lon == 360.0) {
                lon = 0.0;
            }
        }
        out->lon = lon;
    }
    if (v->z == 0.0) {
        out->lat = 0.0;
    } else {
        out->lat = htm_clamp(atan2(v->z, sqrt(d2)) * HTM_DEG_PER_RAD,
                             -90.0, 90.0);
    }
    return HTM_OK;
}


/* ---- Angular separation and distance ---- */

double htm_sc_dist2(const struct htm_sc *p1, const struct htm_sc *p2)
{
    double x, y, z, d2;
    x = sin((p1->lon - p2->lon) * HTM_RAD_PER_DEG * 0.5);
    x *= x;
    y = sin((p1->lat - p2->lat) * HTM_RAD_PER_DEG * 0.5);
    y *= y;
    z = cos((p1->lat + p2->lat) * HTM_RAD_PER_DEG * 0.5);
    z *= z;
    d2 = 4.0 * (x * (z - y) + y);
    return d2 < 0.0 ? 0.0 : (d2 > 4.0 ? 4.0 : d2);
}


double htm_sc_angsep(const struct htm_sc *p1, const struct htm_sc *p2)
{
    double angsep, x;
    x = htm_sc_dist2(p1, p2) * 0.25;
    angsep = 2.0 * HTM_DEG_PER_RAD * asin(sqrt(x));
    return angsep > 180.0 ? 180.0 : angsep;
}


double htm_v3_angsepu(const struct htm_v3 *unit_v1,
                      const struct htm_v3 *unit_v2)
{
    double angsep, x;
    x = htm_v3_dist2(unit_v1, unit_v2) * 0.25;
    angsep = 2.0 * HTM_DEG_PER_RAD * asin(sqrt(x > 1.0 ? 1.0 : x));
    return angsep > 180.0 ? 180.0 : angsep;
}


double htm_v3_angsep(const struct htm_v3 *v1, const struct htm_v3 *v2)
{
    struct htm_v3 n;
    double ss, cs, angsep;
    htm_v3_cross(&n, v1, v2);
    ss = htm_v3_norm(&n);
    cs = htm_v3_dot(v1, v2);
    if (cs == 0.0 && ss == 0.0) {
        return 0.0;
    }
    angsep = atan2(ss, cs) * HTM_DEG_PER_RAD;
    return (angsep > 180.0) ? 180.0 : angsep;
}


double htm_v3_edgedist2(const struct htm_v3 *v,
                        const struct htm_v3 *v1,
                        const struct htm_v3 *v2,
                        const struct htm_v3 *e)
{
    struct htm_v3 c;
    htm_v3_cross(&c, v, e);
    if (htm_v3_dot(&c, v1) > 0.0 && htm_v3_dot(&c, v2) < 0.0) {
        double d = htm_v3_dot(v, e);
        double x = d * d / htm_v3_norm2(e);
        double y;
        /* x is the square of the sin of the minimum angle between v and the
           edge. To map to a square secant distance, compute
           2.0*(1 - sqrt(1 - x)) */
        if (x > 1.0) {
            return 2.0;
        } else if (x < 1.0e-7) {
            /* for small x, use taylor series to compute a result accurate to
               about 1 ulp */
            y = x * x;
            return x + (0.25*y + 0.125*x*y);
        }
        y = 1.0 - sqrt(1.0 - x);
        /* 1 newton-raphson iteration to improve accuracy. */
        return (x - y*y)/(1 - y);
    } else {
        double d1, d2;
        d1 = htm_v3_dist2(v, v1);
        d2 = htm_v3_dist2(v, v2);
        return d1 < d2 ? d1 : d2;
    }
}


/* ---- Spherical Ellipses ---- */

enum htm_errcode htm_s2ellipse_init(struct htm_s2ellipse * const ell,
                                    const struct htm_v3 * const f1,
                                    const struct htm_v3 * const f2,
                                    const double a)
{
    double e, ss, c;
    struct htm_v3 f11, f22, f12;

    e = 0.5 * htm_v3_angsepu(f1, f2);
    if (e > 90.0 - 2.777777777777777778e-6 || a <= e || a >= 180.0 - e) {
        return HTM_EANG;
    }
    htm_v3_add(&ell->cen, f1, f2);
    htm_v3_normalize(&ell->cen, &ell->cen);
    ss = sin(2.0 * HTM_RAD_PER_DEG * a);
    c = cos(2.0 * HTM_RAD_PER_DEG * a);
    ss *= ss;
    htm_v3_cwise_mul(&f11, f1, f1);
    htm_v3_cwise_mul(&f22, f2, f2);
    htm_v3_cwise_mul(&f12, f1, f2);
    ell->xx = ss - f11.x - f22.x + 2.0*c*f12.x;
    ell->yy = ss - f11.y - f22.y + 2.0*c*f12.y;
    ell->zz = ss - f11.z - f22.z + 2.0*c*f12.z;
    ell->xy = c*(f1->x*f2->y + f1->y*f2->x) - f1->x*f1->y - f2->x*f2->y;
    ell->xz = c*(f1->x*f2->z + f1->z*f2->x) - f1->x*f1->z - f2->x*f2->z;
    ell->yz = c*(f1->y*f2->z + f1->z*f2->y) - f1->y*f1->z - f2->y*f2->z;
    ell->a = a;
    return HTM_OK;
}


enum htm_errcode htm_s2ellipse_init2(struct htm_s2ellipse *ellipse,
                                     const struct htm_v3 *cen,
                                     double a,
                                     double b,
                                     double angle)
{
    struct htm_v3 N, E, n, e;
    double s, c;

    if (ellipse == NULL || cen == NULL) {
        return HTM_ENULLPTR;
    }
    if (a <= 0.0 || b <= 0.0 ||
        a > 90.0 - 2.777777777777777778e-6 ||
        b > 90.0 - 2.777777777777777778e-6) {
        return HTM_EANG;
    }
    /* elliptical cone is given by:

       x^2/tan(a)^2 + y^2/tan(b)^2 - z^2 = 0
       
       in the basis [n, e, cen], where n and e are the north/east vectors
       at cen rotated clockwise by posang.

       Let M be the 3x3 orthogonal matrix with rows n, e and cen, and let M'
       be its tranpose. Then the matrix Q of the elliptical cone is:

       Q = M' D M

       where D is the diagonal matrix:

       D = [ 1/tan(a)^2  0            0
             0           1/tan(b)^2   0
             0           0           -1 ]
     */
    ellipse->cen = *cen;
    ellipse->a = a;
    a = tan(HTM_RAD_PER_DEG * a);
    b = tan(HTM_RAD_PER_DEG * b);
    a = 1.0 / (a*a);
    b = 1.0 / (b*b);
    htm_v3_ne(&N, &E, cen);
    /* N, E is the north, east basis at cen */
    s = sin(HTM_RAD_PER_DEG * angle);
    c = cos(HTM_RAD_PER_DEG * angle);
    htm_v3_mul(&n, &N, c);
    htm_v3_mul(&e, &E, s);
    htm_v3_sub(&n, &n, &e); /* n = cos(angle)*N - sin(angle)*E */
    htm_v3_mul(&N, &N, s);
    htm_v3_mul(&E, &E, c);
    htm_v3_add(&e, &N, &E); /* e = sin(angle)*N + cos(angle)*E */
    /* have the elements of M and D, compute and store Q */
    ellipse->xx = a * n.x*n.x  +  b * e.x*e.x  -  cen->x*cen->x;
    ellipse->yy = a * n.y*n.y  +  b * e.y*e.y  -  cen->y*cen->y;
    ellipse->zz = a * n.z*n.z  +  b * e.z*e.z  -  cen->z*cen->z;
    ellipse->xy = a * n.x*n.y  +  b * e.x*e.y  -  cen->x*cen->y;
    ellipse->xz = a * n.x*n.z  +  b * e.x*e.z  -  cen->x*cen->z;
    ellipse->yz = a * n.y*n.z  +  b * e.y*e.z  -  cen->y*cen->z;
    return HTM_OK;
}


/* ---- Convex Spherical Polygons ---- */

HTM_INLINE int _htm_nv_valid(size_t n)
{
    return (n >= 3 &&
            n < (SIZE_MAX - sizeof(struct htm_s2cpoly)) /
                (2 * sizeof(struct htm_v3)));
}

struct htm_s2cpoly * htm_s2cpoly_init(const struct htm_v3 *verts,
                                      size_t n,
                                      enum htm_errcode *err)
{
    struct htm_s2cpoly *out;
    size_t i;
    if (verts == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return NULL;
    }
    if (!_htm_nv_valid(n)) {
        if (err != NULL) {
            *err = HTM_ELEN;
        }
        return NULL;
    }
    out = (struct htm_s2cpoly *) malloc(
        sizeof(struct htm_s2cpoly) + 2 * n * sizeof(struct htm_v3));
    if (out == NULL) {
        if (err != NULL) {
            *err = HTM_ENOMEM;
        }
        return NULL;
    }
    out->n = n;
    out->vsum = verts[n - 1];
    for (i = 0; i < n - 1; ++i) {
        /* the cross product of two consecutive vertices gives a vector
           parallel to the edge plane normal. */
        htm_v3_rcross(&out->ve[i + n], &verts[i], &verts[i + 1]);
        htm_v3_add(&out->vsum, &out->vsum, &verts[i]);
    }
    /* compute last edge plane */
    htm_v3_rcross(&out->ve[2*n - 1], &verts[n - 1], &verts[0]);
    /* if vertices are clockwise, then the dot-product of vsum with
       any edge plane is negative. */
    if (htm_v3_dot(&out->vsum, &out->ve[n]) < 0.0) {
        struct htm_v3 tmp;
        /* reorder and invert edge plane normals */
        for (i = 0; i < n/2; ++i) {
            tmp = out->ve[i + n];
            htm_v3_neg(&out->ve[i + n], &out->ve[2*n - i - 2]);
            htm_v3_neg(&out->ve[2*n - i - 2], &tmp);
        }
        htm_v3_neg(&out->ve[2*n - 1], &out->ve[2*n - 1]);
        for (i = 0; i < n; ++i) {
            out->ve[i] = verts[n - i - 1];
        }
    } else {
        memcpy(out->ve, verts, n * sizeof(struct htm_v3));
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return out;
}


struct htm_s2cpoly * htm_s2cpoly_box(const struct htm_v3 *cen,
                                     double width,
                                     double height,
                                     double angle,
                                     enum htm_errcode *err)
{
    struct htm_v3 verts[4];
    struct htm_v3 edges[4];
    struct htm_v3 north, east;

    if (cen == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return NULL;
    }
    if (width <= 0.0 || height <= 0.0 ||
        width >= HTM_RMAX || height >= HTM_RMAX) {
        if (err != NULL) {
            *err = HTM_EANG;
        }
        return NULL;
    }
    if (htm_v3_norm2(cen) == 0.0) {
        if (err != NULL) {
            *err = HTM_EZERONORM;
        }
        return NULL;
    }
    /* compute N,E at cen */
    htm_v3_ne(&north, &east, cen);
    /* rotate N by +- height/2 around E, and E by +- width/2 around N
       to obtain edge plane normals */
    htm_v3_rot(&edges[0], &east, &north, 0.5*width);
    htm_v3_rot(&edges[2], &east, &north, -0.5*width);
    htm_v3_rot(&edges[1], &north, &east, -0.5*height);
    htm_v3_rot(&edges[3], &north, &east, 0.5*height);
    /* cross products of edge plane normals yields vertices */
    htm_v3_rcross(&verts[0], &edges[0], &edges[1]);
    htm_v3_normalize(&verts[0], &verts[0]);
    htm_v3_rcross(&verts[1], &edges[2], &edges[1]);
    htm_v3_normalize(&verts[1], &verts[1]);
    htm_v3_rcross(&verts[2], &edges[2], &edges[3]);
    htm_v3_normalize(&verts[2], &verts[2]);
    htm_v3_rcross(&verts[3], &edges[0], &edges[3]);
    htm_v3_normalize(&verts[3], &verts[3]);
    /* rotate to desired orientation ... */
    if (angle != 0.0) {
        htm_v3_rot(&verts[0], &verts[0], cen, angle);
        htm_v3_rot(&verts[1], &verts[1], cen, angle);
        htm_v3_rot(&verts[2], &verts[2], cen, angle);
        htm_v3_rot(&verts[3], &verts[3], cen, angle);
    }
    /* ... and build polygon from vertices */
    return htm_s2cpoly_init(verts, 4, err);
}


/*  Utility function to build an N-gon inscribed in the given circle.
 */
struct htm_s2cpoly * htm_s2cpoly_ngon(const struct htm_v3 *cen,
                                      double r,
                                      size_t n,
                                      enum htm_errcode *err)
{
    struct htm_s2cpoly *out;
    struct htm_v3 *verts;
    struct htm_v3 north, east, v;
    double sr, cr;
    size_t i;

    if (cen == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return NULL;
    }
    if (r <= 0.0 || r >= HTM_RMAX) {
        if (err != NULL) {
            *err = HTM_EANG;
        }
        return NULL;
    }
    if (!_htm_nv_valid(n)) {
        if (err != NULL) {
            *err = HTM_ELEN;
        }
        return NULL;
    }
    if (htm_v3_norm2(cen) == 0.0) {
        if (err != NULL) {
            *err = HTM_EZERONORM;
        }
        return NULL;
    }
    verts = (struct htm_v3 *) malloc(n * sizeof(struct htm_v3));
    if (verts == NULL) {
        if (err != NULL) {
            *err = HTM_ENOMEM;
        }
        return NULL;
    }
    htm_v3_ne(&north, &east, cen);
    sr = sin(r * HTM_RAD_PER_DEG);
    cr = cos(r * HTM_RAD_PER_DEG);
    for (i = 0; i < n; ++i) {
        double ang, sa, ca;
        ang = (HTM_RAD_PER_DEG * 360.0 * i) / n;
        sa = sin(ang);
        ca = cos(ang);
        v.x = ca * north.x + sa * east.x;
        v.y = ca * north.y + sa * east.y;
        v.z = ca * north.z + sa * east.z;
        verts[i].x = cr * cen->x + sr * v.x;
        verts[i].y = cr * cen->y + sr * v.y;
        verts[i].z = cr * cen->z + sr * v.z;
        htm_v3_normalize(&verts[i], &verts[i]);
    }
    out = htm_s2cpoly_init(verts, n, err);
    free(verts);
    return out;
}


struct htm_s2cpoly * htm_s2cpoly_line(const struct htm_v3 *v1,
                                      const struct htm_v3 *v2,
                                      double r,
                                      enum htm_errcode *err)
{
    struct htm_v3 verts[4];
    struct htm_v3 edges[4];
    struct htm_v3 axis1, axis2;
    double a;
    enum htm_errcode e = htm_v3_tanrot(&a, v1, v2, r);
    if (e != HTM_OK) {
        if (err != NULL) {
            *err = e;
        }
        return NULL;
    }
    /* compute rotation axes */
    htm_v3_sub(&axis1, v1, v2);
    htm_v3_rcross(&axis2, v1, v2);
    /* compute edge plane normals */
    htm_v3_rot(&edges[0], &axis2, &axis1, a);
    htm_v3_rcross(&edges[1], v1, &axis2);
    htm_v3_rot(&edges[1], &edges[1], &axis2, -r);
    htm_v3_rot(&edges[2], &axis2, &axis1, -a);
    htm_v3_rcross(&edges[3], v2, &axis2);
    htm_v3_rot(&edges[3], &edges[3], &axis2, r);
    /* cross products of edge plane normals yields vertices */
    htm_v3_rcross(&verts[0], &edges[0], &edges[1]);
    htm_v3_normalize(&verts[0], &verts[0]);
    htm_v3_rcross(&verts[1], &edges[2], &edges[1]);
    htm_v3_normalize(&verts[1], &verts[1]);
    htm_v3_rcross(&verts[2], &edges[2], &edges[3]);
    htm_v3_normalize(&verts[2], &verts[2]);
    htm_v3_rcross(&verts[3], &edges[0], &edges[3]);
    htm_v3_normalize(&verts[3], &verts[3]);
    /* ... and build polygon from vertices */
    return htm_s2cpoly_init(verts, 4, err);
}


int htm_s2cpoly_cv3(const struct htm_s2cpoly *cp, const struct htm_v3 *v)
{
    size_t i;
    const size_t n = cp->n;
    for (i = 0; i < n; ++i) {
        if (htm_v3_dot(v, &cp->ve[n + i]) < 0.0) {
            return 0;
        }
    }
    return 1;
}


double htm_s2cpoly_area(const struct htm_s2cpoly *poly)
{
    double asum;
    size_t i, n;

    if (poly == NULL) {
        return 0.0;
    }
    /* see Girard's theorem */
    for (i = 0, n = poly->n, asum = 0.0; i < n; ++i) {
        struct htm_v3 v;
        double sina, cosa;
        size_t j = (i == 0) ? n - 1 : i - 1;
        htm_v3_rcross(&v, &poly->ve[n + j], &poly->ve[n + i]);
        sina = 0.5 * htm_v3_norm(&v);
        cosa = - htm_v3_dot(&poly->ve[n + j], &poly->ve[n + i]);
        asum += atan2(sina, cosa);
    }
    return asum - (poly->n - 2) * M_PI;
}


struct htm_s2cpoly * htm_s2cpoly_clone(const struct htm_s2cpoly *poly)
{
    struct htm_s2cpoly *clone;
    size_t n;

    if (poly == NULL) {
        return NULL;
    }
    n = sizeof(struct htm_s2cpoly) + 2 * poly->n * sizeof(struct htm_v3);
    clone = (struct htm_s2cpoly *) malloc(n);
    if (clone != NULL) {
        memcpy(clone, poly, n);
    }
    return clone;
}


enum htm_errcode htm_s2cpoly_pad(struct htm_s2cpoly *poly, double r)
{
    struct htm_v3 stackbuf[128];
    struct htm_v3 *vecs;
    struct htm_v3 tmp;
    double angle;
    size_t i, n;
    int hemis;
    enum htm_errcode err = HTM_OK;

    if (poly == NULL) {
        return HTM_ENULLPTR;
    } else if (r < 0.0) {
        return HTM_EANG;
    } else if (r == 0.0) {
        return HTM_OK;
    }
    n = poly->n;
    if (n > sizeof(stackbuf) / (2*sizeof(struct htm_v3))) {
        /* allocate scratch space */
        vecs = (struct htm_v3 *) malloc(n * 2 * sizeof(struct htm_v3));
        if (vecs == NULL) {
            return HTM_ENOMEM;
        }
    } else {
        vecs = stackbuf;
    }

    /* rotate edge plane normals outward */
    for (i = 0; i < n; ++i) {
        size_t j = (i == 0 ? n - 1: i - 1);
        err = htm_v3_tanrot(&angle, &poly->ve[j], &poly->ve[i], r);
        if (err != HTM_OK) {
            goto cleanup;
        }
        htm_v3_sub(&tmp, &poly->ve[i], &poly->ve[j]);
        htm_v3_rot(&vecs[j], &poly->ve[n + j], &tmp, angle);
    }

    /* compute new vertices */
    for (i = 0; i < n; ++i) {
        size_t j = (i == 0 ? n - 1: i - 1);
        htm_v3_rcross(&vecs[n + i], &vecs[j], &vecs[i]);
        htm_v3_normalize(&vecs[n + i], &vecs[n + i]);
        if (htm_v3_dot(&vecs[n + i], &poly->ve[i]) < 0.0) {
            htm_v3_neg(&vecs[n + i], &vecs[n + i]);
        }
    }

    /* checks that the union of old and new vertices is hemispherical */
    memcpy(vecs, poly->ve, n * sizeof(struct htm_v3));
    hemis = htm_v3_hemispherical(vecs, 2*n, &err);
    if (err != HTM_OK) {
        goto cleanup;
    } else if (!hemis) {
        err = HTM_EANG;
        goto cleanup;
    }
    /* copy in new vertices and compute edge plane normals */
    memcpy(poly->ve, vecs + n, n * sizeof(struct htm_v3));
    poly->vsum = vecs[2*n - 1];
    for (i = 0; i < n - 1; ++i) {
        htm_v3_rcross(&poly->ve[i + n], &poly->ve[i], &poly->ve[i + 1]);
        htm_v3_add(&poly->vsum, &poly->vsum, &poly->ve[i]);
    }
    /* compute last edge plane */
    htm_v3_rcross(&poly->ve[2*n - 1], &poly->ve[n - 1], &poly->ve[0]);

cleanup:
    if (vecs != stackbuf) {
        free(vecs);
    }
    return err;
}


/* ---- Testing whether a set of points is hemispherical ---- */

/*  This test is used by the convexity test and convex hull algorithm. It is
    implemented using Megiddo's algorithm for linear programming in R2, see:

    Megiddo, N. 1982. Linear-time algorithms for linear programming in R3 and related problems.
    In Proceedings of the 23rd Annual Symposium on Foundations of Computer Science (November 03 - 05, 1982).
    SFCS. IEEE Computer Society, Washington, DC, 329-338. DOI= http://dx.doi.org/10.1109/SFCS.1982.74
 */

/** A pair of doubles.
  */
struct _htm_pair {
    double first;
    double second;
};


/** A list of pairs.
  */
struct _htm_pairs {
    size_t n;                  /**< Number of pairs in list. */
    struct _htm_pair pairs[];  /**< Pair array. */
};


HTM_INLINE void _htm_pairs_append(struct _htm_pairs *pairs,
                                  const struct _htm_pair *p)
{
    size_t n = pairs->n;
    pairs->pairs[n] = *p;
    pairs->n = n + 1;
}


static const double HTM_INF = 1.0 / 0.0;


static void _htm_g(double *out,
                   const struct _htm_pairs *constraints,
                   double x)
{
    size_t i;
    double ai = constraints->pairs[0].first;
    double v = ai * x + constraints->pairs[0].second;
    double amin = ai;
    double amax = ai;
    for (i = 1; i < constraints->n; ++i) {
        double ai = constraints->pairs[i].first;
        double vi = ai * x + constraints->pairs[i].second;
        if (vi == v) {
            if (ai < amin) {
                amin = ai;
            }
            if (ai > amax) {
                amax = ai;
            }
        } else if (vi > v) {
            v = vi;
            amin = ai;
            amax = ai;
        }
    }
    out[0] = v;
    out[1] = amin;
    out[2] = amax;
}


static void _htm_h(double *out,
                   const struct _htm_pairs *constraints,
                   double x)
{
    size_t i;
    double ai = constraints->pairs[0].first;
    double v = ai * x + constraints->pairs[0].second;
    double amin = ai;
    double amax = ai;
    for (i = 1; i < constraints->n; ++i) {
        double ai = constraints->pairs[i].first;
        double vi = ai * x + constraints->pairs[i].second;
        if (vi == v) {
            if (ai < amin) {
                amin = ai;
            }
            if (ai > amax) {
                amax = ai;
            }
        } else if (vi < v) {
            v = vi;
            amin = ai;
            amax = ai;
        }
    }
    out[0] = v;
    out[1] = amin;
    out[2] = amax;
}


static size_t _htm_prune_g(double *intersections,
                           size_t ni,
                           struct _htm_pairs *constraints,
                           const struct _htm_pair *x)
{
    size_t i = 0;
    size_t n = constraints->n - 1;
    while (i < n) {
        double xx;
        double a1 = constraints->pairs[i].first;
        double b1 = constraints->pairs[i].second;
        double a2 = constraints->pairs[i + 1].first;
        double b2 = constraints->pairs[i + 1].second;
        double da = a1 - a2;
        if (fabs(da) < DBL_MIN / DBL_EPSILON) {
            xx = HTM_INF;
        } else {
            xx = (b2 - b1) / da;
        }
        if (HTM_ISSPECIAL(xx)) {
            if (b1 > b2) {
                constraints->pairs[i + 1] = constraints->pairs[n];
            } else {
                constraints->pairs[i] = constraints->pairs[n];
            }
            --n;
        } else {
            if (xx <= x->first) {
                if (a1 > a2) {
                    constraints->pairs[i + 1] = constraints->pairs[n];
                } else {
                    constraints->pairs[i] = constraints->pairs[n];
                }
                --n;
            } else if (xx >= x->second) {
                if (a1 > a2) {
                    constraints->pairs[i] = constraints->pairs[n];
                } else {
                    constraints->pairs[i + 1] = constraints->pairs[n];
                }
                --n;
            } else {
                /* save intersection */
                intersections[ni] = xx;
                ++ni;
                i += 2;
            }
        }
    }
    constraints->n = n + 1;
    return ni;
}


static size_t _htm_prune_h(double *intersections,
                           size_t ni,
                           struct _htm_pairs *constraints,
                           const struct _htm_pair *x)
{
    size_t i = 0;
    size_t n = constraints->n - 1;
    while (i < n) {
        double xx;
        double a1 = constraints->pairs[i].first;
        double b1 = constraints->pairs[i].second;
        double a2 = constraints->pairs[i + 1].first;
        double b2 = constraints->pairs[i + 1].second;
        double da = a1 - a2;
        if (fabs(da) < DBL_MIN / DBL_EPSILON) {
            xx = HTM_INF;
        } else {
            xx = (b2 - b1) / da;
        }
        if (HTM_ISSPECIAL(xx)) {
            if (b1 < b2) {
                constraints->pairs[i + 1] = constraints->pairs[n];
            } else {
                constraints->pairs[i] = constraints->pairs[n];
            }
            --n;
        } else {
            if (xx <= x->first) {
                if (a1 < a2) {
                    constraints->pairs[i + 1] = constraints->pairs[n];
                } else {
                    constraints->pairs[i] = constraints->pairs[n];
                }
                --n;
            } else if (xx >= x->second) {
                if (a1 < a2) {
                    constraints->pairs[i] = constraints->pairs[n];
                } else {
                    constraints->pairs[i + 1] = constraints->pairs[n];
                }
                --n;
            } else {
                /* save intersection */
                intersections[ni] = xx;
                ++ni;
                i += 2;
            }
        }
    }
    constraints->n = n + 1;
    return ni;
}


static int _htm_feasible_2d(struct _htm_pairs *I1,
                            struct _htm_pairs *I2,
                            double *intersections,
                            const struct htm_v3 *points,
                            size_t n,
                            double z)
{
    struct _htm_pair xr;
    size_t i, ni;
    double g[3];
    double h[3];

    xr.first = -HTM_INF;
    xr.second = HTM_INF;
    I1->n = 0;
    I2->n = 0;
    /* transform each constraint of the form x*v.x + y*v.y + z*v.z > 0
       into y op a*x + b or x op c, where op is < or > */
    for (i = 0; i < n; ++i) {
        if (fabs(points[i].y) <= DBL_MIN) {
            if (fabs(points[i].x) <= DBL_MIN) {
                if (z * points[i].z <= 0.0) {
                    /* inequalities trivially lack a solution */
                    return 0;
                }
                /* current inequality is trivially true, skip it */
            } else {
                double xlim = - z * points[i].z / points[i].x;
                if (points[i].x > 0.0) {
                    if (xlim > xr.first) {
                        xr.first = xlim;
                    }
                } else {
                    if (xlim < xr.second) {
                        xr.second = xlim;
                    }
                }
                if (xr.second <= xr.first) {
                    /* inequalities trivially lack a solution */
                    return 0;
                }
            }
        } else {
            /* finite since |z|, |points[i].xyz| <= 1.0 and 1/DBL_MIN
               are finite */
            struct _htm_pair coeffs;
            coeffs.first = - points[i].x / points[i].y;
            coeffs.second = - z * points[i].z / points[i].y;
            if (points[i].y > 0.0) {
                _htm_pairs_append(I1, &coeffs);
            } else {
                _htm_pairs_append(I2, &coeffs);
            }
        }
    }
    /* At this point (xmin, xmax) is non-empty - if either I1 or I2 is empty
       then a solution trivially exists. */
    if (I1->n == 0 || I2->n == 0) {
        return 1;
    }

    /* Check for a feasible solution to the inequalities I1 of the form
       y > a*x + b, I2 of the form y < a*x + b, x > xr.min and x < xr.max */
    while (1) {
        double med;
        ni = _htm_prune_g(intersections, 0, I1, &xr);
        ni = _htm_prune_h(intersections, ni, I2, &xr);
        if (ni == 0) {
            /* I1 and I2 each contain exactly one constraint */
            double a1 = I1->pairs[0].first;
            double b1 = I1->pairs[0].second;
            double a2 = I2->pairs[0].first;
            double b2 = I2->pairs[0].second;
            double xi;
            xi = (b2 - b1) / (a1 - a2);
            if (HTM_ISSPECIAL(xi)) {
                return b1 < b2;
            }
            return (xi > xr.first || a1 < a2) && (xi < xr.second || a1 > a2);
        }
        med = htm_select(intersections, ni, ni >> 1);
        /* If g(x) < h(x), x is a feasible solution. Otherwise, refine the
           search interval by examining the one-sided derivates of g/h. */
        _htm_g(g, I1, med);
        _htm_h(h, I2, med);
        if (g[0] <= h[0]) {
            return 1;
        } else if (g[1] > h[2]) {
            xr.second = med;
        } else if (g[2] < h[1]) {
            xr.first = med;
        } else {
            return 0;
        }
    }
    /* never reached */
    return 0;
}


static int _htm_feasible_1d(const struct htm_v3 *points, size_t n, double y)
{
    double xmin = - HTM_INF;
    double xmax = HTM_INF;
    size_t i;
    for (i = 0; i < n; ++i) {
        double xp = points[i].x;
        double yp = points[i].y;
        if (fabs(xp) <= DBL_MIN) {
            if (y * yp <= 0.0) {
                return 0;
            }
            /* inequality is trivially true, skip it */
        } else {
            double xlim = - y * yp / xp;
            if (xp > 0.0) {
                if (xlim > xmin) {
                    xmin = xlim;
                }
            } else if (xlim < xmax) {
                xmax = xlim;
            }
            if (xmax <= xmin) {
                return 0;
            }
        }
    }
    return 1;
}


int htm_v3_hemispherical(const struct htm_v3 *points,
                         size_t n,
                         enum htm_errcode *err)
{
    unsigned char stackbuf[4096];
    unsigned char *buf;
    struct _htm_pairs *I1;
    struct _htm_pairs *I2;
    double *intersections;
    size_t i, b;
    int pos, neg;

    if (points == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return 0;
    }
    if (n == 0) {
        if (err != NULL) {
            *err = HTM_ELEN;
        }
        return 0;
    }
    /* setup pointers to 2 _htm_pairs structures, each containing at most
       n pairs, and an array of n doubles. Use stack for small inputs,
       perform a single allocation otherwise. */
    i = sizeof(struct _htm_pairs) + n*sizeof(struct _htm_pair);
    b = 2*i + 3*sizeof(struct _htm_pairs) + (n + 2) * sizeof(double);
    if (b <= sizeof(stackbuf)) {
        buf = stackbuf;
    } else {
        buf = (unsigned char *) malloc(b);
        if (buf == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            return 0;
        }
    }
    /* compute b s.t. buf + b is a multiple of sizeof(struct _htm_pairs) */
    b = sizeof(struct _htm_pairs) - ((size_t) buf) % sizeof(struct _htm_pairs);
    /* compute i, a multiple of sizeof(struct _htm_pairs) bytes s.t. an
       _htm_pairs followed by n _htm_pair structures fits in i */
    i += sizeof(struct _htm_pairs) - i % sizeof(struct _htm_pairs);
    /* setup scratch space pointers */
    I1 = (struct _htm_pairs *) (buf + b);
    I2 = (struct _htm_pairs *) (buf + (i + b));
    i = 2 * i + b;
    i += sizeof(double) - ((size_t) (buf + i)) % sizeof(double);
    intersections = (double *) (buf + i);

    /* Check whether the set of linear equations
       x*v[0] + y*v[1] + z*v[2] > 0.0 (for v in points)
       has a solution (x, y, z). If (x,y,z) is a solution (is feasible),
       so is C*(x,y,z), C > 0. Therefore we can fix z to 1, -1 and
       perform 2D feasibility tests. */
    if (_htm_feasible_2d(I1, I2, intersections, points, n, 1.0) != 0) {
        goto feasible;
    }
    if (_htm_feasible_2d(I1, I2, intersections, points, n, -1.0) != 0) {
        goto feasible;
    }
    /* At this point a feasible solution must have z = 0. Fix y to 1, -1 and
       perform 1D feasibility tests. */
    if (_htm_feasible_1d(points, n, 1.0) != 0) {
        goto feasible;
    }
    if (_htm_feasible_1d(points, n, -1.0) != 0) {
        goto feasible;
    }
    /* At this point a feasible solution must have y = z = 0. If all x
       coordinates are non-zero and of the same sign, then there is a
       feasible solution. */
    for (i = 0, pos = 0, neg = 0; i < n; ++i) {
        if (points[i].x > 0.0) {
            if (neg != 0) {
                goto not_feasible;
            }
            pos = 1;
        } else if (points[i].x < 0.0) {
            if (pos != 0) {
                goto not_feasible;
            }
            neg = 1;
        } else {
            goto not_feasible;
        }
    }

feasible:
    if (buf != stackbuf) {
        free(buf);
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return 1;

not_feasible:
    if (buf != stackbuf) {
        free(buf);
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    return 0;
}


/* ---- Convexity test ---- */

/*  Square norm of the robust cross product of 2 cartesian unit vectors must be
    >= HTM_RCROSS_N2MIN, or the edge joining them is considered degenerate.
 */
static const double HTM_RCROSS_N2MIN = 4.0e-16;

/*  Dot product of a unit plane normal and a cartesian unit vector must be
    > HTM_SIN_MIN, or the vector is considered to be on the plane.
 */
static const double HTM_SIN_MIN = 1.0e-10;

/*  Dot product of 2 cartesian unit vectors must be < HTM_COS_MAX
    or they are considered degenerate.
*/
static const double HTM_COS_MAX = 0.999999999999999;


int htm_v3_convex(const struct htm_v3 *points,
                  size_t n,
                  enum htm_errcode *err)
{
    struct htm_v3 cen = { 0.0, 0.0, 0.0 };
    struct htm_v3 plane, p1, p2;
    double d, n2, wind = 0.0;
    size_t end;
    int cw = 0, ccw = 0;
    enum htm_errcode e;

    if (n < 3) {
        if (err != NULL) {
            *err = HTM_ELEN;
        }
        return 0;
    }
    if (!htm_v3_hemispherical(points, n, &e)) {
        if (err != NULL) {
            *err = e;
        }
        return 0;
    }
    htm_v3_centroid(&cen, points, n);
    htm_v3_rcross(&p1, &cen, &points[n - 1]);
    n2 = htm_v3_norm2(&p1);
    if (fabs(n2) < HTM_RCROSS_N2MIN) {
        if (err != NULL) {
            *err = HTM_EDEGEN;
        }
        return 0;
    }
    for (end = 0; end < n; ++end) {
        size_t beg = (end < 2 ? (n - 2) + end : end - 2);
        size_t mid = (end == 0 ? n - 1 : end - 1);
        htm_v3_rcross(&plane, &points[mid], &points[end]);
        n2 = htm_v3_norm2(&plane);
        if (htm_v3_dot(&points[mid], &points[end]) >= HTM_COS_MAX ||
            n2 < HTM_RCROSS_N2MIN) {
            if (err != NULL) {
                *err = HTM_EDEGEN;
            }
            return 0;
        }
        htm_v3_div(&plane, &plane, sqrt(n2));
        d = htm_v3_dot(&plane, &points[beg]);
        if (d > HTM_SIN_MIN) {
            if (cw) {
                if (err != NULL) {
                    *err = HTM_OK;
                } 
                return 0;
            }
            ccw = 1;
        } else if (d < - HTM_SIN_MIN) {
            if (ccw) {
                if (err != NULL) {
                    *err = HTM_OK;
                }
                return 0;
            }
            cw = 1;
        }
        /* check that vertices always wind around cen in the same direction */
        d = htm_v3_dot(&plane, &cen);
        if ((d < HTM_SIN_MIN && ccw) || (d > -HTM_SIN_MIN && cw)) {
            if (err != NULL) {
                *err = HTM_OK;
            }
            return 0;
        }
        /* sum up winding angle for edge (mid, end) */
        htm_v3_rcross(&p2, &cen, &points[end]);
        n2 = htm_v3_norm2(&p2);
        if (fabs(n2) < HTM_RCROSS_N2MIN) {
            if (err != NULL) {
                *err = HTM_EDEGEN;
            }
            return 0;
        }
        wind += htm_v3_angsep(&p1, &p2);
        p1 = p2;
    }
    if (err != NULL) {
        *err = HTM_OK;
    }
    /* for convex polygons, the closest multiple of 360 to the
       total winding angle is 1 */
    if (wind > 180.0 && wind < 540.0) {
       return ccw ? 1 : -1;
    }
    return 0;
}


/* ---- Convex hull algorithm ---- */

/** An (angle, vertex) pair.
  */
struct _htm_av3 {
    double angle;
    struct htm_v3 v;
} HTM_ALIGNED(16);


static void _htm_av3_insertsort(struct _htm_av3 *elts, size_t n)
{
    struct _htm_av3 k;
    size_t i, j;

    for (j = 1; j < n; ++j) {
        k = elts[j];
        for (i = j; i > 0 && elts[i - 1].angle > k.angle; --i) {
            elts[i] = elts[i - 1];
        }
        elts[i] = k;
    }
}


static void _htm_av3_merge(struct _htm_av3 *dst,
                           const struct _htm_av3 *left,
                           size_t nleft,
                           const struct _htm_av3 *right,
                           size_t nright)
{
    while (nleft > 0 && nright > 0) {
        if (left->angle <= right->angle) {
            *dst = *left;
            ++left;
            --nleft;
        } else {
            *dst = *right;
            ++right;
            --nright;
        }
        ++dst;
    }
    if (nleft > 0) {
        memcpy(dst, left, nleft * sizeof(struct _htm_av3));
    } else if (nright > 0) {
        memcpy(dst, right, nright * sizeof(struct _htm_av3));
    }
}


static void _htm_av3_mergesort(struct _htm_av3 *elts, size_t n)
{
    size_t i, ns;
    uint64_t clg2n;
    if (n <= 8) {
        _htm_av3_insertsort(elts, n);
        return;
    }
    /* bottom up merge sort: in-place insertion sort, followed by a tree
       of merges.

       The merge tree is evaluated breadth-first, allowing a simple,
       non-recursive loop-based implementation. Assume that elts has
       space for 2*n entries, providing the scratch space required by the
       mergesort algorithm. Each level of the merge tree m moves the
       first/last n entries to the last/first n entries of elts. In order
       to end up with the sorted array in the first n entries, arrange
       to always perform an even number of merges by adjusting the size
       of the insertion sort performed at the merge tree leaves. */
    clg2n = (uint64_t) (n - 1);
    clg2n |= (clg2n >> 1);
    clg2n |= (clg2n >> 2);
    clg2n |= (clg2n >> 4);
    clg2n |= (clg2n >> 8);
    clg2n |= (clg2n >> 16);
    clg2n |= (clg2n >> 32);
    clg2n = htm_popcount(clg2n);
    if ((clg2n & 1) == 0) {
        ns = 4;
        clg2n -= 2;
    } else {
        ns = 8;
        clg2n -= 3;
    }
    /* insertion sort pass */
    for (i = 0; i < n; i += ns) {
        _htm_av3_insertsort(elts + i, (n - i >= ns) ? ns : n - i);
    }
    /* evaluate merge-tree breadth first */
    for (; clg2n > 0; --clg2n, ns *= 2) {
        struct _htm_av3 *src = elts + ((clg2n & 1) == 0 ? 0 : n);
        struct _htm_av3 *scratch = elts + ((clg2n & 1) == 0 ? n : 0);
        for (i = 0; i + 2*ns < n; i += 2*ns) {
            _htm_av3_merge(&scratch[i], &src[i], ns, &src[i + ns], ns);
        }
        if (n - i > ns) {
            _htm_av3_merge(&scratch[i], &src[i], ns, &src[i + ns], n - i - ns);
        } else {
            memcpy(&scratch[i], &src[i], (n - i) * sizeof(struct _htm_av3));
        }
    }
}


struct htm_s2cpoly * htm_s2cpoly_hull(const struct htm_v3 *points,
                                      size_t n,
                                      enum htm_errcode *err)
{
    struct _htm_av3 stackbuf[128];
    struct _htm_av3 *av;
    struct htm_s2cpoly *poly = NULL;
    struct htm_v3 *anchor;
    struct htm_v3 *edge = NULL;
    struct htm_v3 center, refplane;
    size_t i, extremum, nav, ncv;
    enum htm_errcode e = HTM_OK;

    if (points == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return NULL;
    }
    if (!_htm_nv_valid(n)) {
        if (err != NULL) {
            *err = HTM_ELEN;
        }
        return NULL;
    }
    if (htm_v3_hemispherical(points, n, &e) != 1) {
        if (err != NULL) {
            *err = (e == HTM_OK) ? HTM_EHEMIS : e;
        }
        return NULL;
    }
    /* allocate space for 2n _htm_av3 structs. We need n vertices, n
       scratch entries for merge sort, and space for n edge plane normals. */
    if (n <= sizeof(stackbuf)/(2 * sizeof(struct _htm_av3))) {
        av = stackbuf;
    } else {
        av = (struct _htm_av3 *) malloc(2 * sizeof(struct _htm_av3) * n);
        if (av == NULL) {
            if (err != NULL) {
                *err = HTM_ENOMEM;
            }
            return NULL;
        }
    }
    anchor = &av[0].v;

    /* point furthest from the center is on the hull. */
    htm_v3_centroid(&center, points, n);
    {
        double n2, maxsep = 0.0;
        extremum = 0;
        for (i = 0, maxsep = 0.0; i < n; ++i) {
            double sep = htm_v3_dist2(&points[i], &center);
            if (sep > maxsep) {
                extremum = i;
                maxsep = sep;
            }
        }
        *anchor = points[extremum];
        htm_v3_rcross(&refplane, &center, anchor);
        n2 = htm_v3_dot(&refplane, &refplane);
        if (n2 < HTM_RCROSS_N2MIN) {
            /* vertex and center are too close */
            e = HTM_EDEGEN;
            goto cleanup;
        }
        htm_v3_div(&refplane, &refplane, sqrt(n2));
    }

    /* build (angle,vertex) list */
    av[0].angle = 0.0;
    nav = 1;
    for (i = 0; i < n; ++i) {
        struct htm_v3 plane;
        double n2;
        if (extremum == i) {
            continue;
        }
        htm_v3_rcross(&plane, &center, &points[i]);
        n2 = htm_v3_norm2(&plane);
        /* skip points too close to center */
        if (n2 >= HTM_RCROSS_N2MIN) {
            struct htm_v3 p;
            double sa, angle;

            htm_v3_div(&plane, &plane, sqrt(n2));
            htm_v3_rcross(&p, &refplane, &plane);
            sa = htm_v3_norm(&p);
            if (htm_v3_dot(&p, &center) < 0.0) {
                sa = -sa;
            }
            angle = atan2(sa, htm_v3_dot(&refplane, &plane));
            if (angle < 0.0) {
                angle += 2.0*M_PI;
            }
            av[nav].angle = angle;
            av[nav].v = points[i];
            ++nav;
        }
    }
    if (nav < 3) {
        e = HTM_EDEGEN;
        goto cleanup;
    }

    /* order points by winding angle from the first (extreme) vertex */
    _htm_av3_mergesort(av, nav);
    /* stable: av[0].v == anchor */

    /* Loop over vertices using a Graham scan adapted for spherical geometry.
       Store verticess in

           av[0].v, av[1].v, ...

       and edges in

           av[nav].v av[nav + 1].v ...
     */
    edge = NULL;
    ncv = 1;
    for (i = 1; i < nav;) {
        struct htm_v3 p;
        struct htm_v3 *v = &av[i].v;
        double n2;

        htm_v3_rcross(&p, anchor, v);
        n2 = htm_v3_norm2(&p);

        if (htm_v3_dot(anchor, v) < HTM_COS_MAX && n2 >= HTM_RCROSS_N2MIN) {
            if (ncv == 1) {
                /* compute first edge */
                edge = &av[nav].v;
                htm_v3_div(edge, &p, sqrt(n2));
                anchor = &av[1].v;
                *anchor = *v;
                ++ncv;
            } else {
                double d = htm_v3_dot(v, edge);
                if (d > HTM_SIN_MIN) {
                    /* v is inside the edge defined by the last
                       2 vertices on the hull */
                    edge = &av[nav + ncv - 1].v;
                    htm_v3_div(edge, &p, sqrt(n2));
                    anchor = &av[ncv].v;
                    *anchor = *v;
                    ++ncv;
                } else if (d < - HTM_SIN_MIN) {
                    /* backtrack - the most recently added hull vertex
                       is not actually on the hull. */
                    --ncv;
                    anchor = &av[ncv - 1].v;
                    edge = &av[nav + ncv - 2].v;
                    /* reprocess v to decide whether to add it to the hull
                       or whether further backtracking is necessary. */
                    continue;
                }
                /* v is coplanar with edge, skip it */
            }
        }
        ++i;
    }

    /* handle backtracking necessary for last edge */
    while (1) {
        struct htm_v3 p;
        struct htm_v3 *v = &av[0].v;
        double n2;

        if (ncv < 3) {
            e = HTM_EDEGEN;
            goto cleanup;
        }
        htm_v3_rcross(&p, anchor, v);
        n2 = htm_v3_norm2(&p);
        if (htm_v3_dot(anchor, v) < HTM_COS_MAX && n2 >= HTM_RCROSS_N2MIN) {
            if (htm_v3_dot(v, edge) > HTM_SIN_MIN) {
                htm_v3_div(&av[nav + ncv - 1].v, &p, sqrt(n2));
                break;
            }
        }
        --ncv;
        anchor = &av[ncv - 1].v;
        edge = &av[nav + ncv - 1].v;
    }

    /* allocate and populate convex polygon */
    poly = (struct htm_s2cpoly *) malloc(
        sizeof(struct htm_s2cpoly) + 2 * ncv * sizeof(struct htm_v3));
    if (poly == NULL) {
        e = HTM_ENOMEM;
        goto cleanup;
    }
    center.x = center.y = center.z = 0.0;
    for (i = 0; i < ncv; ++i) {
        htm_v3_add(&center, &center, &av[i].v);
        poly->ve[i] = av[i].v;
        poly->ve[i + ncv] = av[i + nav].v;
    }
    poly->n = ncv;
    poly->vsum = center;

cleanup:
    if (av != stackbuf) {
        free(av);
    }
    if (err != NULL) {
        *err = e;
    }
    return poly;
}


#ifdef __cplusplus
}
#endif

