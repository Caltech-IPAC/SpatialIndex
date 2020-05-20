/** \file
    \brief      Listing HTM triangles overlapping a region.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "tinyhtm/geometry.h"
#include "tinyhtm/htm.h"

/** \cond */

/* Should IDs be decimal encoded? */
static int decimal = 0;
/* Should IDs be output in range form? */
static int ranges = 0;
/* HTM subdivision level */
static int level = 0;
/* Maximum number of ranges */
static size_t maxranges = SIZE_MAX;


static void err(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stdout, "ERROR: ");
    vfprintf(stdout, fmt, ap);
    va_end(ap);
    fprintf(stdout, "\n");
    fflush(stdout);
    exit(EXIT_FAILURE);
}


static void get_level(const char *s)
{
    char *endptr;
    long l = strtol(s, &endptr, 10);
    if (errno != 0 || endptr == s || l < 0 ||
        l > (decimal ? HTM_DEC_MAX_LEVEL : HTM_MAX_LEVEL)) {
        err("HTM subdivision level `%s' is non-integeral, negative or too large", s);
    }
    level = (int) l;
}


static void get_maxranges(const char *s)
{
    char *endptr;
    unsigned long long mr = strtoull(s, &endptr, 10);
    if (errno != 0 || endptr == s || mr < 4 || mr > SIZE_MAX) {
        err("Maximum range count `%s' is non-integeral, less than 4, or too large", s);
    }
    maxranges = (size_t) mr;
}


static double get_double(const char *s)
{
    char *endptr;
    double d = strtod(s, &endptr);
    if (errno != 0 || endptr == s) {
        err("failed to convert argument `%s' to a double", s);
    }
    return d;
}


static void print_ids(const struct htm_ids * const ids)
{
    for (size_t i = 0; i < ids->n; ++i) {
        if (ranges) {
            int64_t min = ids->range[i].min;
            int64_t max = ids->range[i].max;
            printf("%lld %lld\n",
                   (long long) (decimal ? htm_idtodec(min) : min),
                   (long long) (decimal ? htm_idtodec(max) : max));
        } else {
            for (int64_t j = ids->range[i].min; j <= ids->range[i].max; ++j) {
                printf("%lld\n", (long long) (decimal ? htm_idtodec(j) : j));
            }
        }
    }
}


static void circle_list(char **argv)
{
    struct htm_sc sc;
    struct htm_v3 cen;
    double r;
    struct htm_ids *ids = NULL;
    enum htm_errcode ec;

    ec = htm_sc_init(&sc, get_double(argv[0]), get_double(argv[1]));
    if (ec != HTM_OK) {
        err("Invalid circle center coordinates: %s", htm_errmsg(ec));
    }
    ec = htm_sc_tov3(&cen, &sc);
    if (ec != HTM_OK) {
        err("Failed to convert spherical coordinates to a unit vector: %s",
            htm_errmsg(ec));
    }
    r = get_double(argv[2]);
    ids = htm_s2circle_ids(ids, &cen, r, level, maxranges, &ec);
    if (ec != HTM_OK) {
        err("Failed to find HTM triagles overlapping circle: %s",
            htm_errmsg(ec));
    }
    print_ids(ids);
    free(ids);
}


static void ellipse_list(char **argv)
{
    struct htm_s2ellipse ellipse;
    struct htm_sc sc;
    struct htm_v3 cen;
    double a, b, angle;
    struct htm_ids *ids = NULL;
    enum htm_errcode ec;

    ec = htm_sc_init(&sc, get_double(argv[0]), get_double(argv[1]));
    if (ec != HTM_OK) {
        err("Invalid ellipse center coordinates: %s", htm_errmsg(ec));
    }
    ec = htm_sc_tov3(&cen, &sc);
    if (ec != HTM_OK) {
        err("Failed to convert spherical coordinates to a unit vector: %s",
            htm_errmsg(ec));
    }
    a = get_double(argv[2]);
    b = get_double(argv[3]);
    angle = get_double(argv[4]);
    ec = htm_s2ellipse_init2(&ellipse, &cen, a, b, angle);
    if (ec != HTM_OK) {
        err("Invalid ellipse parameters: %s", htm_errmsg(ec));
    }
    ids = htm_s2ellipse_ids(ids, &ellipse, level, maxranges, &ec);
    if (ec != HTM_OK) {
        err("Failed to find HTM triagles overlapping ellipse: %s",
            htm_errmsg(ec));
    }
    print_ids(ids);
    free(ids);
}


static void hull_list(const int argc,
                      char **argv)
{
    struct htm_s2cpoly *poly = NULL;
    struct htm_sc sc;
    struct htm_v3 *verts = NULL;
    struct htm_ids *ids = NULL;
    int i;
    enum htm_errcode ec = HTM_OK;

    verts = (struct htm_v3 *) malloc(sizeof(struct htm_v3) * (size_t) argc/2);
    if (verts == NULL) {
        err(htm_errmsg(HTM_ENOMEM));
    }
    for (i = 0; i < argc/2; ++i) {
        ec = htm_sc_init(&sc, get_double(argv[2*i]), get_double(argv[2*i + 1]));
        if (ec != HTM_OK) {
            free(verts);
            err("Invalid vertex coordinates: %s", htm_errmsg(ec));
        }
        ec = htm_sc_tov3(&verts[i], &sc);
        if (ec != HTM_OK) {
            free(verts);
            err("Failed to convert spherical coordinates to a unit vector: %s",
                htm_errmsg(ec));
        }
    }
    poly = htm_s2cpoly_hull(verts, (size_t) i, &ec);
    free(verts);
    if (ec != HTM_OK) {
        err("Failed to compute convex hull: %s", htm_errmsg(ec));
    }
    ids = htm_s2cpoly_ids(ids, poly, level, maxranges, &ec);
    free(poly);
    if (ec != HTM_OK) {
        err("Failed to find HTM triangles overlapping convex hull: %s",
            htm_errmsg(ec));
    }
    print_ids(ids);
    free(ids);
}
