/** \file
    \brief      Counting points in a region with HTM tree indexes.

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
#include "tinyhtm/tree.h"

/** \cond */

/* Should output be in JSON format or in IPAC SVC format? */
static int json = 0;
/* Should the count be estimated, or determined exactly? */
static int estimate = 0;

/* Performs string escaping for the JSON/IPAC SVC formats. */
static const char * esc(const char *s) {
    static char buf[8192];
    char *e = buf;
    if (s == NULL) {
        return "null";
    } else {
        *e++ = '"';
        for (; *s != '\0' && e < buf + (sizeof(buf) - 2); ++s) {
            switch (*s) {
                case  '"': e[0] = '\\'; e[1] =  '"'; e += 2; break;
                case '\\': e[0] = '\\'; e[1] = '\\'; e += 2; break;
                case '\b': e[0] = '\\'; e[1] =  'b'; e += 2; break;
                case '\f': e[0] = '\\'; e[1] =  'f'; e += 2; break;
                case '\n': e[0] = '\\'; e[1] =  'n'; e += 2; break;
                case '\r': e[0] = '\\'; e[1] =  'r'; e += 2; break;
                case '\t': e[0] = '\\'; e[1] =  't'; e += 2; break;
                default:
                    if (*s > 0x1f && *s < 0x7f) {
                        *e++ = *s;
                    }
                    break;
            }
        }
        if (*s != '\0' && e == buf + (sizeof(buf) - 2)) {
            strcpy(buf + (sizeof(buf) - 6), " ...\"");
        } else {
            e[0] = '"';
            e[1] = '\0';
        }
    }
    return buf;
}


static void err(const char *fmt, ...)
{
    char buf[4095];
    int len;

    va_list ap;
    va_start(ap, fmt);
    len = vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    if (len >= (int) sizeof(buf)) {
        strcpy(buf + (sizeof(buf) - 5), " ...");
    }
    if (json) {
        printf("{\"stat\":\"ERROR\", \"msg\":%s}\n", esc(buf));
    } else {
        printf("[struct stat=\"ERROR\", msg=%s]\n", esc(buf));
    }
    fflush(stdout);
    exit(EXIT_FAILURE);
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


static void print_count(int64_t count)
{
    if (json) {
        printf("{\"stat\":\"OK\", \"count\":%lld}\n", (long long) count);
    } else {
        printf("[struct stat=\"OK\", count=\"%lld\"]\n", (long long) count);
    }
}

static void print_range(const struct htm_range *range)
{
    if (json) {
        printf("{\"stat\":\"OK\", \"min\":%lld, \"max\":%lld}\n",
               (long long) range->min, (long long) range->max);
    } else {
        printf("[struct stat=\"OK\", min=\"%lld\", max=\"%lld\"]\n",
               (long long) range->min, (long long) range->max);
    }
}

static void circle_count(const char * const treefile,
                         const char * const datafile,
                         char **argv)
{
    struct htm_tree tree;
    struct htm_sc sc;
    struct htm_v3 cen;
    double r;
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
    ec = htm_tree_init(&tree, treefile, datafile);
    if (ec != HTM_OK) {
        err("Failed to load tree and/or data file: %s", htm_errmsg(ec));
    }
    if (estimate != 0) {
        struct htm_range range = htm_tree_s2circle_range(&tree, &cen, r, &ec);
        htm_tree_destroy(&tree);
        if (ec != HTM_OK) {
            err("Failed to estimate points in circle: %s", htm_errmsg(ec));
        }
        print_range(&range);
    } else {
        int64_t count = htm_tree_s2circle_count(&tree, &cen, r, &ec);
        htm_tree_destroy(&tree);
        if (ec != HTM_OK) {
            err("Failed to count points in circle: %s", htm_errmsg(ec));
        }
        print_count(count);
    }
}


static void ellipse_count(const char * const treefile,
                          const char * const datafile,
                          char **argv)
{
    struct htm_tree tree;
    struct htm_s2ellipse ellipse;
    struct htm_sc sc;
    struct htm_v3 cen;
    double a, b, angle;
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
    ec = htm_tree_init(&tree, treefile, datafile);
    if (ec != HTM_OK) {
        err("Failed to load tree and/or data file: %s", htm_errmsg(ec));
    }
    if (estimate != 0) {
        struct htm_range range = htm_tree_s2ellipse_range(&tree, &ellipse, &ec);
        htm_tree_destroy(&tree);
        if (ec != HTM_OK) {
            err("Failed to estimate points in ellipse: %s", htm_errmsg(ec));
        }
        print_range(&range);
    } else {
        int64_t count = htm_tree_s2ellipse_count(&tree, &ellipse, &ec); 
        htm_tree_destroy(&tree);
        if (ec != HTM_OK) {
            err("Failed to count points in ellipse: %s", htm_errmsg(ec));
        }
        print_count(count);
    }
}


static void hull_count(const char * const treefile,
                       const char * const datafile,
                       const int argc,
                       char **argv)
{
    struct htm_tree tree;
    struct htm_s2cpoly *poly;
    struct htm_sc sc;
    struct htm_v3 *verts;
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
    if (poly == NULL || ec != HTM_OK) {
        free(verts);
        err("Failed to compute convex hull: %s", htm_errmsg(ec));
    }
    free(verts);
    ec = htm_tree_init(&tree, treefile, datafile);
    if (ec != HTM_OK) {
        free(poly);
        err("Failed to load tree and/or data file: %s", htm_errmsg(ec));
    }
    if (estimate != 0) {
        struct htm_range range = htm_tree_s2cpoly_range(&tree, poly, &ec);
        htm_tree_destroy(&tree);
        free(poly);
        if (ec != HTM_OK) {
            err("Failed to estimate points in hull: %s", htm_errmsg(ec));
        }
        print_range(&range);
    } else {
        int64_t count = htm_tree_s2cpoly_count(&tree, poly, &ec);
        htm_tree_destroy(&tree);
        free(poly);
        if (ec != HTM_OK) {
            err("Failed to count points in hull: %s", htm_errmsg(ec));
        }
        print_count(count);
    }
}


static void test_tree(const char * const treefile,
                      const char * const datafile,
                      char **argv)
{
    struct htm_tree tree;
    double r;
    unsigned long long i;
    enum htm_errcode ec;

    r = get_double(argv[0]);
    ec = htm_tree_init(&tree, treefile, datafile);
    if (ec != HTM_OK) {
        err("Failed to load tree and/or data file: %s", htm_errmsg(ec));
    }
    for (i = 0; i < tree.count; ++i) {
        int64_t c = htm_tree_s2circle_count(&tree, &tree.entries[i].v, r, &ec);
        if (c != 1) {
            err("Circle of radius %g around %llu:(%.18g, %.18g, %.18g) "
                "contains %lld points (expecting 1).", r, i,
                tree.entries[i].v.x, tree.entries[i].v.y, tree.entries[i].v.z,
                (long long) c);
        }
    }
    htm_tree_destroy(&tree);
}
