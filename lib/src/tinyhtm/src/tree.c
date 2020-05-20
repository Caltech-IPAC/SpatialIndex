/** \file
    \brief      HTM tree index implementation.

    For API documentation, see tree.h.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#include "tinyhtm/tree.h"

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

#include "tinyhtm/varint.h"

#ifdef __cplusplus
extern "C" {
#endif


enum htm_errcode htm_tree_init(struct htm_tree *tree,
                               const char * const treefile,
                               const char * const datafile)
{
    struct stat sb;
    const unsigned char *s;
    uint64_t off, count;
    int i;
    const size_t pagesz =  (size_t) sysconf(_SC_PAGESIZE);
    enum htm_errcode err = HTM_OK;

    /* set defaults */
    tree->leafthresh = 0;
    tree->count = 0;
    for (i = 0; i < 8; ++i) {
        tree->root[i] = NULL;
    }
    tree->entries = (const struct htm_tree_entry *) MAP_FAILED;
    tree->index = (const void *) MAP_FAILED;
    tree->indexsz = 0;
    tree->datasz = 0;
    tree->indexfd = -1;
    tree->datafd = -1;

    /* check inputs */
    if (tree == NULL || datafile == NULL) {
        return HTM_ENULLPTR;
    }
    if (stat(datafile, &sb) != 0) {
        return HTM_EIO;
    }
    if (sb.st_size % sizeof(struct htm_tree_entry) != 0 || sb.st_size == 0) {
        return HTM_EINV;
    }
    count = (uint64_t) sb.st_size / sizeof(struct htm_tree_entry);

    /* memory map datafile */
    tree->datasz = (size_t) sb.st_size;
    if (tree->datasz % pagesz != 0) {
        tree->datasz += pagesz - tree->datasz % pagesz;
    }
    tree->datafd = open(datafile, O_RDONLY);
    if (tree->datafd == -1) {
        err = HTM_EIO;
        goto cleanup;
    }
    tree->entries = (const struct htm_tree_entry *) mmap(
        NULL, tree->datasz, PROT_READ, MAP_SHARED | MAP_NORESERVE,
        tree->datafd, 0);
    if ((void *) tree->entries == MAP_FAILED) {
        err = HTM_EMMAN;
        goto cleanup;
    }
    if (madvise((void *) tree->entries, tree->datasz, MADV_RANDOM) != 0) {
        err = HTM_EMMAN;
        goto cleanup;
    }

    /* memory map treefile (if there is one) */
    if (treefile == NULL) {
        tree->count = count;
        return HTM_OK;
    }
    if (stat(treefile, &sb) != 0) {
        err = HTM_EIO;
        goto cleanup;
    }
    tree->indexsz = (size_t) sb.st_size;
    if (tree->indexsz % pagesz != 0) {
        tree->indexsz += pagesz - tree->indexsz % pagesz;
    }
    tree->indexfd = open(treefile, O_RDONLY);
    if (tree->indexfd == -1) {
        err = HTM_EIO;
        goto cleanup;
    }
    tree->index = (const void *) mmap(
        NULL, tree->indexsz, PROT_READ, MAP_SHARED | MAP_NORESERVE,
        tree->indexfd, 0);
    if ((void *) tree->index == MAP_FAILED) {
        err = HTM_EMMAN;
        goto cleanup;
    }
    if (madvise((void *) tree->index, tree->indexsz, MADV_RANDOM) != 0) {
        err = HTM_EMMAN;
        goto cleanup;
    }

    /* parse tree file header */
    s = (const unsigned char *) tree->index;
    tree->leafthresh = htm_varint_decode(s);
    s += 1 + htm_varint_nfollow(*s);
    tree->count = htm_varint_decode(s);
    s += 1 + htm_varint_nfollow(*s);
    if (tree->count != count) {
        /* tree index point count does not agree with data file */
        err = HTM_ETREE;
        goto cleanup;
    }
    for (i = 0; i < 8; ++i) {
        off = htm_varint_decode(s);
        s += 1 + htm_varint_nfollow(*s);
        if (off == 0) {
            tree->root[i] = NULL;
        } else {
            tree->root[i] = s + (off - 1);
        }
    }
    if (s - (const unsigned char *) tree->index >= sb.st_size) {
        /* header overflowed tree file size */
        err = HTM_ETREE;
        goto cleanup;
    }
    return HTM_OK;

cleanup:
    htm_tree_destroy(tree);
    return err;
}


void htm_tree_destroy(struct htm_tree *tree)
{
    int i;
    if (tree == NULL) {
        return;
    }
    /* unmap and close data file */
    if ((void *) tree->entries != MAP_FAILED) {
        munmap((void *)tree->entries, tree->datasz);
        tree->entries = (const struct htm_tree_entry *) MAP_FAILED;
    }
    tree->datasz = 0;
    if (tree->datafd != -1) {
        close(tree->datafd);
        tree->datafd = -1;
    }
    /* unmap and close tree file */
    if (tree->index != MAP_FAILED) {
        munmap((void *)tree->index, tree->indexsz);
        tree->index = (const void *) MAP_FAILED;
    }
    tree->indexsz = 0;
    if (tree->indexfd != -1) {
        close(tree->indexfd);
        tree->indexfd = - 1;
    }
    /* set remaining fields to default values */
    tree->leafthresh = 0;
    tree->count = 0;
    for (i = 0; i < 8; ++i) {
        tree->root[i] = NULL;
    }
}


enum htm_errcode htm_tree_lock(struct htm_tree *tree, size_t datathresh)
{
    if (tree == NULL) {
        return HTM_ENULLPTR;
    }
    if (tree->indexfd != -1) {
        if (mlock(tree->index, tree->indexsz) != 0) {
            return HTM_ENOMEM;
        }
    }
    if (tree->datasz <= datathresh) {
        if (mlock(tree->entries, tree->datasz) != 0) {
            return HTM_ENOMEM;
        }
    }
    return HTM_OK;
}


int64_t htm_tree_s2circle_scan(const struct htm_tree *tree,
                               const struct htm_v3 *center,
                               double radius,
                               enum htm_errcode *err)
{
    double dist2;
    int64_t count;
    uint64_t i;

    if (tree == NULL || center == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return -1;
    }
    if (radius < 0.0) {
        return 0;
    } else if (radius >= 180.0) {
        return (int64_t) tree->count;
    }
    dist2 = sin(radius * 0.5 * HTM_RAD_PER_DEG);
    dist2 = 4.0 * dist2 * dist2;
    count = 0;
    for (i = 0, count = 0; i < tree->count; ++i) {
        if (htm_v3_dist2(center, &tree->entries[i].v) <= dist2) {
            ++count;
        }
    }
    return count;
}


int64_t htm_tree_s2ellipse_scan(const struct htm_tree *tree,
                                const struct htm_s2ellipse *ellipse,
                                enum htm_errcode *err)
{
    int64_t count;
    uint64_t i;

    if (tree == NULL || ellipse == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return -1;
    }
    count = 0;
    for (i = 0, count = 0; i < tree->count; ++i) {
        if (htm_s2ellipse_cv3(ellipse, &tree->entries[i].v) != 0) {
            ++count;
        }
    }
    return count;
}


int64_t htm_tree_s2cpoly_scan(const struct htm_tree *tree,
                              const struct htm_s2cpoly *poly,
                              enum htm_errcode *err)
{
    int64_t count;
    uint64_t i;

    if (tree == NULL || poly == NULL) {
        if (err != NULL) {
            *err = HTM_ENULLPTR;
        }
        return -1;
    }
    count = 0;
    for (i = 0, count = 0; i < tree->count; ++i) {
        if (htm_s2cpoly_cv3(poly, &tree->entries[i].v) != 0) {
            ++count;
        }
    }
    return count;
}


#ifdef __cplusplus
}
#endif

