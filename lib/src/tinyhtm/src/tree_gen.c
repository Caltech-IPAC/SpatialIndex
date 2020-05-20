/** \file
    \brief      HTM tree generation.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */

/** \page htm_tree_gen  Tree index generation with htm_tree_gen

    The \c htm_tree_gen utility takes character separated value files
    containing points (integer ID, longitude angle and latitude angle
    columns are expected), and outputs an HTM tree index which is useful
    for fast point-in-region counts.

    \section usage Usage

    <pre><tt>
    htm_tree_gen [options] &lt;tree file&gt; &lt;input file 1&gt; [&lt;input file 2&gt; ...]
    </tt></pre>

    \section opts Command Line Options

    Running <tt>htm_tree_index --help</tt> provides a list of the supported
    command line options and their descriptions.
    
    \section over Overview

    HTM tree indexes can be used to quickly count or estimate how many points
    fall in a spherical region.

    The indexing utility produces two binary files - a data file consisting
    of points sorted in ascending HTM ID order, and a tree file which contains
    an HTM tree over the points. Tree nodes contain a count of the number of
    points inside the node and an index into the data file; empty nodes are
    omitted.

    This leads to the following simple and fast counting algorithm:

    Given a geometric region, visit all non-empty nodes overlapping the region.
    For nodes fully covered by the region, add the point count to a running
    total. For partially covered internal nodes, visit the non-empty
    children. For partially covered leaves, read points from the data
    file and check whether they are contained in the region.

    \section data Data File Algorithm

    Producing the sorted data file is conceptually simple; all that is
    required is an external sorting routine. The implementation uses
    quicksort on in-memory blocks followed by a series of k-way merges.

    \section tree Tree File Algorithm

    More difficult is an external algorithm for producing the tree file,
    especially if one wishes to optimize with respect to cache-line and
    memory page size.

    Alstrup's Split-and-Refine algorithm is used for layout. It produces
    a cache oblivious layout via repeated application of an arbitrary black-box
    algorithm, which is required to produce an optimal layout for a specified
    block size. The black-box is Clark and Munro's greedy bottom-up
    method for worst-cast optimal tree layout.

    These algorithms are detailed in the papers below:

    <pre>
    David R. Clark and J. Ian Munro.
    Efficient suffix trees on secondary storage.
    In Proceedings of the 7th Annual ACM-SIAM Symposium on Discrete Algorithms,
    Pages 383-391, Atlanta, January 1996.
    </pre>

    <pre>
    Stephen Alstrup, Michael A. Bender, Erik D. Demaine, Martin Farach-Colton,
    Theis Rauhe, and Mikkel Thorup.
    Efficient tree layout in a multilevel memory hierarchy.
    arXiv:cs.DS/0211010, 2004.
    http://www.arXiv.org/abs/cs.DS/0211010.
    </pre>

    The tree generation phase produces an HTM tree with the following
    properties:

        - Each node N corresponds to an HTM triangle, and stores
          a count of the number of entries inside N, count(N), as well
          as the index (in the data file) of its first entry.
        - If level(N) == 20, N is a leaf - no children are generated.
        - If count(N) < K for some fixed threshold K,
          N is a leaf - no children are generated.
        - If count(N) >= K and level(N) < 20, N is an internal node - the
          non-empty children of N are generated.

    Given the data file (points sorted by HTM ID), it is possible to emit
    tree nodes satisfying the above in post-order with a single sequential
    scan.

    Iterating over nodes in post-order means that a parent node is processed
    only after all children have been processed. Because Clark and Munro's
    layout algorithm is bottom-up, one can simultaneously perform layout for
    all the block-sizes used by Split-and-Refine.

    The result is a string of block IDs for each node - a unique node ID.
    Note that block IDs are handed out sequentially for a given block size.
    Since tree nodes are visited in post-order, the block ID of a parent
    must be greater than or equal to the block ID of a child. Once block IDs 
    for a node have been computed at every block size, the node can be written
    to disk. If the largest block size is B, the in-memory node size is M,
    and the on disk node size is N, then 8BM/N bytes of RAM are required
    in the worst case.

    Next, the tree node file is sorted by node ID, yielding the desired
    cache-oblivious layout (per Alstrup). The sort is in ascending
    lexicographic block ID string order, always placing children before
    parents.

    At this stage, nodes contain very space-heavy node IDs, and child pointers
    are also node IDs! These are necessary for layout, but are not useful for
    tree searches. In fact, mapping a node ID to a file offset is non-trivial.

    Therefore, the sorted node file is converted to a much more space
    efficient representation:

        - Leaf nodes store a position count and the data file index of the
          first position inside the node. The index is stored as an offset
          relative to the index of the parent (this results in smaller
          numbers towards the leaves, allowing variable length coding to
          squeeze them into fewer bytes - see below).
        - Internal nodes additionally store relative tree file offsets for
          4 children. Tree file offsets of empty children are set to 0.
        - Note that internal nodes can be distinguished from leaves simply
          by comparing their position count to the leaf threshold K.

    Variable length coding is used for counts and data file indexes. Using
    such an encoding for tree file offsets is hard. Why? Because with variable
    length encoding, the size of a child offset depends on the offset value,
    which in turn depends on the size of the child offsets of nodes between
    the parent and child - in other words, it is node order dependent.
    This invalidates the block size computations required by Clark and
    Munro's algorithm, which has no a priori knowledge of the final node
    order - the reason for invoking it in the first place is to determine
    that order! Nevertheless, offsets are also variable length coded since
    it significantly reduces tree file size. The layout algorithm simply
    treats child offsets as having some fixed (hopefully close to average)
    on-disk size.

    The main difficulty during compression is the computation of tree file
    offsets from node IDs. Note that in the sorted node file a child will
    always occur before a parent. Therefore, the algorithm operates as
    follows: first, an empty hash-table that maps node IDs to file offsets
    is initialized. Next, the sorted node file is scanned sequentially.
    Finally, for each node:

        - the offsets of its children are looked up in the hashtable and the
          corresponding hashtable entries are removed (a node has exactly one
          parent).
        - the compressed byte string corresponding to the node is written
          out in reverse, and the size of the node is added to the running
          total file size.
        - an entry mapping the node id to the file size is added to the
          hashtable. To be precise: the parent of a node can subtract the
          file size just after a child was written from the current file
          size to obtain the child's offset.

    This yields a a compressed tree file T in reverse byte order. The last
    step consists of simply reversing the order of the bytes in T, giving
    the final tree file (with parent nodes always occurring before children).

    \section notes Notes

    The Clark and Munro layout algorithm is invoked only with a limited set
    of block sizes. The small number chosen (5) results in various internal
    structures having reasonable power-of-2 sizes, whilst still covering the
    block sizes that commonly occur in the multi-level memory hierarchies of
    modern 64-bit x86 systems. The largest block size used is 2 MiB, which
    corresponds to the size of a large page on x86-64. Using larger sizes
    seems to be of questionable utility and would increase the memory
    requirements for the tree generation phase, as nodes can only be written
    once a block ID for every block size has been assigned.

    The code is targeted at 64-bit machines with plenty of virtual address
    space. For example, huge files are mmapped in read-only mode (with
    MAP_NORESERVE). If system admins are restricting per-process address
    space, one may have to run `ulimit -v unlimited` prior to htm_tree_gen.

    The largest table that will have to be dealt with in the near future will
    contain around 10 billion points; at 32 bytes per record, this will
    involve asking for 320GB of address space. Should this be problematic,
    fancier IO strategies may have to be employed. Note however that e.g.
    RHEL 6 (and other recent x86-64 Linux distros) have a per process virtual
    address space limit of 128TB, and that in any case, practical tree
    generation at larger scales will require parallelization, possibly even
    across machines.

    While such parallelization is conceptually straightforward, even for
    tree generation/compression, it adds very significant implementation
    complexity, and so the code is mostly serial for now.
  */
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <getopt.h>
#include <pthread.h>

#include "tinyhtm/geometry.h"
#include "tinyhtm/htm.h"
#include "tinyhtm/tree.h"
#include "tinyhtm/varint.h"


/** \cond */

/* ================================================================ */
/*                     Utilities and plumbing                       */
/* ================================================================ */


/*  Prints an error message to stderr and exits.
 */
static void err(const char *fmt, ...)
{
    va_list ap;
    fprintf(stderr, "ERROR: ");
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, "\n");
    fflush(stderr);
    exit(EXIT_FAILURE);
}

/*  Prints a message to stdout and exits.
 */
static void msg(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);
    fflush(stdout);
}

/*  Returns the number of seconds that have elapsed since the epoch.
 */
static double now() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return ((double) t.tv_sec) + ((double) t.tv_usec) / 1.0e6;
}


/* ---- Memory and IO related parameters ---- */

struct mem_params {
    size_t memsz;   /* Max memory usage in bytes */
    size_t sortsz;  /* Size of blocks operated on by in-memory sorts */
    size_t ioblksz; /* Size of IO blocks */
    size_t k;       /* Number of merge segments in one multi-way merge pass */
};

static void mem_params_init(struct mem_params *mem,
                            size_t total,
                            size_t ioblksz)
{
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);
    if (total % (2*pagesz) != 0) {
        total += 2*pagesz - total % (2*pagesz);
    }
    if (ioblksz % pagesz != 0) {
        ioblksz += pagesz - ioblksz % pagesz;
    }
    if (total < 6*ioblksz) {
        total = 6*ioblksz;
    }
    mem->memsz = total;
    mem->sortsz = total / 2;
    mem->ioblksz = ioblksz;
    mem->k = (total - 2*ioblksz) / (2*ioblksz);
}


/* ---- Asynchronous block writer ---- */

enum bk_write_state {
    BK_WRITE_START = 0,
    BK_WRITE_READY,
    BK_WRITE_BLK,
    BK_WRITE_EXIT,
    BK_WRITE_ERR
};


struct blk_writer {
    size_t n;           /* Number of bytes per block */
    size_t i;           /* Number of bytes in current block */
    unsigned char *buf; /* Current block pointer */
    void *mem;          /* Space for 2 blocks of n tree entries */
    size_t memsz;       /* Size of mem in bytes */
    int fd;             /* File descriptor of file being written to */
    enum bk_write_state state;
    unsigned char *wrbuf;
    size_t wrbytes;
    pthread_attr_t attr;
    pthread_mutex_t mtx;
    pthread_cond_t cv;
    pthread_t thr;
};


/*  Writes data blocks in the background.
 */
static void * bk_write(void *arg)
{
    struct blk_writer *w = (struct blk_writer *) arg;
    pthread_mutex_lock(&w->mtx);
    while (1) {
        unsigned char *buf;
        size_t n;
        ssize_t b;
        /* signal readiness for a command */
        w->state = BK_WRITE_READY;
        pthread_cond_signal(&w->cv);
        /* wait for a command */
        pthread_cond_wait(&w->cv, &w->mtx);
        if (w->state == BK_WRITE_READY) {
            continue; /* nothing to do */
        } else if (w->state == BK_WRITE_EXIT) {
            break; /* exit background writer thread */
        }
        /* write a block */
        buf = w->wrbuf;
        n = w->wrbytes;
        pthread_mutex_unlock(&w->mtx);
        while (n > 0) {
            b = write(w->fd, buf, n);
            if (b < 0) {
                if (errno != EINTR) {
                    /* error - exit background writer thread */
                    pthread_mutex_lock(&w->mtx);
                    w->state = BK_WRITE_ERR;
                    goto done;
                }
            } else {
                buf += b;
                n -= (size_t) b;
            }
        }
        pthread_mutex_lock(&w->mtx);
    }
done:
    pthread_cond_signal(&w->cv);
    pthread_mutex_unlock(&w->mtx);
    return NULL;
}


/*  Creates a new block writer with the given block size.
 */
static struct blk_writer * blk_writer_init(const char * const file,
                                           size_t blksz)
{
    struct blk_writer *w;
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);

    w = (struct blk_writer *) malloc(sizeof(struct blk_writer));
    if (w == NULL) {
        err("malloc() failed");
    }
    w->n = blksz;
    w->i = 0;
    w->fd = open(file, O_CREAT | O_TRUNC | O_APPEND | O_WRONLY,
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
    if (w->fd == -1) {
        err("failed to open file %s for writing", file);
    }
    w->memsz = 2 * blksz;
    if (w->memsz % pagesz != 0) {
        w->memsz += pagesz - w->memsz % pagesz;
    }
    w->mem = mmap(NULL, w->memsz, PROT_READ | PROT_WRITE,
                  MAP_ANON | MAP_PRIVATE, -1, 0);
    if (w->mem == MAP_FAILED) {
        err("write buffer allocation via mmap() failed");
    }
    w->buf = (unsigned char *) w->mem;
    w->state = BK_WRITE_START;
    w->wrbuf = NULL;
    w->wrbytes = 0;
    pthread_attr_init(&w->attr);
    pthread_attr_setdetachstate(&w->attr, PTHREAD_CREATE_JOINABLE);
    pthread_mutex_init(&w->mtx, NULL);
    pthread_cond_init(&w->cv, NULL);
    pthread_create(&w->thr, &w->attr, &bk_write, (void *) w);
    return w;
}


/*  Issues (asynchronous) write of one block.
 */
static void blk_writer_issue(struct blk_writer * const w,
                             void (*sortfn)(void *, size_t))
{
    if (sortfn != NULL) {
        (*sortfn)(w->buf, w->i);
    }
    pthread_mutex_lock(&w->mtx);
    /* wait until background writer thread is ready for a command */
    while (w->state != BK_WRITE_READY && w->state != BK_WRITE_ERR) {
        pthread_cond_wait(&w->cv, &w->mtx);
    }
    if (w->state == BK_WRITE_ERR) {
        pthread_mutex_unlock(&w->mtx);
        err("background thread failed to write() disk block");
    }
    /* issue write for the current block */
    w->state = BK_WRITE_BLK;
    w->wrbuf = w->buf;
    w->wrbytes = w->i;
    pthread_cond_signal(&w->cv);
    pthread_mutex_unlock(&w->mtx);
    /* flip write buffer */
    w->i = 0;
    if (w->buf == (unsigned char *) w->mem) {
        w->buf = ((unsigned char *) w->mem) + w->n;
    } else {
        w->buf = ((unsigned char *) w->mem);
    }
}


/*  Closes writer w; all as yet unwritten data is flushed to disk first.

    If a non-NULL sorting function pointer is supplied, unwritten
    data is sorted prior to being written.
 */
static void blk_writer_close(struct blk_writer * const w,
                             void (*sortfn)(void *, size_t))
{
    blk_writer_issue(w, sortfn);
    pthread_mutex_lock(&w->mtx);
    /* wait until background writer thread is ready for a command */
    while (w->state != BK_WRITE_READY && w->state != BK_WRITE_ERR) {
        pthread_cond_wait(&w->cv, &w->mtx);
    }
    if (w->state == BK_WRITE_ERR) {
        pthread_mutex_unlock(&w->mtx);
        err("background thread failed to write() disk block");
    }
    /* issue thread exit command ... */
    w->state = BK_WRITE_EXIT;
    pthread_cond_signal(&w->cv);
    pthread_mutex_unlock(&w->mtx);
    /* ... and wait for writer thread to terminate */
    pthread_join(w->thr, NULL);
    /* clean up */
    if (munmap(w->mem, w->memsz) != 0) {
        err("munmap() of write buffers failed");
    }
    if (fsync(w->fd) != 0) {
        err("fsync() failed");
    }
    if (close(w->fd) != 0) {
        err("file close() failed");
    }
    pthread_cond_destroy(&w->cv);
    pthread_mutex_destroy(&w->mtx);
    pthread_attr_destroy(&w->attr);
    free(w);
}


/*  Appends data to the given block writer.
 */
HTM_INLINE void blk_writer_append(struct blk_writer * const w,
                                  const void * const data,
                                  const size_t nbytes,
                                  void (*sortfn)(void *, size_t))
{
    size_t i = w->i;
    if (i + nbytes > w->n) {
        blk_writer_issue(w, sortfn);
        i = 0;
    }
    assert(i + nbytes <= w->n);
    memcpy(w->buf + i, data, nbytes);
    w->i = i + nbytes;
}


/* ---- Multi-way external merge sort ---- */

/*  A contiguous sequence of sorted items, dubbed a multi-way merge segment.
 */
struct mrg_seg {
    const void *cur;
    const void *end;
    const void *blk;
};


/*  Initializes a merge segment.
 */
static void mrg_seg_init(struct mrg_seg * const s,
                         const void * start,
                         const void * end,
                         const size_t blksz)
{
    size_t n;
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);
    size_t nbytes = (const char *) end - (const char *) start;
    if (end <= start || blksz % pagesz != 0) {
        err("Invalid merge segment");
    }
    s->cur = start;
    s->end = end;
    n = ((size_t) start) % pagesz;
    if (n != 0) {
        start = (const char *) start - n;
        nbytes += n;
    }
    n = (nbytes > 2*blksz ? 2*blksz : nbytes);
    if (n % pagesz != 0) {
        n += pagesz - n % pagesz;
    }
    s->blk = (const char *) start + blksz;
    if (madvise((void *) start, n, MADV_WILLNEED) != 0) {
        err("madvise() failed");
    }
}


/*  Handles prefetch and merge segment exhaustion.
 */
static int mrg_seg_advance(struct mrg_seg * const s,
                           const size_t blksz,
                           const void * const cur)
{
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);

    if (cur == s->end) {
        void *start = (char *) s->blk - blksz;
        size_t n = (char *) s->end - (char *) start;
        if (n % pagesz != 0) {
            n += pagesz - n % pagesz;
        }
        if (madvise(start, n, MADV_DONTNEED) != 0) {
            err("madvise() failed");
        }
        return 0;
    }
    assert(cur >= s->blk && cur < s->end);
    if (madvise((char *) s->blk - blksz, blksz, MADV_DONTNEED) != 0) {
        err("madvise() failed");
    }
    s->cur = cur;
    s->blk = (const char *) s->blk + blksz;
    if (s->blk < s->end) {
        size_t n = (char *) s->end - (char *) s->blk;
        if (n >= blksz) {
            n = blksz;
        } else if (n % pagesz != 0) {
            n += pagesz - n % pagesz;
        }
        if (madvise((void *) s->blk, n, MADV_WILLNEED) != 0) {
            err("madvise() failed");
        }
    }
    return 1;
}


/*  Consumes one item in the given merge segment; returns 0 if there
    are no more items in the segment.
 */
HTM_INLINE int mrg_seg_consume(struct mrg_seg * const s,
                               const size_t blksz,
                               const size_t itemsz)
{
    const void *cur = (const char *) s->cur + itemsz;
    if (cur < s->end && cur < s->blk) {
        s->cur = cur;
        return 1;
    }
    return mrg_seg_advance(s, blksz, cur);
}


/*  Adds segs[n] to the min-heap segs[0], segs[1], ..., segs[n - 1].
 */
static void heap_up(struct mrg_seg *segs,
                    size_t n,
                    int (*cmpfn)(const void *, const void *))
{
    struct mrg_seg tmp;
    size_t p;

    while (n != 0) {
        p = (n - 1) / 2;
        if ((*cmpfn)(segs[p].cur, segs[n].cur) != 0) {
            break;
        }
        tmp = segs[p];
        segs[p] = segs[n];
        segs[n] = tmp;
        n = p;
    }
}


/*  Fix possible min-heap violation by segs[0].
 */
static void heap_down(struct mrg_seg *segs,
                      const size_t n,
                      int (*cmpfn)(const void *, const void *))
{
    struct mrg_seg tmp;
    size_t i;

    if (n > 1) {
        i = 0;
        while (1) {
            size_t left = 2*i + 1;
            size_t right = 2*i + 2;
            size_t least = i;
            if (left < n && (*cmpfn)(segs[left].cur, segs[i].cur) != 0) {
                least = left;
            }
            if (right < n && (*cmpfn)(segs[right].cur, segs[least].cur) != 0) {
                least = right;
            }
            if (least == i) {
                break;
            }
            tmp = segs[i];
            segs[i] = segs[least];
            segs[least] = tmp;
            i = least;
        }
    }
}


/*  Performs one multi-way merge pass.
 */
static void mrg_pass(struct blk_writer * const w,
                     const void * const data,
                     const struct mem_params * const mem,
                     const size_t filesz,
                     const size_t sortsz,
                     const size_t itemsz,
                     int (*cmpfn)(const void *, const void *))
{
    struct mrg_seg *segs;
    size_t start;

    /* allocate merge segments */
    segs = (struct mrg_seg *) malloc(
        mem->k * sizeof(struct mrg_seg));
    if (segs == NULL) {
        err("malloc() failed");
    }
    for (start = 0; start < filesz;) {
        size_t ns, end;
        /* initialize up to mem->k merge segments */
        for (ns = 0; ns < mem->k && start < filesz; ++ns, start = end) {
            end = (start + sortsz > filesz ? filesz : start + sortsz);
            mrg_seg_init(&segs[ns],
                         (const char *) data + start,
                         (const char *) data + end,
                         mem->ioblksz);
            heap_up(segs, ns, cmpfn);
        }
        /* merge ns segments */
        while (ns > 0) {
            /* write minimum value from all ns merge segments to disk */
            blk_writer_append(w, segs->cur, itemsz, NULL);
            if (mrg_seg_consume(segs, mem->ioblksz, itemsz) == 0) {
                segs[0] = segs[ns - 1];
                --ns;
            }
            heap_down(segs, ns, cmpfn);
        }
    }
    free(segs);
}


/*  Computes the number of k-way merge passes required to sort n items,
    i.e. the ceiling of the base-k logarithm of n.
 */
static int mrg_npasses(size_t n, size_t k) {
    int m = 1;
    size_t b = k;
    while (b < n) {
        b *= k;
        ++m;
    }
    return m;
}


/*  Sorts the given file using multi-way external merge sort.

    @param[in] file     File to sort. Runs of size mem->sortb are assumed
                        to be sorted already.
    @param[in] tmpl     Temporary file name template.
    @param[in] mem      Memory parameters.
    @param[in] itemsz   Size of a single item in bytes.
    @param[in] nitems   Number of items in file.
 */
static void ext_sort(const char *file,
                     const char *scratch,
                     const struct mem_params *mem,
                     size_t itemsz,
                     size_t nitems,
                     int (*cmpfn)(const void *, const void *))
{
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);
    const size_t filesz = nitems * itemsz;
    size_t sortsz = mem->sortsz - mem->sortsz % itemsz;
    const size_t nblk = filesz / sortsz + (filesz % sortsz != 0 ? 1u : 0u);
    const int nmp = mrg_npasses(nblk, mem->k);
    int mp;
    double ttot;

    if (nblk < 2) {
        msg("Skipping multi-way merge step (%s already sorted)\n", file);
        return;
    }
    msg("Multi-way external merge sort of %s\n", file);
    ttot = now();

    /* Perform k-way merge passes; the sorted block grows by a
       factor of k on every pass */
    for (mp = 0; mp < nmp; ++mp, sortsz *= mem->k) {
        const char *inf;
        const char *outf;
        const void *data;
        struct blk_writer *w;
        double t;
        size_t nmap;
        int infd;

        msg("\t- merge pass %d/%d ... ", mp + 1, nmp);
        t = now();

        /* memory map input file and create writer for output file */
        if (mp % 2 == 0) {
            inf = file;
            outf = scratch;
        } else {
            inf = scratch;
            outf = file;
        }
        infd = open(inf, O_RDONLY);
        if (infd == -1) {
            err("failed to open file %s for reading", inf);
        }
        if (filesz % pagesz != 0) {
            nmap = filesz + (pagesz - filesz % pagesz);
        } else {
            nmap = filesz;
        }
        data = mmap(NULL, nmap, PROT_READ,
                    MAP_SHARED | MAP_NORESERVE, infd, 0);
        if (data == MAP_FAILED) {
            err("failed to mmap() file %s", inf);
        }
        if (madvise((void *) data, filesz, MADV_DONTNEED) != 0) {
            err("madvise() failed");
        }
        w = blk_writer_init(outf, mem->ioblksz);
        /* perform merges */
        mrg_pass(w, data, mem, filesz, sortsz, itemsz, cmpfn);
        /* cleanup */
        blk_writer_close(w, NULL);
        if (munmap((void *) data, nmap) != 0) {
            err("failed to munmap() file %s", inf);
        }
        if (close(infd) != 0) {
            err("failed to close() file %s", inf);
        }
        msg("%.3f sec\n", now() - t);
    }
    /* make sure sorted results are in data file and delete scratch file. */
    if (nmp % 2 == 1) {
        if (rename(scratch, file) != 0) {
            err("failed to rename file %s to %s", scratch, file);
        }
    } else {
        if (unlink(scratch) != 0) {
            err("failed to delete file %s", scratch);
        }
    }
    msg("\t%.3f sec total\n\n", now() - ttot);
}


/* ---- Fast constant size memory allocator ---- */

/*  Undefine FAST_ALLOC (or define it as 0) to allocate nodes with malloc()
    instead. Useful when checking memory safety, e.g. with valgrind.
 */
#define FAST_ALLOC 1

#if FAST_ALLOC
#   define ARENA_SEGSZ 16*1024*1024

/*  A contiguous memory segment belonging to an arena.
 */
struct arena_seg {
    struct arena_seg *prev;
    void *mem;
};

static struct arena_seg * arena_seg_init(struct arena_seg * const prev,
                                         const size_t itemsz)
{
    unsigned char *node, *end;
    struct arena_seg *seg = (struct arena_seg *) malloc(sizeof(struct arena_seg));
    if (seg == NULL) {
        err("malloc() failed");
    }
    seg->prev = prev;
    seg->mem = mmap(NULL, ARENA_SEGSZ, PROT_READ | PROT_WRITE,
                    MAP_ANON | MAP_PRIVATE, -1, 0);
    if (seg->mem == MAP_FAILED) {
        err("mmap() failed");
    }
    /* initialize free list */
    node = (unsigned char *) seg->mem;
    end = node + (ARENA_SEGSZ/itemsz - 1)*itemsz;
    while (node < end) {
        *((void **) node) = node + itemsz;
        node += itemsz;
    }
    *((void **) node) = NULL;
    return seg;
}

/*  A (non-shrinkable) memory arena. Memory nodes can be freed one
    at a time, or en-masse.
 */
struct arena {
    struct arena_seg *tail; /* reverse linked list of memory segments */
    struct mem_node *head;  /* head of linked list of free memory locations */
    size_t itemsz;
    size_t nseg;
};

static void arena_init(struct arena * const a, const size_t itemsz)
{
    a->tail = arena_seg_init(NULL, itemsz);
    a->head = a->tail->mem;
    a->itemsz = itemsz;
    a->nseg = 1;
}

static void arena_destroy(struct arena * const a)
{
    struct arena_seg *seg = a->tail;
    while (seg != NULL) {
        struct arena_seg *prev = seg->prev;
        if (munmap(seg->mem, ARENA_SEGSZ) != 0) {
            err("munmap() failed");
        }
        seg->prev = NULL;
        seg->mem = NULL;
        free(seg);
        seg = prev;
    }
    a->tail = NULL;
    a->head = NULL;
    a->nseg = 0;
}

HTM_INLINE void * arena_alloc(struct arena * const a)
{
    void *item;
    if (a->head == NULL) {
        a->tail = arena_seg_init(a->tail, a->itemsz);
        a->head = a->tail->mem;
        ++a->nseg;
    }
    item = a->head;
    a->head = *((void **) item);
    return item;
}

HTM_INLINE void arena_free(struct arena * const a, void * const n)
{
    *((void **) n) = a->head;
    a->head = n;
}

#endif /* FAST_ALLOC */


/* ================================================================ */
/*           Phase 1: Produce data file from ASCII inputs           */
/* ================================================================ */

/*  An entry in an HTM tree.
 */
struct tree_entry {
    int64_t htmid;
    int64_t rowid;
    struct htm_sc sc;
} HTM_ALIGNED(16);

/*  Tests whether tree entry e1 is less than e2; this is
    the same as testing whether the HTM ID of e1 is less than e2.
    Row IDs are used to break ties.
  */
HTM_INLINE int tree_entry_lt(const struct tree_entry *e1,
                             const struct tree_entry *e2)
{
    return (e1->htmid < e2->htmid ||
            (e1->htmid == e2->htmid && e1->rowid < e2->rowid));
}

/*  Returns 1 if the first tree entry is less than the second.
 */
static int tree_entry_cmp(const void *e1, const void *e2)
{
    return tree_entry_lt((const struct tree_entry *) e1,
                         (const struct tree_entry *) e2);
}

/* ---- Quicksort of tree entries ---- */

static void tree_entry_isort(struct tree_entry *entries, size_t n)
{
    size_t i, j, k;
    for (i = 0; i < n; ++i) {
        k = i;
        for (j = i + 1; j < n; ++j) {
            if (tree_entry_lt(&entries[j], &entries[k])) {
                k = j;
            }
        }
        if (k != i) {
            struct tree_entry tmp = entries[k];
            entries[k] = entries[i];
            entries[i] = tmp;
        }
    }
}

static void tree_entry_qsort(struct tree_entry *entries,
                             size_t left,
                             size_t right)
{
    struct tree_entry pivot;
    size_t l, r, mid;

    while (1) {

        if (right <= left + 8) {
            tree_entry_isort(entries + left, right - left + 1);
            return;
        }

        /* find median-of-3 */
        mid = left + ((right - left) >> 1);
        if (tree_entry_lt(entries + left, entries + mid)) {
            if (tree_entry_lt(entries + right, entries + left)) {
                mid = left; /* right, left, mid */
            } else if (tree_entry_lt(entries + right, entries + mid)) {
                mid = right; /* left, right, mid */
            } /* left, mid, right */
        } else {
            if (tree_entry_lt(entries + left, entries + right)) {
                mid = left; /* mid, left, right */
            } else if (tree_entry_lt(entries + mid, entries + right)) {
                mid = right; /* mid, right, left */
            } /* right, mid, left */
        }

        /* set pivot to median-of-3 and store at left-most array location */
        pivot = entries[mid];
        if (mid != left) {
            entries[mid] = entries[left];
            entries[left] = pivot;
        }

        /* partition around pivot */
        l = left;
        r = right;
        do {
            while (tree_entry_lt(&pivot, &entries[right])) {
                --right;
                if (left >= right) {
                    goto end;
                }
            }
            entries[left] = entries[right];
            do {
                ++left;
                if (left >= right) {
                    left = right;
                    goto end;
                }
            } while (tree_entry_lt(&entries[left], &pivot));
            entries[right] = entries[left];
            --right;
        } while (left < right);
end:
        entries[left] = pivot;

        /* recurse on smaller partition */
        if (2*left <= r + l) {
            if (l < left) {
                tree_entry_qsort(entries, l, left - 1);
            }
            ++left;
            right = r;
        } else {
            if (r > left) {
                tree_entry_qsort(entries, left + 1, r);
            }
            right = left - 1;
            left = l;
        }
        /* process larger partition without recursing to limit stack usage */
    }
}


/*  Sorts a memory block containing tree entries.
 */
static void tree_entry_sort(void *data, size_t nbytes)
{
    size_t n;
    assert(((size_t) data) % 16 == 0 &&
           "block pointer not aligned to a multiple of 16 bytes");
    assert(nbytes % sizeof(struct tree_entry) == 0 &&
           "block size not a multiple of sizeof(struct tree_entry)");
    n = nbytes / sizeof(struct tree_entry);
    if (n > 1) {
        tree_entry_qsort((struct tree_entry *) data, 0, n - 1);
    }
}


/* ---- Convert text inputs to a block-sorted binary tree entry file ---- */

static char * eat_delim(char *s, char delim, const char *fname, size_t lineno)
{
    for (; isspace(*s) && *s != delim; ++s) { }
    if (*s != delim) {
        err("[%s:%llu] - invalid/truncated record",
            fname, (unsigned long long) lineno);
    }
    return s + 1;
}

static char * eat_ws(char *s, char delim, const char *fname, size_t lineno)
{
    for (; isspace(*s) && *s != delim; ++s) { }
    if (*s == delim || *s == '\0') {
        err("[%s:%llu] - invalid/truncated record",
            fname, (unsigned long long) lineno);
    }
    return s;
}


/*  Converts ASCII input files to a block-sorted binary file.
 */
static size_t blk_sort_ascii(char **infile,
                             const int nfile,
                             const char *outfile,
                             const char delim,
                             const struct mem_params * const mem)
{
    char line[16384];
    size_t nentries;
    struct blk_writer *out;
    double ttot;
    int i;

    if (mem->sortsz % sizeof(struct tree_entry) != 0) {
        err("Write block size is not a multiple of "
            "sizeof(struct tree_entry)");
    }
    msg("Creating block-sorted tree entry file %s from ASCII file(s)\n",
        outfile);
    ttot = now();
    nentries = 0;
    out = blk_writer_init(outfile, mem->sortsz);

    /* For each input file... */
    for (i = 0; i < nfile; ++i) {
        FILE *f;
        size_t lineno;
        struct tree_entry entry;
        struct htm_v3 v;
        double t;

        msg("\t- processing %s ... ", infile[i]);
        t = now();

        /* open input file */
        f = fopen(infile[i], "r");
        if (f == NULL) {
            err("Failed to open file %s for reading", infile[i]);
        }
        for (lineno = 1; fgets(line, sizeof(line), f) != NULL; ++lineno) {
            char *s, *endptr;
            double lon, lat;
            size_t len = strlen(line);

            if (line[len - 1] != '\n') {
                if (feof(f) == 0) {
                    err("Line %llu of file %s is too long (> %d characters)",
                        (unsigned long long) lineno, infile[i],
                        (int) sizeof(line));
                }
            }
            s = eat_ws(line, delim, infile[i], lineno);
            entry.rowid = (int64_t) strtoll(s, &endptr, 0);
            if (endptr == s || endptr == NULL || errno != 0) {
                err("[%s:%llu] - failed to convert row_id to an integer",
                    infile[i], (unsigned long long) lineno);
            }
            s = eat_delim(endptr, delim, infile[i], lineno);
            s = eat_ws(s, delim, infile[i], lineno);
            lon = strtod(s, &endptr);
            if (endptr == s || endptr == NULL || errno != 0) {
                err("[%s:%llu] - failed to convert right ascension/longitude "
                    "to a double", infile[i], (unsigned long long) lineno);
            }
            s = eat_delim(endptr, delim, infile[i], lineno);
            s = eat_ws(s, delim, infile[i], lineno);
            lat = strtod(s, &endptr);
            if (endptr == s || endptr == NULL || errno != 0) {
                err("[%s:%llu] - failed to convert declination/latitude "
                    "to a double", infile[i], (unsigned long long) lineno);
            }
            s = endptr;
            if (*s != delim && *s != '\0' && !isspace(*s)) {
                err("[%s:%llu] - invalid record",
                    infile[i], (unsigned long long) lineno);
            }
            /* compute and store tree_entry for line */
            if (htm_sc_init(&entry.sc, lon, lat) != HTM_OK) {
                err("[%s:%llu] - invalid spherical coordinates",
                    infile[i], (unsigned long long) lineno);
            }
            if (htm_sc_tov3(&v, &entry.sc) != HTM_OK) {
                err("[%s:%llu] - failed to convert spherical coordinates "
                    "to a unit vector", infile[i], (unsigned long long) lineno);
            }
            entry.htmid = htm_v3_id(&v, 20);
            if (entry.htmid == 0) {
                err("[%s:%llu] - failed to compute HTM ID for spherical "
                    "coordinates", infile[i], (unsigned long long) lineno);
            }
            blk_writer_append(out, &entry, sizeof(struct tree_entry),
                              &tree_entry_sort);
        }
        if (ferror(f) != 0) {
            err("failed to read file %s", infile[i]);
        }
        if (fclose(f) != 0) {
            err("failed to close file %s", infile[i]);
        }
        nentries += lineno - 1;
        msg("%llu records in %.3f sec\n",
            (unsigned long long) lineno, now() - t);
        /* advance to next input file */
    }

    /* flush and close block writer */
    blk_writer_close(out, &tree_entry_sort);
    msg("\t%.3f sec for %llu records total\n\n",
        now() - ttot, (unsigned long long) nentries);
    return nentries;
}


/* ================================================================ */
/*               Phase 2: Tree generation and layout                */
/* ================================================================ */

/*  Number of levels-of-detail used by Split-and-Refine.
 */
#define NLOD 5

/*  Tree layout block sizes in bytes, from largest to smallest.
 */
static const uint32_t layout_size[NLOD] = {
    2097152,    /* 2MiB: large page size for modern x86 processors. */
    65536,      /* Between small and large page size. Chosen in hopes
                   of improving OS disk prefetch effectiveness. */
    4096,       /* 4KiB: default page size for modern x86 processors */
    256,        /* Between cache line and small page sizes. Chosen in hopes
                   of improving HW cache-line prefetch effectiveness. */
    64          /* L1/L2/L3 cache-line size for modern x86 processors. */
};

/*  A node ID, consisting of NLOD block IDs (one per layout block size)
    and the index of the node in a post-order traversal. This last index
    makes node IDs unique (more than one node might fit in a block at
    the finest LOD).
 */
struct node_id {
    uint64_t block[NLOD + 1];
} HTM_ALIGNED(16);

/*  On-disk representation of a tree node.
 */
struct disk_node {
    struct node_id id;
    uint64_t count;
    uint64_t index;
    struct node_id child[4];
} HTM_ALIGNED(16);

/*  Returns 1 if the given node ID corresponds to an empty child
    (all block IDs are zero).
 */
HTM_INLINE int node_empty(const struct node_id * const id)
{
    int i;
    for (i = 0; i < NLOD + 1; ++i) {
        if (id->block[i] != 0) {
            return 0;
        }
    }
    return 1;
}

/*  Returns 1 if the given node IDs are equal.
 */
HTM_INLINE int node_id_eq(const struct node_id * const id1,
                          const struct node_id * const id2)
{
    int i;
    for (i = 0; i < NLOD + 1; ++i) {
        if (id1->block[i] != id2->block[i]) {
            return 0;
        }
    }
    return 1;
}

/*  Returns 1 if the first node ID is less than the second. We say that
    a node N1 is less than N2 if the string of block IDs for N1 is
    lexicographically less than the string for N2.
 */
HTM_INLINE int node_id_lt(const struct node_id * const id1,
                          const struct node_id * const id2)
{
    int i;
    for (i = 0; i < NLOD + 1; ++i) {
        if (id1->block[i] < id2->block[i]) {
            return 1;
        } else if (id1->block[i] > id2->block[i]) {
            break;
        }
    }
    return 0;
}

/*  Compares nodes by ID.
 */
HTM_INLINE int disk_node_lt(const struct disk_node *n1,
                            const struct disk_node *n2)
{
    return node_id_lt(&n1->id, &n2->id);
}

static int disk_node_cmp(const void *n1, const void * n2)
{
    return disk_node_lt((const struct disk_node *) n1,
                        (const struct disk_node *) n2);
}


/* ---- Quicksort of disk nodes ---- */

static void disk_node_isort(struct disk_node *nodes, size_t n)
{
    size_t i, j, k;
    for (i = 0; i < n; ++i) {
        k = i;
        for (j = i + 1; j < n; ++j) {
            if (disk_node_lt(&nodes[j], &nodes[k])) {
                k = j;
            }
        }
        if (k != i) {
            struct disk_node tmp = nodes[k];
            nodes[k] = nodes[i];
            nodes[i] = tmp;
        }
    }
}

static void disk_node_qsort(struct disk_node *nodes, size_t left, size_t right)
{
    struct disk_node pivot;
    size_t l, r, mid;

    while (1) {

        if (right <= left + 8) {
            disk_node_isort(nodes + left, right - left + 1);
            return;
        }

        /* find median-of-3 */
        mid = left + ((right - left) >> 1);
        if (disk_node_lt(nodes + left, nodes + mid)) {
            if (disk_node_lt(nodes + right, nodes + left)) {
                mid = left; /* right, left, mid */
            } else if (disk_node_lt(nodes + right, nodes + mid)) {
                mid = right; /* left, right, mid */
            } /* left, mid, right */
        } else {
            if (disk_node_lt(nodes + left, nodes + right)) {
                mid = left; /* mid, left, right */
            } else if (disk_node_lt(nodes + mid, nodes + right)) {
                mid = right; /* mid, right, left */
            } /* right, mid, left */
        }

        /* set pivot to median-of-3 and store at left-most array location */
        pivot = nodes[mid];
        if (mid != left) {
            nodes[mid] = nodes[left];
            nodes[left] = pivot;
        }

        /* partition around pivot */
        l = left;
        r = right;
        do {
            while (disk_node_lt(&pivot, &nodes[right])) {
                --right;
                if (left >= right) {
                    goto end;
                }
            }
            nodes[left] = nodes[right];
            do {
                ++left;
                if (left >= right) {
                    left = right;
                    goto end;
                }
            } while (disk_node_lt(&nodes[left], &pivot));
            nodes[right] = nodes[left];
            --right;
        } while (left < right);
end:
        nodes[left] = pivot;

        /* recurse on smaller partition */
        if (2*left <= r + l) {
            if (l < left) {
                disk_node_qsort(nodes, l, left - 1);
            }
            ++left;
            right = r;
        } else {
            if (r > left) {
                disk_node_qsort(nodes, left + 1, r);
            }
            right = left - 1;
            left = l;
        }
        /* process larger partition without recursing to limit stack usage */
    }
}


/*  Sorts a memory block containing disk nodes.
 */
static void disk_node_sort(void *data, size_t nbytes)
{
    size_t n;
    assert(((size_t) data) % 16 == 0 &&
           "block pointer not aligned to a multiple of 16 bytes");
    assert(nbytes % sizeof(struct disk_node) == 0 &&
           "block size not a multiple of sizeof(struct disk_node)");
    n = nbytes / sizeof(struct disk_node);
    if (n > 1) {
        disk_node_qsort((struct disk_node *) data, 0, n - 1);
    }
}


/* ---- In-memory node representation ---- */

/*  Status of an in-memory tree node.
 */
enum node_status {
    NODE_INIT = 0,  /* node is minty fresh */
    NODE_EMITTED,   /* node has been processed by emit_node() */
    NODE_LAID_OUT,  /* node has been processed by layout_node() */
};

/*  In-memory representation of a tree node.  Note that block size/depth
    is packed into a single 32 bit integer; as a result, nodes occupy
    exactly 256 bytes and are nicely cache aligned.
 */
struct mem_node {
    struct node_id id;        /* Hierarchical ID for node */
    int64_t htmid;            /* HTM ID of node */
    uint64_t index;           /* File offset of first point in node */
    uint64_t count;           /* Number of points in node */
    enum node_status status;
    uint32_t blockinfo[NLOD]; /* Clark & Munro: block depth (8 MSBs) and
                                 block size (24 LSBs) for each LOD. */
    struct mem_node *child[4];
} HTM_ALIGNED(16);


/* get/set block size/depth from 32 bit blockinfo */
HTM_INLINE uint32_t get_block_size(uint32_t blockinfo) {
    return blockinfo & 0xffffff;
}
HTM_INLINE uint8_t get_block_depth(uint32_t blockinfo) {
    return blockinfo >> 24;
}
HTM_INLINE uint32_t make_block_info(uint32_t size, uint8_t depth) {
    return (size & 0xffffff) | (((uint32_t) depth) << 24);
}


/* ---- Tree generation ---- */

/*  Tree generation context.
 */
struct tree_gen_context {
#if FAST_ALLOC
    struct arena ar;        /* node memory allocator */
#endif
    size_t   nnodes;        /* number of nodes in the tree */
    uint64_t leafthresh;    /* maximum # of points per leaf */
    uint64_t poidx;         /* next post-order tree traversal index */
    uint64_t blockid[NLOD]; /* index of next block ID to assign for each LOD */
    struct blk_writer *wr;  /* node writer */
};

static void tree_gen_context_init(struct tree_gen_context * const ctx,
                                  const uint64_t leafthresh,
                                  const char * const file,
                                  const size_t blksz)
{
    int i;
    ctx->nnodes = 0;
    ctx->leafthresh = leafthresh;
    ctx->poidx = 0;
    for (i = 0; i < NLOD; ++i) {
        ctx->blockid[i] = 0;
    }
#if FAST_ALLOC
    arena_init(&ctx->ar, sizeof(struct mem_node));
#endif
    ctx->wr = blk_writer_init(file, blksz);
}

static void tree_gen_context_destroy(struct tree_gen_context * const ctx)
{
#if FAST_ALLOC
    arena_destroy(&ctx->ar);
#endif
    blk_writer_close(ctx->wr, &disk_node_sort);
    ctx->wr = NULL;
}


/*  Assigns a block ID at the specified level-of-detail to all nodes in
    the sub-tree rooted at n that do not already have a block ID at that
    level-of-detail. When assigning a block ID to a node, check whether all
    block IDs are now valid and if so write the node out to disk.
    Children of nodes that are written are destroyed.
 */
static void assign_block(struct tree_gen_context * const ctx,
                         struct mem_node * const n,
                         const uint64_t blockid,
                         const int lod)
{
    int i;
    if (n->id.block[lod] != 0) {
        return;
    }
    /* visit children */
    for (i = 0; i < 4; ++i) {
        if (n->child[i] != NULL) {
            assign_block(ctx, n->child[i], blockid, lod);
        }
    }
    n->id.block[lod] = blockid;
    for (i = 0; i < NLOD; ++i) {
        if (n->id.block[i] == 0) {
            return;
        }
    }
    /* write node to disk */
    {
        struct disk_node d;
        d.id = n->id;
        d.count = n->count;
        d.index = n->index;
        for (i = 0; i < 4; ++i) {
            struct mem_node *tmp = n->child[i];
            if (tmp != NULL) {
                /* copy id of child to disk node, throw away child */
                d.child[i] = tmp->id;
#if FAST_ALLOC
                arena_free(&ctx->ar, tmp);
#else
                free(tmp);
#endif
                n->child[i] = NULL;
            } else {
                memset(&d.child[i], 0, sizeof(struct node_id));
            }
        }
        blk_writer_append(ctx->wr, &d, sizeof(struct disk_node), &disk_node_sort);
        ++ctx->nnodes;
    }
}


/*  Estimates the compressed on-disk size of a tree node.
 */
static uint32_t estimate_node_size(const struct mem_node * const node,
                                   const uint32_t nchild)
{
    uint32_t sz = htm_varint_len(node->index) + htm_varint_len(node->count);
    if (nchild > 0) {
        /* There is no way to compute size of a child offset accurately
           without knowing the final node layout, so use a guess of 4 bytes
           per non-empty child. The offset for an empty child (0) will
           will occupy exactly 1 byte, regardless of layout. */
        sz += nchild*3 + 4;
    }
    return sz;
}

struct child_info {
    struct mem_node *node;
    uint32_t size;
    uint8_t depth;
    int8_t idx;
};

/*  Insertion sort of (at most 8) child_info structs; they are sorted
    by block depth and size.
 */
static void child_info_isort(struct child_info * const c, const int n)
{
    if (n < 2) {
        return;
    } else {
        int i, j, k;
        for (i = 0; i < n; ++i) {
            k = i;
            for (j = i + 1; j < n; ++j) {
                if (c[j].depth < c[k].depth ||
                    (c[j].depth == c[k].depth && c[j].size < c[k].size)) {
                    k = j;
                }
            }
            if (k != i) {
                struct child_info tmp = c[k];
                c[k] = c[i];
                c[i] = tmp;
            }
        }
    }
}

static void layout_node(struct mem_node * const node,
                        struct tree_gen_context * const ctx)
{
    struct child_info cinfo[4];
    int c, nchild;

    if (node->status > NODE_EMITTED) {
        return;
    }

    /* visit children */
    for (c = 0, nchild = 0; c < 4; ++c) {
        if (node->child[c] != NULL) {
            layout_node(node->child[c], ctx);
            cinfo[nchild].node = node->child[c];
            cinfo[nchild].idx = c;
            ++nchild;
        }
    }

    /* update status and assign post-order index. */
    node->status = NODE_LAID_OUT;
    node->id.block[5] = ++ctx->poidx;

    /* Clark & Munro at every level-of-detail */
    if (nchild == 0) {
        /* leaf */
        int lod;
        const uint32_t nodesz = estimate_node_size(node, nchild);
        const uint32_t info = make_block_info(nodesz, 1u);
        for (lod = 0; lod < NLOD; ++lod) {
            node->blockinfo[lod] = info;
            if (nodesz > layout_size[lod]) {
                uint64_t blockid = ++ctx->blockid[lod];
                assign_block(ctx, node, blockid, lod);
            }
        }
    } else {
        /* internal node */
        int lod;
        const uint32_t nodesz = estimate_node_size(node, nchild);
        for (lod = 0; lod < NLOD; ++lod) {
            uint64_t blockid;
            uint32_t totsz = nodesz;
            int close = 0, endclose = nchild;

            for (c = 0; c < nchild; ++c) {
                uint32_t s;
                struct mem_node *tmp = cinfo[c].node;
                s = get_block_size(tmp->blockinfo[lod]);
                cinfo[c].size = s;
                cinfo[c].depth = get_block_depth(tmp->blockinfo[lod]);
                totsz += s;
            }
            child_info_isort(cinfo, nchild);

            if (cinfo[0].depth == cinfo[nchild - 1].depth) {
                /* all children have the same block depth */
                if (totsz <= layout_size[lod]) {
                    /* children and parent all fit in one block; merge them */
                    node->blockinfo[lod] = make_block_info(totsz, cinfo[0].depth);
                    continue;
                } else {
                    /* cannot fit family into one block: scan children
                       from smallest to largest, placing as many as possible
                       in the same block as the parent. */
                    totsz = nodesz;
                    for (close = 0; close < nchild - 1; ++close) {
                        if (totsz + cinfo[close].size > layout_size[lod]) {
                            break;
                        }
                        totsz += cinfo[close].size;
                    }
                    /* increase block depth of parent by 1 */
                    node->blockinfo[lod] = make_block_info(totsz, cinfo[0].depth + 1);
                }
            } else {
                /* nchild > 1, not all children have the same block depth */
                totsz = nodesz;
                for (endclose = nchild - 1; endclose > 0; --endclose) {
                    totsz += cinfo[endclose].size;
                    if (cinfo[endclose - 1].depth != cinfo[nchild - 1].depth) {
                        break;
                    }
                }
                if (totsz < layout_size[lod]) {
                    /* merge parent and largest-depth children */
                    node->blockinfo[lod] = make_block_info(totsz, cinfo[nchild - 1].depth);
                } else {
                    /* fresh block for parent, increase block depth by 1 */
                    node->blockinfo[lod] = make_block_info(nodesz, cinfo[nchild - 1].depth + 1);
                    endclose = nchild;
                }
            }
            /* scan remaining children from smallest to largest, merging
               runs of children into a single block where possible. */
            totsz = cinfo[close].size;
            for (c = close + 1; c < endclose; ++c) {
                if (totsz + cinfo[c].size > layout_size[lod]) {
                    blockid = ++ctx->blockid[lod];
                    for (; close < c; ++close) {
                        assign_block(ctx, cinfo[close].node, blockid, lod);
                    }
                    totsz = cinfo[c].size;
                } else {
                    totsz += cinfo[c].size;
                }
            }
            blockid = ++ctx->blockid[lod];
            for (; close < endclose; ++close) {
                assign_block(ctx, cinfo[close].node, blockid, lod);
            }
        }
    }
}

/*  Called on a node when all points belonging to it have been accounted for.
 */
static void emit_node(struct mem_node * const node,
                      struct tree_gen_context * const ctx)
{
    int c;
    if (node->status > NODE_INIT) {
        return;
    }
    /* visit children */
    for (c = 0; c < 4; ++c) {
        if (node->child[c] != NULL) {
            emit_node(node->child[c], ctx);
        }
    }
    if (node->count < ctx->leafthresh) {
        /* if the node point count is too small,
           make it a leaf by deleting all children. */
        for (c = 0; c < 4; ++c) {
            if (node->child[c] != NULL) {
#if FAST_ALLOC
                arena_free(&ctx->ar, node->child[c]);
#else
                free(node->child[c]);
#endif
                node->child[c] = NULL;
            }
        }
        node->status = NODE_EMITTED;
    } else {
        /* otherwise, layout the subtree rooted at node */
        layout_node(node, ctx);
    }
}

/*  Adds a new node to root (one of S0,S1,S2,S3,N0,N1,N2,N3).
 */
static void add_node(struct mem_node * const root,
                     struct tree_gen_context * const ctx,
                     int64_t htmid,
                     int64_t count,
                     int64_t index)
{
    struct mem_node *node;
    int lvl = 0;

    for (lvl = 0, node = root; lvl < 20; ++lvl) {
        /* keep subdividing */
        int i, c;
        node->count += count;
        c = (htmid >> 2*(19 - lvl)) & 3;
        for (i = 0; i < c; ++i) {
            struct mem_node *tmp = node->child[i];
            if (tmp != NULL && tmp->status == NODE_INIT) {
                emit_node(node->child[i], ctx);
            }
        }
        index -= node->index; /* relativize index */
        if (node->child[c] != NULL) {
            node = node->child[c];
        } else {
            /* create child node */
#if FAST_ALLOC
            struct mem_node *child = (struct mem_node *) arena_alloc(&ctx->ar);
#else
            struct mem_node *child =
                (struct mem_node *) malloc(sizeof(struct mem_node));
            if (child == NULL) {
                err("malloc() failed");
            }
#endif
            memset(child, 0, sizeof(struct mem_node));
            child->index = index;
            child->htmid = node->htmid*4 + c;
            node->child[c] = child;
            node = child;
        }
    }
    assert(node->htmid == htmid);
    node->count = count;
}


/*  Container for the 8 level-0 HTM tree nodes.
 */
struct tree_root {
    uint64_t count; /* Total number of points in tree */
    struct mem_node *child[8];
    struct node_id childid[8];
};


/*  Assigns block IDs to the HTM root nodes.
 */
static void finish_root(struct tree_root * const super,
                        struct tree_gen_context * const ctx)
{
    struct child_info cinfo[8];
    int c, close, nchild, lod;

    for (c = 0, nchild = 0; c < 8; ++c) {
        struct mem_node *tmp = super->child[c];
        if (tmp != NULL) {
            super->count += tmp->count;
            cinfo[nchild].node = tmp;
            cinfo[nchild].idx = c;
            ++nchild;
        }
    }
    /* Clark & Munro for each level of detail */
    for (lod = 0; lod < NLOD; ++lod) {
        uint64_t blockid;
        uint32_t totsz;
        for (c = 0; c < nchild; ++c) {
            struct mem_node *tmp = cinfo[c].node;
            cinfo[c].size = get_block_size(tmp->blockinfo[lod]);
            cinfo[c].depth = get_block_depth(tmp->blockinfo[lod]);
        }
        child_info_isort(cinfo, nchild);
        /* scan children from smallest to largest, merging
           as many as possible into blocks. */
        close = 0;
        totsz = cinfo[0].size;
        for (c = 1; c < nchild; ++c) {
            if (totsz + cinfo[c].size > layout_size[lod]) {
                blockid = ++ctx->blockid[lod];
                for (; close < c; ++close) {
                    assign_block(ctx, cinfo[close].node, blockid, lod);
                }
                totsz = cinfo[c].size;
            }
        }
        blockid = ++ctx->blockid[lod];
        for (; close < nchild; ++close) {
            assign_block(ctx, cinfo[close].node, blockid, lod);
        }
    }
    /* At this point, all nodes are guaranteed to have been written to disk.
       Copy HTM root node IDs to the super root and then throw them away. */
    for (c = 0; c < nchild; ++c) {
        super->childid[cinfo[c].idx] = cinfo[c].node->id;
#if FAST_ALLOC
        arena_free(&ctx->ar, cinfo[c].node);
#else
        free(cinfo[c].node);
#endif
        super->child[cinfo[c].idx] = NULL;
    }
}


/*  Tree generation driver function.
 */
static size_t tree_gen(const char * const datafile,
                       const char * const treefile,
                       const struct mem_params * const mem,
                       struct tree_root * const super,
                       const uint64_t leafthresh,
                       const size_t npoints)
{
    struct tree_gen_context ctx;
    const struct tree_entry *data;
    void *behind;
    int64_t htmid;
    uint64_t index, count;
    size_t i, nseg;
    int r, fd;
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);
    size_t mapsz = npoints * sizeof(struct tree_entry);
    double t;

    if (npoints == 0) {
        err("no input points");
    }
    msg("Generating block sorted tree node file %s from %s\n",
        treefile, datafile);
    t = now();
    fd = open(datafile, O_RDONLY);
    if (fd == -1) {
        err("failed to open file %s for reading", datafile);
    }
    if (mapsz % pagesz != 0) {
        mapsz += pagesz - mapsz % pagesz;
    }
    behind = mmap(NULL, mapsz, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);
    if (behind == MAP_FAILED) {
        err("mmap() file %s for reading", datafile);
    }
    if (madvise(behind, mapsz, MADV_SEQUENTIAL) != 0) {
        err("madvise() failed on mmap for file %s", datafile);
    }
    data = (const struct tree_entry *) behind;
    behind = ((unsigned char *) behind) + mem->ioblksz;
    tree_gen_context_init(&ctx, leafthresh, treefile, mem->sortsz);
    memset(super, 0, sizeof(struct tree_root));

    /* walk over tree entries, adding tree nodes. */
    if (data[0].htmid == 0) {
        err("invalid HTM ID");
    }
    htmid = 0;
    index = 0;
    count = 0;
    r = -1;
    for (i = 0; i < npoints; ++i) {
        if ((void *) &data[i] > behind) {
            void *ptr = ((unsigned char *) behind) - mem->ioblksz;
            if (madvise(ptr, mem->ioblksz, MADV_DONTNEED) != 0) {
                err("madvise() failed");
            }
            behind = ((unsigned char *) behind) + mem->ioblksz;
        }
        if (data[i].htmid == htmid) {
            /* increase point count for the current htmid */
            ++count;
        } else {
            int r2;
            assert(data[i].htmid > htmid);
            if (r >= 0) {
                /* add previous node if there is one */
                add_node(super->child[r], &ctx, htmid, count, (uint64_t) index);
            }
            /* reset index, count, and htmid */
            count = 1;
            index = (uint64_t) i;
            htmid = data[i].htmid;
            r2 = (int) (htmid >> 40) - 8;
            if (r2 < 0 || r2 > 7) {
                err("invalid HTM ID");
            }
            if (r != r2) {
                /* need a new HTM root node */
                if (r >= 0) {
                    /* emit and layout the previous root if there is one */
                    emit_node(super->child[r], &ctx);
                    layout_node(super->child[r], &ctx);
                }
                r = r2;
                /* create new root */
#if FAST_ALLOC
                super->child[r] = (struct mem_node *) arena_alloc(&ctx.ar);
#else
                super->child[r] =
                    (struct mem_node *) malloc(sizeof(struct mem_node));
                if (super->child[r] == NULL) {
                    err("malloc() failed");
                }
#endif
                memset(super->child[r], 0, sizeof(struct mem_node));
                super->child[r]->htmid = r + 8;
                super->child[r]->index = index;
            }
        }
    }
    /* add last node, emit and layout last root */
    add_node(super->child[r], &ctx, htmid, count, (uint64_t) index);
    emit_node(super->child[r], &ctx);
    layout_node(super->child[r], &ctx);

    /* assign block IDs to the HTM roots */
    finish_root(super, &ctx);
    if (super->count != npoints) {
        err("bug in tree generation phase - points not all accounted for");
    }
    /* cleanup */
    i = ctx.nnodes;
    nseg = ctx.ar.nseg;
    tree_gen_context_destroy(&ctx);
    if (munmap((void *) data, mapsz) != 0) {
        err("munmap() failed");
    }
    if (close(fd) != 0) {
        err("close() failed");
    }
    msg("\t%.3f sec total (%llu tree nodes, %llu MiB memory)\n\n",
        now() - t, (unsigned long long) i,
        (unsigned long long) (nseg*ARENA_SEGSZ)/(1024*1024));
    return i;
}


/* ================================================================ */
/*                    Phase 3: Tree compression                     */
/* ================================================================ */

/*  Mapping from a node ID to a relative file offset.
 */
struct id_off {
    struct node_id id;
    uint64_t off;
    struct id_off *next;
} HTM_ALIGNED(16);

/*  Simple hash table that maps node IDs to relative file offsets. The
    implementation uses a power-of-2 sized backing array and chains on
    collision. The hash-code of a node ID is taken to be its post-order
    index (which is unique).
 */
struct hash_table {
    size_t n;
    size_t cap;
    struct id_off **array;
#if FAST_ALLOC
    struct arena ar;
#endif
};


static void hash_table_init(struct hash_table *ht)
{
    ht->cap = 65536; /* must be a power of 2 */
    ht->n = 0;
    ht->array = (struct id_off **) malloc(ht->cap * sizeof(struct id_off *));
    if (ht->array == NULL) {
        err("malloc() failed");
    }
    memset(ht->array, 0, ht->cap * sizeof(struct id_off *));
#if FAST_ALLOC
    arena_init(&ht->ar, sizeof(struct id_off));
#endif
}


static void hash_table_destroy(struct hash_table *ht)
{
#if FAST_ALLOC
    arena_destroy(&ht->ar);
#else
    size_t i;
    struct id_off *ptr;
    for (i = 0; i < ht->cap; ++i) {
        ptr = ht->array[i];
        while (ptr != NULL) {
            struct id_off *next = ptr->next;
            free(ptr);
            ptr = next;
        }
    }
#endif
    free(ht->array);
    ht->array = NULL;
    ht->n = 0;
    ht->cap = 0;
}


/*  Grows hash table capacity by a factor of two.
 */
static void hash_table_grow(struct hash_table *ht) {
    size_t i, cap = ht->cap;
    struct id_off **array = ht->array;

    ht->cap = 2*cap;
    ht->array = (struct id_off **) malloc(ht->cap * sizeof(struct id_off *));
    if (ht->array == NULL) {
        err("malloc() failed");
    }
    memset(ht->array, 0, ht->cap * sizeof(struct id_off *));

    /* add previous hash table entries */
    for (i = 0; i < cap; ++i) {
        struct id_off *e = array[i];
        while (e != NULL) {
            struct id_off *next = e->next;
            size_t j = (2*cap - 1) & (size_t) e->id.block[NLOD];
            e->next = ht->array[j];
            ht->array[j] = e;
            e = next;
        }
    }
}


/*  Adds an id to offset mapping to the given hash table.
 */
static void hash_table_add(struct hash_table * const ht,
                           const struct node_id * const id,
                           uint64_t off)
{
    struct id_off *e;
    size_t i;
    if ((ht->n*3)/4 > ht->cap) {
        hash_table_grow(ht);
    }
    i = (ht->cap - 1) & (size_t) id->block[NLOD];
#if FAST_ALLOC
    e = (struct id_off *) arena_alloc(&ht->ar);
#else
    e = (struct id_off *) malloc(sizeof(struct id_off));
    if (e == NULL) {
        err("malloc() failed");
    }
#endif
    e->id = *id;
    e->off = off;
    e->next = ht->array[i];
    ht->array[i] = e;
    ++ht->n;
}


/*  Returns the offset of the node with the given ID and removes
    the corresponding hash table entry.
 */
static uint64_t hash_table_get(struct hash_table * const ht,
                               const struct node_id * const id)
{
    struct id_off *e, *prev;
    const size_t i = (ht->cap - 1) & (size_t) id->block[NLOD];
    prev = NULL;
    e = ht->array[i];
    /* child must occur before parent */
    while (1) {
        if (e == NULL) {
            err("tree generation bug: parent node written before child");
        }
        if (node_id_eq(id, &e->id)) {
            uint64_t off = e->off;
            if (prev == NULL) {
                ht->array[i] = e->next;
            } else {
                prev->next = e->next;
            }
#if FAST_ALLOC
            arena_free(&ht->ar, e);
#else
            free(e);
#endif
            --ht->n;
            return off;
        }
        prev = e;
        e = e->next;
    }
    /* never reached */
}


/*  Writes out a tree node (byte reversed).
 */
static uint64_t compress_node(struct hash_table * const ht,
                              struct blk_writer * const wr,
                              const struct disk_node * const n,
                              const uint64_t filesz,
                              const uint64_t leafthresh)
{
    unsigned char buf[48];
    unsigned char *s = buf;
    uint64_t sz = filesz;
    unsigned int v;
    int c, leaf;

    /* write out child offsets (from child 3 to child 0) */
    for (c = 3, leaf = 1; c >= 0; --c) {
        if (node_empty(&n->child[c])) {
            *s = 0;
            ++s;
            ++sz;
        } else {
            /* this is tricky - child 3 of n can be laid out immediately after
               n, yielding a child offset of 0. But 0 also means "empty child",
               so instead encode the actual offset + 1. */
            v = htm_varint_rencode(s, sz + 1 - hash_table_get(ht, &n->child[c]));
            s += v;
            sz += v;
            leaf = 0;
        }
    }
    if (leaf != 0) {
        /* n is a leaf: don't store child offsets */
        s -= 4;
        sz -= 4;
    } else if (n->count < leafthresh) {
        err("tree generation bug: internal node contains too few points");
    }
    /* write out relative index, then count */
    v = htm_varint_rencode(s, n->index);
    s += v;
    sz += v;
    v = htm_varint_rencode(s, n->count);
    s += v;
    sz += v;
    /* write out byte reversed node, add node id to hashtable */
    blk_writer_append(wr, buf, (size_t) (s - buf), NULL);
    hash_table_add(ht, &n->id, sz);
    return sz;
}


/*  Writes out a tree header (byte reversed).
 */
static uint64_t write_tree_header(struct hash_table * const ht,
                                  struct blk_writer * const wr,
                                  const struct tree_root * const super,
                                  const uint64_t filesz,
                                  const uint64_t leafthresh)
{
    unsigned char buf[96];
    unsigned char *s;
    uint64_t sz;
    int r;
    unsigned int v;

    /* write offsets of N3, N2, N1, N0, S3, S2, S1, S0 */
    for (r = 7, s = buf, sz = filesz; r >= 0; --r) {
        if (node_empty(&super->childid[r])) {
            *s = 0;
            ++s;
            ++sz;
        } else {
            /* N3 could be laid out immediately after super root,
               yielding a child offset of 0, which means "empty child".
               Therefore encode 1 + actual offset. */
            v = htm_varint_rencode(
                s, sz + 1 - hash_table_get(ht, &super->childid[r]));
            s += v;
            sz += v;
        }
    }
    /* write total number of points in tree */
    v = htm_varint_rencode(s, super->count);
    s += v;
    sz += v;
    /* write leaf threshold */
    v = htm_varint_rencode(s, leafthresh);
    s += v;
    sz += v;
    blk_writer_append(wr, buf, (size_t) (s - buf), NULL);
    if (ht->n != 0) {
        err("tree compression bug: node id hash table non-empty");
    }
    return sz;
}


/*  Tree compression driver function.
 */
static uint64_t tree_compress(const char * const treefile,
                              const char * const scratchfile,
                              const struct mem_params * const mem,
                              const struct tree_root * const super,
                              const size_t nnodes,
                              const uint64_t leafthresh)
{
    struct hash_table ht;
    struct blk_writer *wr;
    const struct disk_node *data;
    void *behind;
    double t;
    size_t i;
    uint64_t filesz;
    int fd;
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);
    size_t mapsz = nnodes * sizeof(struct disk_node);

    if (nnodes == 0) {
        err("no input nodes");
    }
    t = now();
    msg("Generating reversed compressed tree node file %s from %s\n",
        scratchfile, treefile);
    fd = open(treefile, O_RDONLY);
    if (fd == -1) {
        err("failed to open file %s for reading", treefile);
    }
    if (mapsz % pagesz != 0) {
        mapsz += pagesz - mapsz % pagesz;
    }
    behind = mmap(NULL, mapsz, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);
    if (behind == MAP_FAILED) {
        err("mmap() file %s for reading", treefile);
    }
    if (madvise(behind, mapsz, MADV_SEQUENTIAL) != 0) {
        err("madvise() failed on mmap for file %s", treefile);
    }
    data = (const struct disk_node *) behind;
    behind = ((unsigned char *) behind) + mem->ioblksz;
    hash_table_init(&ht);
    wr = blk_writer_init(scratchfile, mem->ioblksz);

    /* write nodes */
    for (i = 0, filesz = 0; i < nnodes; ++i) {
        if (i > 0) {
            assert(disk_node_lt(&data[i - 1], &data[i]) &&
                   "tree node file not sorted");
        }
        if ((void *) &data[i] > behind) {
            void *ptr = ((unsigned char *) behind) - mem->ioblksz;
            if (madvise(ptr, mem->ioblksz, MADV_DONTNEED) != 0) {
                err("madvise() failed");
            }
            behind = ((unsigned char *) behind) + mem->ioblksz;
        }
        filesz = compress_node(&ht, wr, &data[i], filesz, leafthresh);
    }
    /* and tree header */
    filesz = write_tree_header(&ht, wr, super, filesz, leafthresh);
    /* cleanup */
    blk_writer_close(wr, NULL);
    hash_table_destroy(&ht);
    if (munmap((void *) data, mapsz) != 0) {
        err("munmap() failed");
    }
    if (close(fd) != 0) {
        err("close() failed");
    }
    msg("\t%.3f sec total\n\n", now() - t);
    return filesz;
}


/*  Byte reverses an input file to an output file and removes the input file.
 */
static void reverse_file(const char * const infile,
                         const char * const outfile,
                         const struct mem_params * const mem,
                         const uint64_t filesz)
{
    const unsigned char *data;
    struct blk_writer *wr;
    double t;
    uint64_t i, j;
    int fd;
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);
    size_t mapsz = filesz;

    if (filesz == 0) {
        err("cannot reverse an empty file");
    }
    msg("Reversing file %s to produce %s\n", infile, outfile);
    t = now();
    fd = open(infile, O_RDONLY);
    if (fd == -1) {
        err("failed to open file %s for reading", infile);
    }
    if (mapsz % pagesz != 0) {
        mapsz += pagesz - mapsz % pagesz;
    }
    data = (const unsigned char *) mmap(NULL, mapsz, PROT_READ,
                                        MAP_SHARED | MAP_NORESERVE, fd, 0);
    if ((void *) data == MAP_FAILED) {
        err("failed to mmap() file %s", infile);
    }
    wr = blk_writer_init(outfile, mem->ioblksz);
    j = mem->ioblksz * (filesz / mem->ioblksz);
    if (j == filesz) {
        j -= mem->ioblksz;
    }
    for (i = filesz - 1; i > 0; --i) {
        blk_writer_append(wr, &data[i], 1, NULL);
        if (i == j) {
            size_t sz = filesz - j;
            if (sz > mem->ioblksz) {
                sz = mem->ioblksz;
            }
            madvise((void *) &data[i], sz, MADV_DONTNEED);
            j -= mem->ioblksz;
        }
    }
    blk_writer_append(wr, data, 1, NULL);
    blk_writer_close(wr, NULL);
    if (munmap((void *) data, mapsz) != 0) {
        err("munmap() failed");
    }
    if (close(fd) != 0) {
        err("close() failed");
    }
    if (unlink(infile) != 0) {
        err("failed to delete file %s", infile);
    }
    msg("\t%.3f sec total\n\n", now() - t);
}


/* ================================================================ */
/*        Phase 4: Convert spherical coords to unit vectors         */
/* ================================================================ */

static void spherical_to_vec(const char * const datafile,
                             const char * const scratchfile,
                             const struct mem_params * const mem,
                             const size_t npoints)
{
    const struct tree_entry *data;
    struct blk_writer *wr;
    void *behind;
    double t;
    size_t i;
    int fd;
    const size_t pagesz = (size_t) sysconf(_SC_PAGESIZE);
    size_t mapsz = sizeof(struct tree_entry) * npoints;

    if (mapsz == 0) {
        err("Empty data file");
    }
    msg("Converting spherical coordinates in %s to unit vectors\n", datafile);
    t = now();
    fd = open(datafile, O_RDONLY);
    if (fd == -1) {
        err("failed to open file %s for reading", datafile);
    }
    if (mapsz % pagesz != 0) {
        mapsz += pagesz - mapsz % pagesz;
    }
    behind = mmap(NULL, mapsz, PROT_READ, MAP_SHARED | MAP_NORESERVE, fd, 0);
    if (behind == MAP_FAILED) {
        err("failed to mmap() file %s", datafile);
    }
    if (madvise(behind, mapsz, MADV_SEQUENTIAL) != 0) {
        err("madvise() failed");
    }
    data = (const struct tree_entry *) behind;
    behind = ((unsigned char *) behind) + mem->ioblksz;
    wr = blk_writer_init(scratchfile, mem->ioblksz);
    for (i = 0; i < npoints; ++i) {
        struct htm_tree_entry e;
        if ((void *) &data[i] > behind) {
            void *ptr = ((unsigned char *) behind) - mem->ioblksz;
            if (madvise(ptr, mem->ioblksz, MADV_DONTNEED) != 0) {
                err("madvise() failed");
            }
            behind = ((unsigned char *) behind) + mem->ioblksz;
        }
        htm_sc_tov3(&e.v, &data[i].sc);
        e.rowid = data[i].rowid;
        blk_writer_append(wr, &e, sizeof(struct htm_tree_entry), NULL);
    }
    blk_writer_close(wr, NULL);
    if (munmap((void *) data, mapsz) != 0) {
        err("munmap() failed");
    }
    if (close(fd) != 0) {
        err("close() failed");
    }
    if (rename(scratchfile, datafile) != 0) {
        err("failed to rename file %s to %s", scratchfile, datafile);
    }
    msg("\t%.3f sec total\n\n", now() - t);
}
