/** \file
    \brief      Variable length integer encoding/decoding.

    The scheme employed here is as follows: an unsigned integer is mapped
    to a sequence of bytes. The 3 MSBs of the first byte contain the
    number of subsequent bytes in the sequence, and the 5 LSBs contain
    the 5 MSBs of the value. Bytes after the first contain the remaining bits
    of the value, from most to least significant. The smallest representable
    integer is 0, the largest 2^61 - 1.

    \authors    Serge Monkewitz
    \copyright  IPAC/Caltech
  */
#ifndef HTM_VARINT_H
#define HTM_VARINT_H

#include "common.h"

#ifdef __cplusplus
extern "C" {
#endif

/** \defgroup varint Variable length integer encoding/decoding
    @{
  */

/** Returns the number of bytes required to represent val.
  */
HTM_INLINE unsigned int htm_varint_len(uint64_t val)
{
#ifdef __GNUC__
    return val == 0 ? 1 : 1 + (66 - __builtin_clzll(val))/8;
#else
    unsigned int n = 3;
    uint64_t v = val;
    while ((v >>= 1) != 0) { ++n; }
    return 1 + n/8;
#endif
}

/** Returns the number of bytes following the specified leading byte
    (in a variable length coding of an unsigned integer).
  */
HTM_INLINE unsigned int htm_varint_nfollow(unsigned char leadbyte)
{
    return leadbyte >> 5;
}

/** Decodes the variable length integer encoded at buf.
  */
HTM_INLINE uint64_t htm_varint_decode(const unsigned char *buf)
{
    uint64_t val = buf[0] & 0x1f;
    unsigned int nfollow = htm_varint_nfollow(buf[0]);
    switch (buf[0] >> 5) {
        case 7: val = (val << 8) | buf[nfollow - 6];
        case 6: val = (val << 8) | buf[nfollow - 5];
        case 5: val = (val << 8) | buf[nfollow - 4];
        case 4: val = (val << 8) | buf[nfollow - 3];
        case 3: val = (val << 8) | buf[nfollow - 2];
        case 2: val = (val << 8) | buf[nfollow - 1];
        case 1: val = (val << 8) | buf[nfollow];
        case 0: break;
    }
    return val;
}

/** Variable-length encodes val into the given buffer, which must point to
    a buffer of at least 8 bytes. The number of bytes actually used in the
    encoding is returned. 
  */
HTM_INLINE unsigned int htm_varint_encode(unsigned char *buf,
                                          const uint64_t val)
{
    unsigned int n = htm_varint_len(val) - 1;
    buf[0] = (n << 5) | ((val >> 8*n) & 0xff);
    switch (n) {
        case 7: buf[n - 6] = (val >> 48) & 0xff;
        case 6: buf[n - 5] = (val >> 40) & 0xff;
        case 5: buf[n - 4] = (val >> 32) & 0xff;
        case 4: buf[n - 3] = (val >> 24) & 0xff;
        case 3: buf[n - 2] = (val >> 16) & 0xff;
        case 2: buf[n - 1] = (val >>  8) & 0xff;
        case 1: buf[n] = val & 0xff;
        case 0: break;
    }
    return n + 1;
}

/** Variable-length encodes val into the given buffer, which must point to
    a buffer of at least 8 bytes. Bytes are written in the reverse order
    of htm_varint_encode(). The number of bytes actually used in the
    encoding is returned.
  */
HTM_INLINE unsigned int htm_varint_rencode(unsigned char *buf,
                                           const uint64_t val)
{
    const unsigned int m = htm_varint_len(val);
    unsigned int n = m - 1;
    buf[n] = (n << 5) | ((val >> 8*n) & 0xff);
    switch (n) {
        case 7: --n; buf[n] = (val >> 48) & 0xff;
        case 6: --n; buf[n] = (val >> 40) & 0xff;
        case 5: --n; buf[n] = (val >> 32) & 0xff;
        case 4: --n; buf[n] = (val >> 24) & 0xff;
        case 3: --n; buf[n] = (val >> 16) & 0xff;
        case 2: --n; buf[n] = (val >>  8) & 0xff;
        case 1: --n; buf[n] = val & 0xff;
        case 0: break;
    }
    return m;
}


/** @}
  */

#ifdef __cplusplus
}
#endif

#endif /* HTM_VARINT_H */

