#ifndef __CLCOMMON_H_
#define __CLCOMMON_H_

typedef unsigned char uint8_t;
typedef unsigned int uint32_t;
typedef unsigned long uint64_t;

typedef struct {
  uint index;
  uint hashid;
  uchar origin;
  uchar chainpos;
  uchar type;
  uchar reserved;
} fermat_t;


typedef struct {
  uint N_;
  uint SIZE_;
  uint STRIPES_;
  uint WIDTH_;
  uint PCOUNT_;
  uint TARGET_;
  uint LIMIT13_;
  uint LIMIT14_;
  uint LIMIT15_;
} config_t;

#define N 12

#if defined(__gfx900__) ||\
    defined(__gfx901__) ||\
    defined(__gfx902__) ||\
    defined(__gfx903__) ||\
    defined(__gfx904__) ||\
    defined(__gfx905__) ||\
    defined(__gfx906__) ||\
    defined(__gfx907__)
#define __AMDVEGA
#endif

#endif //__CLCOMMON_H_
