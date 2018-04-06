#define N 12
#define SCOUNT PCOUNT

__constant__ uint32_t pow2[9] = {1, 2, 4, 8, 16, 32, 64, 128, 256};

__constant__ uint32_t binvert_limb_table[128] = {
  0x01, 0xAB, 0xCD, 0xB7, 0x39, 0xA3, 0xC5, 0xEF,
  0xF1, 0x1B, 0x3D, 0xA7, 0x29, 0x13, 0x35, 0xDF,
  0xE1, 0x8B, 0xAD, 0x97, 0x19, 0x83, 0xA5, 0xCF,
  0xD1, 0xFB, 0x1D, 0x87, 0x09, 0xF3, 0x15, 0xBF,
  0xC1, 0x6B, 0x8D, 0x77, 0xF9, 0x63, 0x85, 0xAF,
  0xB1, 0xDB, 0xFD, 0x67, 0xE9, 0xD3, 0xF5, 0x9F,
  0xA1, 0x4B, 0x6D, 0x57, 0xD9, 0x43, 0x65, 0x8F,
  0x91, 0xBB, 0xDD, 0x47, 0xC9, 0xB3, 0xD5, 0x7F,
  0x81, 0x2B, 0x4D, 0x37, 0xB9, 0x23, 0x45, 0x6F,
  0x71, 0x9B, 0xBD, 0x27, 0xA9, 0x93, 0xB5, 0x5F,
  0x61, 0x0B, 0x2D, 0x17, 0x99, 0x03, 0x25, 0x4F,
  0x51, 0x7B, 0x9D, 0x07, 0x89, 0x73, 0x95, 0x3F,
  0x41, 0xEB, 0x0D, 0xF7, 0x79, 0xE3, 0x05, 0x2F,
  0x31, 0x5B, 0x7D, 0xE7, 0x69, 0x53, 0x75, 0x1F,
  0x21, 0xCB, 0xED, 0xD7, 0x59, 0xC3, 0xE5, 0x0F,
  0x11, 0x3B, 0x5D, 0xC7, 0x49, 0x33, 0x55, 0xFF
};


typedef struct {
  uint32_t index;
  uint32_t hashid;
  uint8_t origin;
  uint8_t chainpos;
  uint8_t type;
  uint8_t reserved;
} fermat_t;

typedef struct {
  uint32_t N_;
  uint32_t SIZE_;
  uint32_t STRIPES_;
  uint32_t WIDTH_;
  uint32_t PCOUNT_;
  uint32_t TARGET_;
  uint32_t LIMIT13_;
  uint32_t LIMIT14_;
  uint32_t LIMIT15_;
} config_t;

__global__ void getconfig(config_t *conf)
{
  config_t c;
  c.N_ = N;
  c.SIZE_ = SIZE;
  c.STRIPES_ = STRIPES;
  c.WIDTH_ = WIDTH;
  c.PCOUNT_ = PCOUNT;
  c.TARGET_ = TARGET;
  c.LIMIT13_ = LIMIT13;
  c.LIMIT14_ = LIMIT14;
  c.LIMIT15_ = LIMIT15;
  *conf = c;
}

__device__ void shl32(uint32_t *data, unsigned size)
{
  #pragma unroll
  for (int j = size-1; j > 0; j--)
    data[j] = data[j-1];
  data[0] = 0;
}

__device__ void shr32(uint32_t *data, unsigned size)
{
  #pragma unroll
  for (int j = 1; j < size; j++)
    data[j-1] = data[j];
  data[size-1] = 0;
}

__device__ void shl(uint32_t *data, unsigned size, unsigned bits)
{
  #pragma unroll
  for(int i = size-1; i > 0; i--)
    data[i] = (data[i] << bits) | (data[i-1] >> (32-bits));
  
  data[0] = data[0] << bits;
}

__device__ void shr(uint32_t *data, unsigned size, unsigned bits)
{
  #pragma unroll
  for(int i = 0; i < size-1; i++)
    data[i] = (data[i] >> bits) | (data[i+1] << (32-bits));
  data[size-1] = data[size-1] >> bits;
}

__device__ void shlreg(uint32_t *data, unsigned size, unsigned bits)
{
  for (unsigned i = 0, ie = bits/32; i < ie; i++)
    shl32(data, size);
  
  if (bits%32)
    shl(data, size, bits%32);
}


__device__ void shrreg(uint32_t *data, unsigned size, unsigned bits)
{
  for (unsigned i = 0, ie = bits/32; i < ie; i++)
    shr32(data, size);
  if (bits%32)
    shr(data, size, bits%32);
}

__device__ uint32_t add128(uint4 *A, uint4 B)
{
//   *A += B; 
  A->x += B.x;
  A->y += B.y;
  A->z += B.z;
  A->w += B.w;
//   uint4 carry = -convert_uint4((*A) < B);
  uint4 carry = { -(A->x < B.x), -(A->y < B.y), -(A->z < B.z), -(A->w < B.w) };
  
  (*A).y += carry.x; carry.y += ((*A).y < carry.x);
  (*A).z += carry.y; carry.z += ((*A).z < carry.y);
  (*A).w += carry.z;
  return carry.w + ((*A).w < carry.z); 
}


__device__ uint32_t add128Carry(uint4 *A, uint4 B, uint32_t externalCarry)
{
//   *A += B;
  A->x += B.x;
  A->y += B.y;
  A->z += B.z;
  A->w += B.w;  
//   uint4 carry = -convert_uint4((*A) < B);
  uint4 carry = { -(A->x < B.x), -(A->y < B.y), -(A->z < B.z), -(A->w < B.w) };
  
  (*A).x += externalCarry; carry.x += ((*A).x < externalCarry);
  (*A).y += carry.x; carry.y += ((*A).y < carry.x);
  (*A).z += carry.y; carry.z += ((*A).z < carry.y);
  (*A).w += carry.z;
  return carry.w + ((*A).w < carry.z); 
}

__device__ uint32_t add256(uint4 *a0, uint4 *a1, uint4 b0, uint4 b1)
{
  return add128Carry(a1, b1, add128(a0, b0));
}

__device__ uint32_t add384(uint4 *a0, uint4 *a1, uint4 *a2, uint4 b0, uint4 b1, uint4 b2)
{
  return add128Carry(a2, b2, add128Carry(a1, b1, add128(a0, b0)));
}

__device__ uint32_t add512(uint4 *a0, uint4 *a1, uint4 *a2, uint4 *a3, uint4 b0, uint4 b1, uint4 b2, uint4 b3)
{
  return add128Carry(a3, b3, add128Carry(a2, b2, add128Carry(a1, b1, add128(a0, b0))));
}

__device__ uint32_t sub64Borrow(uint2 *A, uint2 B, uint32_t externalBorrow)
{
//   uint2 borrow = -convert_uint2((*A) < B);
  uint2 borrow = { -(A->x < B.x), -(A->y < B.y) };
//   *A -= B;
  A->x -= B.x;
  A->y -= B.y;
  
  borrow.x += (*A).x < externalBorrow; (*A).x -= externalBorrow;
  borrow.y += (*A).y < borrow.x; (*A).y -= borrow.x;
  return borrow.y;
}

__device__ uint32_t sub96Borrow(uint4 *A, uint4 B, uint32_t externalBorrow)
{
  //   uint2 borrow = -convert_uint2((*A) < B);
  uint4 borrow = {
    (*A).x < B.x,
      (*A).y < B.y,
      (*A).z < B.z,
      0
  };
  (*A).x -= B.x;
  (*A).y -= B.y;
  (*A).z -= B.z;
  
  borrow.x += (*A).x < externalBorrow; (*A).x -= externalBorrow;
  borrow.y += (*A).y < borrow.x; (*A).y -= borrow.x;
  borrow.z += (*A).z < borrow.y; (*A).z -= borrow.y;
  
  return borrow.z;
}

__device__ uint32_t sub128(uint4 *A, uint4 B)
{
  uint4 borrow = {
    (*A).x < B.x,
      (*A).y < B.y,
      (*A).z < B.z,
      (*A).w < B.w
  };
  (*A).x -= B.x;
  (*A).y -= B.y;
  (*A).z -= B.z;
  (*A).w -= B.w;  
  
  borrow.y += (*A).y < borrow.x; (*A).y -= borrow.x;
  borrow.z += (*A).z < borrow.y; (*A).z -= borrow.y;
  borrow.w += (*A).w < borrow.z; (*A).w -= borrow.z;
  return borrow.w;
}

__device__ uint32_t sub128Borrow(uint4 *A, uint4 B, uint32_t externalBorrow)
{
//   uint4 borrow = -convert_uint4((*A) < B);
  uint4 borrow = { -(A->x < B.x), -(A->y < B.y), -(A->z < B.z), -(A->w < B.w) };  
//   *A -= B;
  A->x -= B.x;
  A->y -= B.y;
  A->z -= B.z;
  A->w -= B.w;  
  
  borrow.x += (*A).x < externalBorrow; (*A).x -= externalBorrow;
  borrow.y += (*A).y < borrow.x; (*A).y -= borrow.x;
  borrow.z += (*A).z < borrow.y; (*A).z -= borrow.y;
  borrow.w += (*A).w < borrow.z; (*A).w -= borrow.z;
  return borrow.w;
}

__device__ uint32_t sub256(uint4 *a0, uint4 *a1, uint4 b0, uint4 b1)
{
  return sub128Borrow(a1, b1, sub128(a0, b0));
}

__device__ uint32_t sub320(uint4 *a0, uint4 *a1, uint2 *a2, uint4 b0, uint4 b1, uint2 b2)
{
  return sub64Borrow(a2, b2, sub128Borrow(a1, b1, sub128(a0, b0)));
}

__device__ uint32_t sub352(uint4 *a0, uint4 *a1, uint4 *a2, uint4 b0, uint4 b1, uint4 b2)
{
  return sub96Borrow(a2, b2, sub128Borrow(a1, b1, sub128(a0, b0)));
}

__device__ uint32_t sub384(uint4 *a0, uint4 *a1, uint4 *a2, uint4 b0, uint4 b1, uint4 b2)
{
  return sub128Borrow(a2, b2, sub128Borrow(a1, b1, sub128(a0, b0)));
}

__device__ uint32_t sub448(uint4 *a0, uint4 *a1, uint4 *a2, uint2 *a3, uint4 b0, uint4 b1, uint4 b2, uint2 b3)
{
  return sub64Borrow(a3, b3, sub128Borrow(a2, b2, sub128Borrow(a1, b1, sub128(a0, b0))));
}

__device__ uint32_t invert_limb(uint32_t limb)
{
  uint32_t inv = binvert_limb_table[(limb/2) & 0x7F];
  inv = 2*inv - inv*inv*limb;
  inv = 2*inv - inv*inv*limb;
  return -inv;
}
