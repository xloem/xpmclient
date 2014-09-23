
void mul320schoolBook_v3(uint4 op1l0, uint4 op1l1, uint2 op1l2, // low --> hi
                         uint4 op2l0, uint4 op2l1, uint2 op2l2, // low --> hi
                         uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3, uint4 *rl4)
{ 
#define b1 op2l2.y
#define b2 op2l2.x
#define b3 op2l1.w
#define b4 op2l1.z
#define b5 op2l1.y
#define b6 op2l1.x
#define b7 op2l0.w
#define b8 op2l0.z
#define b9 op2l0.y
#define b10 op2l0.x

  ulong R0x = op1l0.x * b10;
  ulong R0y = op1l0.y * b10;
  ulong R0z = op1l0.z * b10;
  ulong R0w = op1l0.w * b10;
  ulong R1x = op1l1.x * b10;
  ulong R1y = op1l1.y * b10;
  ulong R1z = op1l1.z * b10;
  ulong R1w = op1l1.w * b10;  
  ulong R2x = op1l2.x * b10;
  ulong R2y = op1l2.y * b10;
  ulong R2z = mul_hi(op1l0.x, b1);
  ulong R2w = mul_hi(op1l0.y, b1);
  ulong R3x = mul_hi(op1l0.z, b1);
  ulong R3y = mul_hi(op1l0.w, b1);
  ulong R3z = mul_hi(op1l1.x, b1);
  ulong R3w = mul_hi(op1l1.y, b1);
  ulong R4x = mul_hi(op1l1.z, b1);
  ulong R4y = mul_hi(op1l1.w, b1);
  ulong R4z = mul_hi(op1l2.x, b1);
  ulong R4w = mul_hi(op1l2.y, b1);
  
  mul320round_v3(op1l0, op1l1, op1l2, b10, b9, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z);
  mul320round_v3(op1l0, op1l1, op1l2, b9, b8, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w);
  mul320round_v3(op1l0, op1l1, op1l2, b8,  b7, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x);
  mul320round_v3(op1l0, op1l1, op1l2,  b7,  b6, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y);
  mul320round_v3(op1l0, op1l1, op1l2,  b6,  b5, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);
  mul320round_v3(op1l0, op1l1, op1l2,  b5,  b4, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w);
  mul320round_v3(op1l0, op1l1, op1l2,  b4,  b3, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x);
  mul320round_v3(op1l0, op1l1, op1l2,  b3,  b2, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y);
  mul320round_v3(op1l0, op1l1, op1l2,  b2,  b1, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z);

  union {
    uint2 v32;
    ulong v64;
  } Int;
  
  Int.v64 = R2z; R2w += Int.v32.y;
  Int.v64 = R2w; R3x += Int.v32.y;  
  Int.v64 = R3x; R3y += Int.v32.y;
  Int.v64 = R3y; R3z += Int.v32.y;
  Int.v64 = R3z; R3w += Int.v32.y;
  Int.v64 = R3w; R4x += Int.v32.y;
  Int.v64 = R4x; R4y += Int.v32.y;
  Int.v64 = R4y; R4z += Int.v32.y;
  Int.v64 = R4z; R4w += Int.v32.y;
  
  *rl0 = (uint4){R0x, R0y, R0z, R0w};
  *rl1 = (uint4){R1x, R1y, R1z, R1w};
  *rl2 = (uint4){R2x, R2y, R2z, R2w};
  *rl3 = (uint4){R3x, R3y, R3z, R3w};
  *rl4 = (uint4){R4x, R4y, R4z, R4w};
  
#undef b1
#undef b2
#undef b3
#undef b4
#undef b5
#undef b6
#undef b7
#undef b8
#undef b9
#undef b10
}

void mul352schoolBook_v3(uint4 op1l0, uint4 op1l1, uint4 op1l2, // low --> hi
                         uint4 op2l0, uint4 op2l1, uint4 op2l2, // low --> hi
                         uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3, uint4 *rl4, uint4 *rl5)
{ 
  #define b1 op2l2.z
  #define b2 op2l2.y
  #define b3 op2l2.x
  #define b4 op2l1.w
  #define b5 op2l1.z
  #define b6 op2l1.y
  #define b7 op2l1.x
  #define b8 op2l0.w
  #define b9 op2l0.z
  #define b10 op2l0.y
  #define b11 op2l0.x  
  ulong R0x = op1l0.x * b11;
  ulong R0y = op1l0.y * b11;
  ulong R0z = op1l0.z * b11;
  ulong R0w = op1l0.w * b11;
  ulong R1x = op1l1.x * b11;
  ulong R1y = op1l1.y * b11;
  ulong R1z = op1l1.z * b11;
  ulong R1w = op1l1.w * b11;  
  ulong R2x = op1l2.x * b11;
  ulong R2y = op1l2.y * b11;
  ulong R2z = op1l2.z * b11;
  ulong R2w = mul_hi(op1l0.x, b1);
  ulong R3x = mul_hi(op1l0.y, b1);
  ulong R3y = mul_hi(op1l0.z, b1);
  ulong R3z = mul_hi(op1l0.w, b1);
  ulong R3w = mul_hi(op1l1.x, b1);
  ulong R4x = mul_hi(op1l1.y, b1);
  ulong R4y = mul_hi(op1l1.z, b1);
  ulong R4z = mul_hi(op1l1.w, b1);
  ulong R4w = mul_hi(op1l2.x, b1);
  ulong R5x = mul_hi(op1l2.y, b1);
  ulong R5y = mul_hi(op1l2.z, b1);
  
  mul352round_v3(op1l0, op1l1, op1l2, b11, b10, &R0y, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w);
  mul352round_v3(op1l0, op1l1, op1l2, b10, b9, &R0z, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x);
  mul352round_v3(op1l0, op1l1, op1l2, b9,  b8, &R0w, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y);
  mul352round_v3(op1l0, op1l1, op1l2,  b8,  b7, &R1x, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z);
  mul352round_v3(op1l0, op1l1, op1l2,  b7,  b6, &R1y, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w);
  mul352round_v3(op1l0, op1l1, op1l2,  b6,  b5, &R1z, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x);
  mul352round_v3(op1l0, op1l1, op1l2,  b5,  b4, &R1w, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y);
  mul352round_v3(op1l0, op1l1, op1l2,  b4,  b3, &R2x, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z);
  mul352round_v3(op1l0, op1l1, op1l2,  b3,  b2, &R2y, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w);
  mul352round_v3(op1l0, op1l1, op1l2,  b2,  b1, &R2z, &R2w, &R3x, &R3y, &R3z, &R3w, &R4x, &R4y, &R4z, &R4w, &R5x);
  
  union {
    uint2 v32;
    ulong v64;
  } Int;
  
  Int.v64 = R2w; R3x += Int.v32.y;  
  Int.v64 = R3x; R3y += Int.v32.y;
  Int.v64 = R3y; R3z += Int.v32.y;
  Int.v64 = R3z; R3w += Int.v32.y;
  Int.v64 = R3w; R4x += Int.v32.y;
  Int.v64 = R4x; R4y += Int.v32.y;
  Int.v64 = R4y; R4z += Int.v32.y;
  Int.v64 = R4z; R4w += Int.v32.y;
  Int.v64 = R4w; R5x += Int.v32.y;
  Int.v64 = R5x; R5y += Int.v32.y;
  Int.v64 = R5y;
  
  *rl0 = (uint4){R0x, R0y, R0z, R0w};
  *rl1 = (uint4){R1x, R1y, R1z, R1w};
  *rl2 = (uint4){R2x, R2y, R2z, R2w};
  *rl3 = (uint4){R3x, R3y, R3z, R3w};
  *rl4 = (uint4){R4x, R4y, R4z, R4w};
  *rl5 = (uint4){R5x, R5y, 0, 0};
  
  #undef b1
  #undef b2
  #undef b3
  #undef b4
  #undef b5
  #undef b6
  #undef b7
  #undef b8
  #undef b9
  #undef b10
  #undef b11
  #undef b12
}

void sqr320(uint4 op1l0, uint4 op1l1, uint2 op1l2,
            uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3, uint4 *rl4)
{
#define a1 op1l2.y
#define a2 op1l2.x
#define a3 op1l1.w
#define a4 op1l1.z
#define a5 op1l1.y
#define a6 op1l1.x
#define a7 op1l0.w
#define a8 op1l0.z
#define a9 op1l0.y
#define a10 op1l0.x
  
  uint64_t R0x = 0, R0y, R0z, R0w, R1x, R1y, R1z, R1w, R2x,
           R2y, R2z, R2w, R3x, R3y, R3z, R3w, R4x,
           R4y, R4z, R4w;

  m_add(a10, a10, 1, &R0x, &R0y, 1);
  m_add(a10, a9, 2, &R0y, &R0z, 1);
  m_dbladd(a10, a8, 2, a9, a9, 1, &R0z, &R0w, 1);
  m_dbladd(a10, a7, 2, a9, a8, 2, &R0w, &R1x, 1);
  m_dbladd(a10, a6, 2, a9, a7, 2, &R1x, &R1y, 1);
      m_add(a8, a8, 1, &R1x, &R1y, 0);
  m_dbladd(a10, a5, 2, a9, a6, 2, &R1y, &R1z, 1);
      m_add(a8, a7, 2, &R1y, &R1z, 0);
  m_dbladd(a10, a4, 2, a9, a5, 2, &R1z, &R1w, 1);
      m_dbladd(a8, a6, 2, a7, a7, 1, &R1z, &R1w, 0);
  m_dbladd(a10, a3, 2, a9, a4, 2, &R1w, &R2x, 1);
      m_dbladd(a8, a5, 2, a7, a6, 2, &R1w, &R2x, 0);
  m_dbladd(a10, a2, 2, a9, a3, 2, &R2x, &R2y, 1);
      m_dbladd(a8, a4, 2, a7, a5, 2, &R2x, &R2y, 0);
      m_add(a6, a6, 1, &R2x, &R2y, 0);
  m_dbladd(a10, a1, 2, a9, a2, 2, &R2y, &R2z, 1);
      m_dbladd(a8, a3, 2, a7, a4, 2, &R2y, &R2z, 0);
      m_add(a6, a5, 2, &R2y, &R2z, 0);
  m_dbladd(a9, a1, 2, a8, a2, 2, &R2z, &R2w, 1);
      m_dbladd(a7, a3, 2, a6, a4, 2, &R2z, &R2w, 0);
      m_add(a5, a5, 1, &R2z, &R2w, 0);
  m_dbladd(a8, a1, 2, a7, a2, 2, &R2w, &R3x, 1);
      m_dbladd(a6, a3, 2, a5, a4, 2, &R2w, &R3x, 0);
  m_dbladd(a7, a1, 2, a6, a2, 2, &R3x, &R3y, 1);
      m_dbladd(a5, a3, 2, a4, a4, 1, &R3x, &R3y, 0);
  m_dbladd(a6, a1, 2, a5, a2, 2, &R3y, &R3z, 1);
      m_add(a4, a3, 2, &R3y, &R3z, 0);
  m_dbladd(a5, a1, 2, a4, a2, 2, &R3z, &R3w, 1);
      m_add(a3, a3, 1, &R3z, &R3w, 0);
  m_dbladd(a4, a1, 2, a3, a2, 2, &R3w, &R4x, 1);
  m_dbladd(a3, a1, 2, a2, a2, 1, &R4x, &R4y, 1);
  m_add(a2, a1, 2, &R4y, &R4z, 1);
  m_add(a1, a1, 1, &R4z, &R4w, 1);
  
  union {
    uint2 v32;
    ulong v64;
  } Int;
  
  Int.v64 = R0y; R0z += Int.v32.y;
  Int.v64 = R0z; R0w += Int.v32.y;
  Int.v64 = R0w; R1x += Int.v32.y;
  Int.v64 = R1x; R1y += Int.v32.y;
  Int.v64 = R1y; R1z += Int.v32.y;
  Int.v64 = R1z; R1w += Int.v32.y;
  Int.v64 = R1w; R2x += Int.v32.y;
  Int.v64 = R2x; R2y += Int.v32.y;
  Int.v64 = R2y; R2z += Int.v32.y;
  Int.v64 = R2z; R2w += Int.v32.y;
  Int.v64 = R2w; R3x += Int.v32.y;
  Int.v64 = R3x; R3y += Int.v32.y;
  Int.v64 = R3y; R3z += Int.v32.y;
  Int.v64 = R3z; R3w += Int.v32.y;
  Int.v64 = R3w; R4x += Int.v32.y;
  Int.v64 = R4x; R4y += Int.v32.y;
  Int.v64 = R4y; R4z += Int.v32.y;
  Int.v64 = R4z; R4w += Int.v32.y;
  
  *rl0 = (uint4){R0x, R0y, R0z, R0w};
  *rl1 = (uint4){R1x, R1y, R1z, R1w};
  *rl2 = (uint4){R2x, R2y, R2z, R2w};
  *rl3 = (uint4){R3x, R3y, R3z, R3w};
  *rl4 = (uint4){R4x, R4y, R4z, R4w};
  
#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef a7
#undef a8
#undef a9
#undef a10
}


void sqr352(uint4 op1l0, uint4 op1l1, uint4 op1l2,
            uint4 *rl0, uint4 *rl1, uint4 *rl2, uint4 *rl3, uint4 *rl4, uint4 *rl5)
{
#define a1 op1l2.z
#define a2 op1l2.y
#define a3 op1l2.x
#define a4 op1l1.w
#define a5 op1l1.z
#define a6 op1l1.y
#define a7 op1l1.x
#define a8 op1l0.w
#define a9 op1l0.z
#define a10 op1l0.y
#define a11 op1l0.x
  
  uint64_t R0 = 0, R1, R2, R3, R4, R5, R6, R7, R8,
           R9, R10, R11, R12, R13, R14, R15, R16,
           R17, R18, R19, R20, R21;

  m_add(a11, a11, 1, &R0, &R1, 1);
  m_add(a11, a10, 2, &R1, &R2, 1);
  m_dbladd(a11, a9, 2, a10, a10, 1, &R2, &R3, 1);
  m_dbladd(a11, a8, 2, a10, a9, 2, &R3, &R4, 1);
  m_dbladd(a11, a7, 2, a10, a8, 2, &R4, &R5, 1);
      m_add(a9, a9, 1, &R4, &R5, 0);
  m_dbladd(a11, a6, 2, a10, a7, 2, &R5, &R6, 1);
      m_add(a9, a8, 2, &R5, &R6, 0);
  m_dbladd(a11, a5, 2, a10, a6, 2, &R6, &R7, 1);
      m_dbladd(a9, a7, 2, a8, a8, 1, &R6, &R7, 0);
  m_dbladd(a11, a4, 2, a10, a5, 2, &R7, &R8, 1);
      m_dbladd(a9, a6, 2, a8, a7, 2, &R7, &R8, 0);
  m_dbladd(a11, a3, 2, a10, a4, 2, &R8, &R9, 1);
      m_dbladd(a9, a5, 2, a8, a6, 2, &R8, &R9, 0);
      m_add(a7, a7, 1, &R8, &R9, 0);
  m_dbladd(a11, a2, 2, a10, a3, 2, &R9, &R10, 1);
      m_dbladd(a9, a4, 2, a8, a5, 2, &R9, &R10, 0);
      m_add(a7, a6, 2, &R9, &R10, 0);
  m_dbladd(a11, a1, 2, a10, a2, 2, &R10, &R11, 1);
      m_dbladd(a9, a3, 2, a8, a4, 2, &R10, &R11, 0);
      m_dbladd(a7, a5, 2, a6, a6, 1, &R10, &R11, 0);
  m_dbladd(a10, a1, 2, a9, a2, 2, &R11, &R12, 1);
      m_dbladd(a8, a3, 2, a7, a4, 2, &R11, &R12, 0);
      m_add(a6, a5, 2, &R11, &R12, 0);
  m_dbladd(a9, a1, 2, a8, a2, 2, &R12, &R13, 1);
      m_dbladd(a7, a3, 2, a6, a4, 2, &R12, &R13, 0);
      m_add(a5, a5, 1, &R12, &R13, 0);
  m_dbladd(a8, a1, 2, a7, a2, 2, &R13, &R14, 1);
      m_dbladd(a6, a3, 2, a5, a4, 2, &R13, &R14, 0);
  m_dbladd(a7, a1, 2, a6, a2, 2, &R14, &R15, 1);
      m_dbladd(a5, a3, 2, a4, a4, 1, &R14, &R15, 0);
  m_dbladd(a6, a1, 2, a5, a2, 2, &R15, &R16, 1);
      m_add(a4, a3, 2, &R15, &R16, 0);
  m_dbladd(a5, a1, 2, a4, a2, 2, &R16, &R17, 1);
      m_add(a3, a3, 1, &R16, &R17, 0);
  m_dbladd(a4, a1, 2, a3, a2, 2, &R17, &R18, 1);
  m_dbladd(a3, a1, 2, a2, a2, 1, &R18, &R19, 1);
  m_add(a2, a1, 2, &R19, &R20, 1);
  m_add(a1, a1, 1, &R20, &R21, 1);
  
  union {
    uint2 v32;
    ulong v64;
  } Int;
  
  Int.v64 = R1; R2 += Int.v32.y;
  Int.v64 = R2; R3 += Int.v32.y;
  Int.v64 = R3; R4 += Int.v32.y;
  Int.v64 = R4; R5 += Int.v32.y;
  Int.v64 = R5; R6 += Int.v32.y;
  Int.v64 = R6; R7 += Int.v32.y;
  Int.v64 = R7; R8 += Int.v32.y;
  Int.v64 = R8; R9 += Int.v32.y;
  Int.v64 = R9; R10 += Int.v32.y;
  Int.v64 = R10; R11 += Int.v32.y;
  Int.v64 = R11; R12 += Int.v32.y;
  Int.v64 = R12; R13 += Int.v32.y;
  Int.v64 = R13; R14 += Int.v32.y;
  Int.v64 = R14; R15 += Int.v32.y;
  Int.v64 = R15; R16 += Int.v32.y;
  Int.v64 = R16; R17 += Int.v32.y;
  Int.v64 = R17; R18 += Int.v32.y;
  Int.v64 = R18; R19 += Int.v32.y;
  Int.v64 = R19; R20 += Int.v32.y;
  Int.v64 = R20; R21 += Int.v32.y;
  
  *rl0 = (uint4){R0, R1, R2, R3};
  *rl1 = (uint4){R4, R5, R6, R7};
  *rl2 = (uint4){R8, R9, R10, R11};
  *rl3 = (uint4){R12, R13, R14, R15};
  *rl4 = (uint4){R16, R17, R18, R19};
  *rl5 = (uint4){R20, R21, 0, 0};
  
#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef a7
#undef a8
#undef a9
#undef a10
#undef a11
#undef a12
}


__kernel void squareBenchmark320(__global uint32_t *m1,
                                 __global uint32_t *out,
                                 unsigned elementsNum)
{
#define OperandSize 10
#define GmpOperandSize 10
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint32_t op1[OperandSize];
    for (unsigned j = 0; j < OperandSize; j++)
      op1[j] = m1[i*GmpOperandSize + j];  
    
    uint4 result[6];    
    
    uint4 op1v[3] = {
      (uint4){op1[0], op1[1], op1[2], op1[3]}, 
      (uint4){op1[4], op1[5], op1[6], op1[7]}, 
      (uint4){op1[8], op1[9], 0, 0}
    };

    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      sqr320(op1v[0], op1v[1], op1v[2].xy,
             &result[0], &result[1], &result[2], &result[3], &result[4]);
      op1v[0] = result[0];
      op1v[1] = result[1];
      op1v[2].xy = result[2].xy;
    }
    
    uint32_t *pResult = (uint32_t*)result;
    for (unsigned j = 0; j < OperandSize*2; j++)
      out[i*OperandSize*2 + j] = pResult[j];
  }
#undef GmpOperandSize
#undef OperandSize
}

__kernel void squareBenchmark352(__global uint32_t *m1,
                                 __global uint32_t *out,
                                 unsigned elementsNum)
{
#define OperandSize 11
#define GmpOperandSize 12  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint32_t op1[OperandSize];
    for (unsigned j = 0; j < OperandSize; j++)
      op1[j] = m1[i*GmpOperandSize + j];  
    
    uint4 result[6];    
    
    uint4 op1v[3] = {
      (uint4){op1[0], op1[1], op1[2], op1[3]}, 
      (uint4){op1[4], op1[5], op1[6], op1[7]}, 
      (uint4){op1[8], op1[9], op1[10], 0}
    };

    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      sqr352(op1v[0], op1v[1], op1v[2],
             &result[0], &result[1], &result[2], &result[3], &result[4], &result[5]);
      op1v[0] = result[0];
      op1v[1] = result[1];
      op1v[2].xyz = result[2].xyz;
    }
    
    uint32_t *pResult = (uint32_t*)result;
    for (unsigned j = 0; j < OperandSize*2; j++)
      out[i*OperandSize*2 + j] = pResult[j];
  }
#undef GmpOperandSize
#undef OperandSize
}


__kernel void multiplyBenchmark320(__global uint32_t *m1,
                                   __global uint32_t *m2,
                                   __global uint32_t *out,
                                   unsigned elementsNum)
{
#define OperandSize 10
#define GmpOperandSize 10  
  unsigned globalSize = get_global_size(0);

  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint32_t op1[OperandSize];
    uint32_t op2[OperandSize];
    for (unsigned j = 0; j < OperandSize; j++)
      op1[j] = m1[i*GmpOperandSize + j];
    for (unsigned j = 0; j < OperandSize; j++)
      op2[j] = m2[i*GmpOperandSize + j];    
    
    uint4 result[6];    
    
    uint4 op1v[3] = {
      (uint4){op1[0], op1[1], op1[2], op1[3]}, 
      (uint4){op1[4], op1[5], op1[6], op1[7]}, 
      (uint4){op1[8], op1[9], 0, 0}
    };
    
    uint4 op2v[3] = {
      (uint4){op2[0], op2[1], op2[2], op2[3]}, 
      (uint4){op2[4], op2[5], op2[6], op2[7]}, 
      (uint4){op2[8], op2[9], 0, 0}
    };   

    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      mul320schoolBook_v3(op1v[0], op1v[1], op1v[2].xy, op2v[0], op2v[1], op2v[2].xy,
                          &result[0], &result[1], &result[2], &result[3], &result[4]);
      op1v[0] = result[0];
      op1v[1] = result[1];
      op1v[2].xy = result[2].xy;
    }
    
    uint32_t *pResult = (uint32_t*)result;
    for (unsigned j = 0; j < OperandSize*2; j++)
      out[i*OperandSize*2 + j] = pResult[j];
  }
#undef GmpOperandSize
#undef OperandSize
}

__kernel void multiplyBenchmark352(__global uint32_t *m1,
                                   __global uint32_t *m2,
                                   __global uint32_t *out,
                                   unsigned elementsNum)
{
#define OperandSize 11
#define GmpOperandSize 12  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint32_t op1[OperandSize];
    uint32_t op2[OperandSize];
    for (unsigned j = 0; j < OperandSize; j++)
      op1[j] = m1[i*GmpOperandSize + j];
    for (unsigned j = 0; j < OperandSize; j++)
      op2[j] = m2[i*GmpOperandSize + j];    
    
    uint4 result[6];    
    
    uint4 op1v[3] = {
      (uint4){op1[0], op1[1], op1[2], op1[3]}, 
      (uint4){op1[4], op1[5], op1[6], op1[7]}, 
      (uint4){op1[8], op1[9], op1[10], 0}
    };
    
    uint4 op2v[3] = {
      (uint4){op2[0], op2[1], op2[2], op2[3]}, 
      (uint4){op2[4], op2[5], op2[6], op2[7]}, 
      (uint4){op2[8], op2[9], op2[10], 0}
    };   

    for (unsigned repeatNum = 0; repeatNum < 512; repeatNum++) {
      mul352schoolBook_v3(op1v[0], op1v[1], op1v[2], op2v[0], op2v[1], op2v[2],
                          &result[0], &result[1], &result[2], &result[3], &result[4], &result[5]);
      op1v[0] = result[0];
      op1v[1] = result[1];
      op1v[2].xyz = result[2].xyz;
    }
    
    uint32_t *pResult = (uint32_t*)&result;
    for (unsigned j = 0; j < OperandSize*2; j++)
      out[i*OperandSize*2 + j] = pResult[j];
  }
  
#undef GmpOperandSize
#undef OperandSize
}


__kernel void fermatTestBenchMark320(__global uint32_t *restrict numbers,
                                     __global uint32_t *restrict out,
                                     unsigned elementsNum)
{
#define OperandSize 10  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint32_t lNumbers[OperandSize];
    for (unsigned j = 0; j < OperandSize; j++)
      lNumbers[j] = numbers[i*OperandSize+j];

    uint4 result[3];
    
    uint4 lNumbersv[3] = {
      (uint4){lNumbers[0], lNumbers[1], lNumbers[2], lNumbers[3]}, 
      (uint4){lNumbers[4], lNumbers[5], lNumbers[6], lNumbers[7]}, 
      (uint4){lNumbers[8], lNumbers[9], 0, 0}
    };    

    FermatTest320(lNumbersv, &result[0], &result[1], &result[2]);
      
    uint32_t *pResult = (uint32_t*)result;    
    for (unsigned j = 0; j < OperandSize; j++)
      out[i*OperandSize + j] = pResult[j];  
  }
#undef OperandSize
}


__kernel void fermatTestBenchMark352(__global uint32_t *restrict numbers,
                                     __global uint32_t *restrict out,
                                     unsigned elementsNum)
{
#define OperandSize 11  
  unsigned globalSize = get_global_size(0);
  for (unsigned i = get_global_id(0); i < elementsNum; i += globalSize) {
    uint32_t lNumbers[OperandSize];
    for (unsigned j = 0; j < OperandSize; j++)
      lNumbers[j] = numbers[i*OperandSize+j];

    uint4 result[3];
    
    uint4 lNumbersv[3] = {
      (uint4){lNumbers[0], lNumbers[1], lNumbers[2], lNumbers[3]}, 
      (uint4){lNumbers[4], lNumbers[5], lNumbers[6], lNumbers[7]}, 
      (uint4){lNumbers[8], lNumbers[9], lNumbers[10], 0}
    };    

    FermatTest352(lNumbersv, &result[0], &result[1], &result[2]);
      
    uint32_t *pResult = (uint32_t*)result;    
    for (unsigned j = 0; j < OperandSize; j++)
      out[i*OperandSize + j] = pResult[j];  
  }
#undef OperandSize
}
