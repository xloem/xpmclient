#include "../baseclient.h"
#include "../uint256.h"
#include "../opencl.h"

#include <list>


#define PARAM_N       200
#define PARAM_K       9
#define PREFIX                          (PARAM_N / (PARAM_K + 1))
#define NR_INPUTS                       (1 << PREFIX)
// Approximate log base 2 of number of elements in hash tables
#define APX_NR_ELMS_LOG                 (PREFIX + 1)
// Number of rows and slots is affected by this. 20 offers the best performance
// but occasionally misses ~1% of solutions.
#define NR_ROWS_LOG                     20

// Make hash tables OVERHEAD times larger than necessary to store the average
// number of elements per row. The ideal value is as small as possible to
// reduce memory usage, but not too small or else elements are dropped from the
// hash tables.
//
// The actual number of elements per row is closer to the theoretical average
// (less variance) when NR_ROWS_LOG is small. So accordingly OVERHEAD can be
// smaller.
//
// Even (as opposed to odd) values of OVERHEAD sometimes significantly decrease
// performance as they cause VRAM channel conflicts.
#if NR_ROWS_LOG == 16
#define OVERHEAD                        3
#elif NR_ROWS_LOG == 18
#define OVERHEAD                        3
#elif NR_ROWS_LOG == 19
#define OVERHEAD                        5
#elif NR_ROWS_LOG == 20
#define OVERHEAD                        9
#endif

#define NR_ROWS                         (1 << NR_ROWS_LOG)
#define NR_SLOTS            ((1 << (APX_NR_ELMS_LOG - NR_ROWS_LOG)) * OVERHEAD)
// Length of 1 element (slot) in bytes
#define SLOT_LEN                        32
// Total size of hash table
#define HT_SIZE       (NR_ROWS * NR_SLOTS * SLOT_LEN)
// Length of Zcash block header and nonce
#define ZCASH_BLOCK_HEADER_LEN    140
#define ZCASH_NONCE_LEN     32
// Number of bytes Zcash needs out of Blake
#define ZCASH_HASH_LEN                  50
// Number of wavefronts per SIMD for the Blake kernel.
// Blake is ALU-bound (beside the atomic counter being incremented) so we need
// at least 2 wavefronts per SIMD to hide the 2-clock latency of integer
// instructions. 10 is the max supported by the hw.
#define BLAKE_WPS                 10
#define MAX_SOLS      2000


#define COLLISION_BIT_LENGTH (WN / (WK+1))
#define COLLISION_BYTE_LENGTH ((COLLISION_BIT_LENGTH+7)/8)
#define FINAL_FULL_WIDTH (2*COLLISION_BYTE_LENGTH+sizeof(uint32_t)*(1 << (WK)))


#define NDIGITS   (WK+1)
#define DIGITBITS (WN/(NDIGITS))
#define PROOFSIZE (1u<<WK)
#define COMPRESSED_PROOFSIZE ((COLLISION_BIT_LENGTH+1)*PROOFSIZE*4/(8*sizeof(uint32_t)))

#define HASHESPERBLAKE (512/WN)
#define HASHOUT (HASHESPERBLAKE*WN/8)

// Optional features
#undef ENABLE_DEBUG

/*
** Return the offset of Xi in bytes from the beginning of the slot.
*/
#define xi_offset_for_round(round)  (8 + ((round) / 2) * 4)

// An (uncompressed) solution stores (1 << PARAM_K) 32-bit values
#define SOL_SIZE      ((1 << PARAM_K) * 4)

typedef struct  sols_s
{
    uint32_t  nr;
    uint32_t  likely_invalidss;
    uint8_t valid[MAX_SOLS];
    uint32_t  values[MAX_SOLS][(1 << PARAM_K)];
}   sols_t;


typedef struct  debug_s
{
    uint32_t    dropped_coll;
    uint32_t    dropped_stor;
}               debug_t;


struct stats_t {
  unsigned id;
  double sols;
  
  stats_t(){
    id = 0;
    sols = 0;
  }
  
};

struct SHA_256;

struct MinerInstance {
  cl_context _context;
  cl_program _program;
  
  cl_command_queue queue;
  clBuffer<uint8_t> buf_ht0;
  clBuffer<uint8_t> buf_ht1;  
  clBuffer<sols_t> buf_sols;
  clBuffer<debug_t> buf_dbg;
  
  cl_kernel k_init_ht;
  cl_kernel k_rounds[PARAM_K];
  cl_kernel k_sols;
 
  uint256 nonce;
  
  MinerInstance() {}
  bool init(cl_context context, cl_program program, cl_device_id dev, unsigned threadsNum, unsigned threadsPerBlock);
};

class ZCashMiner {
private:
  struct solsPerSecond {
    time_t time;
    unsigned sols;
  };  
  
private:
  unsigned mID;
  unsigned pipelinesNum;
  MinerInstance *miners;
  unsigned _threadsNum;
  unsigned _threadsPerBlocksNum;  

  
  std::list<solsPerSecond> _stats;

  void pushStats(unsigned solsNum) {
    auto currentTime = time(0);
    if (_stats.empty()) {
      solsPerSecond s = {currentTime, solsNum};
      _stats.push_back(s);
    } else if (_stats.back().time == currentTime) {
      _stats.back().sols += solsNum;
    } else {
      solsPerSecond s = {currentTime, solsNum};
      _stats.push_back(s);
    }
  }
  
  void cleanupStats() {
    auto currentTime = time(0);
    while (!_stats.empty() &&_stats.front().time - currentTime >= 60)
      _stats.pop_front();
  }
  
  double calcStats() {
    auto currentTime = time(0);
    auto It = _stats.rbegin();
    auto last = currentTime;
    unsigned sum = 0;
    while (It != _stats.rend() && (It->time - currentTime) < 60) {
      last = It->time;
      sum += It->sols;
      ++It;
    }
    
    auto timeDiff = currentTime - last;
    return timeDiff ? (double)sum / timeDiff : 0;
  }
  

public:
#pragma pack(push, 1)
  struct CBlockHeader {
    // header
    static const size_t HEADER_SIZE=4+32+32+32+4+4+32; // excluding Equihash solution
    static const int32_t CURRENT_VERSION=4;
    
    struct {
      int32_t nVersion;
      uint256 hashPrevBlock;
      uint256 hashMerkleRoot;
      uint256 hashReserved;
      uint32_t nTime;
      uint32_t nBits;
    } data;
    
    uint256 nNonce;    
    
#pragma pack(pop)
    std::vector<unsigned char> nSolution;
  };
  
  bool _cancel;
  
public:
  ZCashMiner(unsigned id);
  
  bool CheckEquihashSolution(const CBlockHeader *pblock, const uint8_t *proof, size_t size);
  
  bool Initialize(cl_context context,
                  cl_program program,
                  cl_device_id dev,
                  unsigned pipelines,
                  unsigned threadsNum,
                  unsigned threadsPerBlock);
  static void InvokeMining(void *args, zctx_t *ctx, void *pipe);
  void Mining(zctx_t *ctx, void *pipe);
  void cancel() { _cancel = true; }
};

class ZCashGPUClient : public BaseClient {
public:
  ZCashGPUClient(zctx_t *ctx) : BaseClient(ctx) {};
  virtual ~ZCashGPUClient();
  
  bool Initialize(Configuration* cfg, bool benchmarkOnly = false);
  void NotifyBlock(const proto::Block& block);
  void TakeWork(const proto::Work& work);
  int GetStats(proto::ClientStats& stats);
  void Toggle();
  void setup_adl();  
  
private:
  std::vector<std::pair<ZCashMiner*, void*> > mWorkers;
};
