/*
 * xpmclient.h
 *
 *  Created on: 01.05.2014
 *      Author: mad
 */

#ifndef XPMCLIENT_H_
#define XPMCLIENT_H_


#include <gmp.h>
#include <gmpxx.h>

#include "baseclient.h"
#include "opencl.h"
#include "uint256.h"
#include "sha256.h"

#define FERMAT_PIPELINES 2

#define PW 128      // Pipeline width (number of hashes to store)
#define SW 5      // number of sieves in one iteration
#define MSO 128*1024    // max sieve output
#define MFS 2*SW*MSO  // max fermat size


extern unsigned gPrimes[96*1024];
extern std::vector<unsigned> gPrimes2;

extern cl_program gProgram;





struct stats_t {
	
	unsigned id;
	unsigned errors;
	unsigned fps;
	double primeprob;
	double cpd;
	
	stats_t(){
		id = 0;
		errors = 0;
		fps = 0;
		primeprob = 0;
		cpd = 0;
	}
	
};


struct config_t {
	
	cl_uint N;
	cl_uint SIZE;
	cl_uint STRIPES;
	cl_uint WIDTH;
	cl_uint PCOUNT;
	cl_uint TARGET;
	
};



class PrimeMiner {
public:
	
	struct block_t {
		
		static const int CURRENT_VERSION = 2;
		
		int version;
		uint256 hashPrevBlock;
		uint256 hashMerkleRoot;
		unsigned int time;
		unsigned int bits;
		unsigned int nonce;
		
	};
	
	struct search_t {
		
		clBuffer<cl_uint> midstate;
		clBuffer<cl_uint> found;
    clBuffer<cl_uint> primorialBitField;
		clBuffer<cl_uint> count;
		
	};
	
	struct hash_t {
		
		unsigned iter;
		unsigned nonce;
		unsigned time;
		uint256 hash;
		mpz_class shash;
    mpz_class primorial;
	};
	
	struct fermat_t {
		
		cl_uint index;
		cl_uchar origin;
		cl_uchar chainpos;
		cl_uchar type;
		cl_uchar hashid;
		
	};
	
	struct info_t {
		
		clBuffer<fermat_t> info;
		clBuffer<cl_uint> count;
		
	};
	
	struct pipeline_t {
		unsigned current;
		unsigned bsize;
		clBuffer<cl_uint> input;
		clBuffer<cl_uchar> output;
		info_t buffer[2];
	};
  
  struct sieve_t {
    info_t cunningham1[1];
    info_t cunningham2[1];
  };
	
	
	PrimeMiner(unsigned id, unsigned threads, unsigned hashprim, unsigned prim, unsigned depth);
	~PrimeMiner();
	
	bool Initialize(cl_device_id dev);
	
	static void InvokeMining(void *args, zctx_t *ctx, void *pipe);
	
	bool MakeExit;
	
private:
  void FermatInit(pipeline_t &fermat, unsigned mfs);
  void FermatDispatch(pipeline_t &fermat,
                      clBuffer<fermat_t>  sieveBuffers[SW][FERMAT_PIPELINES][2],
                      clBuffer<cl_uint> candidatesCountBuffers[SW][2],
                      unsigned pipelineIdx,
                      int ridx,
                      int widx,
                      uint64_t &testCount,
                      uint64_t &fermatCount,
                      cl_kernel fermatKernel);
	void Mining(zctx_t *ctx, void *pipe);
	
	
	unsigned mID;
	unsigned mThreads;
	
	config_t mConfig;
	unsigned mPrimorial;
	unsigned mHashPrimorial;
	unsigned mBlockSize;
	cl_uint mDepth;
	
	cl_command_queue mSmall;
	cl_command_queue mBig;
	
	cl_kernel mHashMod;
	cl_kernel mSieveSetup;
	cl_kernel mSieve;
	cl_kernel mSieveSearch;
	cl_kernel mFermatSetup;
	cl_kernel mFermatKernel352;
  cl_kernel mFermatKernel320;  
	cl_kernel mFermatCheck;
	
  info_t final;	
};




class XPMClient {
public:
	
	XPMClient(zctx_t* ctx);
	~XPMClient();
	
	bool Initialize(Configuration* cfg);
	
	void NotifyBlock(const proto::Block& block);
	
	void TakeWork(const proto::Work& work);
	
	int GetStats(proto::ClientStats& stats);
	
	void Toggle();
	
	void setup_adl();
	
private:
	
	zctx_t* mCtx;
	
	std::vector<std::pair<PrimeMiner*, void*> > mWorkers;
	std::map<int,int> mDeviceMap;
	std::map<int,int> mDeviceMapRev;
	
	void* mBlockPub;
	void* mWorkPub;
	void* mStatsPull;
	
	unsigned mNumDevices;
	unsigned mStatCounter;
	bool mPaused;
	
	std::vector<int> mCoreFreq;
	std::vector<int> mMemFreq;
	std::vector<int> mPowertune;
        std::vector<int> mFanSpeed;        
	
	
	
};



extern XPMClient* gClient;













#endif /* XPMCLIENT_H_ */
