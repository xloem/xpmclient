/*
 * xpmclient.cpp
 *
 *  Created on: 01.05.2014
 *      Author: mad
 */



#include "xpmclient.h"
#include "prime.h"
#include "benchmarks.h"

extern "C" {
	#include "adl.h"
}

#include <fstream>
#include <set>
#include <memory>
#include <chrono>
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (__cplusplus < 201103L)
#define steady_clock monotonic_clock
#endif  

#include <math.h>

cl_platform_id gPlatform = 0;

std::vector<unsigned> gPrimes2;

double GetPrimeDifficulty(unsigned int nBits)
{
    return ((double) nBits / (double) (1 << nFractionalBits));
}

BaseClient *createClient(void *ctx)
{
  return new XPMClient(ctx);
}

PrimeMiner::PrimeMiner(unsigned id, unsigned threads, unsigned sievePerRound, unsigned depth, unsigned LSize) {
	
	mID = id;
	mThreads = threads;

  mSievePerRound = sievePerRound;
	mDepth = depth;
  mLSize = LSize;  
	
	mBlockSize = 0;
	mConfig = {0};
	
  _context = 0;
	mBig = 0;
	mSmall = 0;
	mHashMod = 0;
	mSieveSetup = 0;
	mSieve = 0;
	mSieveSearch = 0;
	mFermatSetup = 0;
	mFermatKernel352 = 0;
  mFermatKernel320 = 0;  
	mFermatCheck = 0;
	
	MakeExit = false;
	
}

PrimeMiner::~PrimeMiner() {
	
  if (_context) OCL(clReleaseContext(_context));
	if(mBig) OCL(clReleaseCommandQueue(mBig));
	if(mSmall) OCL(clReleaseCommandQueue(mSmall));
	
	if(mHashMod) OCL(clReleaseKernel(mHashMod));
	if(mSieveSetup) OCL(clReleaseKernel(mSieveSetup));
	if(mSieve) OCL(clReleaseKernel(mSieve));
	if(mSieveSearch) OCL(clReleaseKernel(mSieveSearch));
	if(mFermatSetup) OCL(clReleaseKernel(mFermatSetup));
	if(mFermatKernel352) OCL(clReleaseKernel(mFermatKernel352));
  if(mFermatKernel320) OCL(clReleaseKernel(mFermatKernel320));  
	if(mFermatCheck) OCL(clReleaseKernel(mFermatCheck));
	
}


bool PrimeMiner::Initialize(cl_context context, cl_program program, cl_device_id dev) {
	
	cl_int error;
  _context = context;
  
  mHashMod = clCreateKernel(program, "bhashmodUsePrecalc", &error);
	mSieveSetup = clCreateKernel(program, "setup_sieve", &error);
	mSieve = clCreateKernel(program, "sieve", &error);
	mSieveSearch = clCreateKernel(program, "s_sieve", &error);
	mFermatSetup = clCreateKernel(program, "setup_fermat", &error);
	mFermatKernel352 = clCreateKernel(program, "fermat_kernel", &error);
  mFermatKernel320 = clCreateKernel(program, "fermat_kernel320", &error);  
	mFermatCheck = clCreateKernel(program, "check_fermat", &error);
	OCLR(error, false);
	
	mBig = clCreateCommandQueue(context, dev, 0, &error);
	mSmall = clCreateCommandQueue(context, dev, 0, &error);
	OCLR(error, false);
	
	{
		clBuffer<config_t> config;
		config.init(context, 1);
		
		cl_kernel getconf = clCreateKernel(program, "getconfig", &error);
		OCLR(error, false);
                size_t globalSize = 1;
                size_t localSize = 1;
		OCLR(clSetKernelArg(getconf, 0, sizeof(cl_mem), &config.DeviceData), false);
		OCLR(clEnqueueNDRangeKernel(mSmall, getconf, 1, 0, &globalSize, &localSize, 0, 0, 0), false);
		config.copyToHost(mSmall, true);
		
		mConfig = *config.HostData;
		
		OCLR(clReleaseKernel(getconf), false);
	}
	
	printf("N=%d SIZE=%d STRIPES=%d WIDTH=%d PCOUNT=%d TARGET=%d\n",
			mConfig.N, mConfig.SIZE, mConfig.STRIPES, mConfig.WIDTH, mConfig.PCOUNT, mConfig.TARGET);
	
	cl_uint numCU;
	OCLR(clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &numCU, 0), false);
	mBlockSize = numCU * 4 * 64;
	printf("GPU %d: has %d CUs\n", mID, numCU);
	
	return true;
	
}

void PrimeMiner::InvokeMining(void *args, void *ctx, void *pipe) {
	
	((PrimeMiner*)args)->Mining(ctx, pipe);
	
}

void PrimeMiner::FermatInit(pipeline_t &fermat, unsigned mfs)
{
  fermat.current = 0;
  fermat.bsize = 0;
  fermat.input.init(_context, mfs*mConfig.N, CL_MEM_HOST_NO_ACCESS);
  fermat.output.init(_context, mfs, CL_MEM_HOST_NO_ACCESS);

  for(int i = 0; i < 2; ++i){
    fermat.buffer[i].info.init(_context, mfs, CL_MEM_HOST_NO_ACCESS);
    fermat.buffer[i].count.init(_context, 1, CL_MEM_ALLOC_HOST_PTR);
  }
}

void PrimeMiner::FermatDispatch(pipeline_t &fermat,
                                clBuffer<fermat_t> sieveBuffers[SW][FERMAT_PIPELINES][2],
                                clBuffer<cl_uint> candidatesCountBuffers[SW][2],
                                unsigned pipelineIdx,
                                int ridx,
                                int widx,
                                uint64_t &testCount,
                                uint64_t &fermatCount,
                                cl_kernel fermatKernel,
                                unsigned sievePerRound)
{ 
  // fermat dispatch
  {
    cl_uint& count = fermat.buffer[ridx].count[0];
    
    cl_uint left = fermat.buffer[widx].count[0] - fermat.bsize;
    if(left > 0){
      OCL(clEnqueueCopyBuffer(  mBig,
                                fermat.buffer[widx].info.DeviceData,
                                fermat.buffer[ridx].info.DeviceData,
                                fermat.bsize*sizeof(fermat_t), count*sizeof(fermat_t),
                                left*sizeof(fermat_t), 0, 0, 0));
      count += left;
    }
    
    for(int i = 0; i < sievePerRound; ++i){
      cl_uint& avail = (candidatesCountBuffers[i][ridx])[pipelineIdx];
      if(avail){
        OCL(clEnqueueCopyBuffer(mBig,
                                sieveBuffers[i][pipelineIdx][ridx].DeviceData,
                                fermat.buffer[ridx].info.DeviceData,
                                0, count*sizeof(fermat_t),
                                avail*sizeof(fermat_t), 0, 0, 0));
        count += avail;
        testCount += avail;
        fermatCount += avail;
        avail = 0;
      }
    }
    
    fermat.buffer[widx].count[0] = 0;
    fermat.buffer[widx].count.copyToDevice(mBig, false);
    
    fermat.bsize = 0;
    if(count > mBlockSize){                 
      fermat.bsize = count - (count % mBlockSize);
      size_t globalSize[] = { fermat.bsize, 1, 1 };
      size_t localSize[] = { 64, 1, 1 };
      OCL(clSetKernelArg(mFermatSetup, 0, sizeof(cl_mem), &fermat.input.DeviceData));      
      OCL(clSetKernelArg(mFermatSetup, 1, sizeof(cl_mem), &fermat.buffer[ridx].info.DeviceData));
      OCL(clEnqueueNDRangeKernel(mBig, mFermatSetup, 1, 0, globalSize, 0, 0, 0, 0));
      OCL(clSetKernelArg(fermatKernel, 0, sizeof(cl_mem), &fermat.output.DeviceData));
      OCL(clSetKernelArg(fermatKernel, 1, sizeof(cl_mem), &fermat.input.DeviceData));      
      OCL(clEnqueueNDRangeKernel(mBig, fermatKernel, 1, 0, globalSize, localSize, 0, 0, 0));
      OCL(clSetKernelArg(mFermatCheck, 0, sizeof(cl_mem), &fermat.buffer[widx].info.DeviceData));
      OCL(clSetKernelArg(mFermatCheck, 1, sizeof(cl_mem), &fermat.buffer[widx].count.DeviceData));
      OCL(clSetKernelArg(mFermatCheck, 4, sizeof(cl_mem), &fermat.output.DeviceData));      
      OCL(clSetKernelArg(mFermatCheck, 5, sizeof(cl_mem), &fermat.buffer[ridx].info.DeviceData));
      OCL(clEnqueueNDRangeKernel(mBig, mFermatCheck, 1, 0, globalSize, 0, 0, 0, 0));
      fermat.buffer[widx].count.copyToHost(mBig, false);
    } else {
//       printf(" * warning: no enough candidates available (pipeline %u)\n", pipelineIdx);
    }
    //printf("fermat: total of %d infos, bsize = %d\n", count, fermat.bsize);
  }
}

void PrimeMiner::Mining(void *ctx, void *pipe) {
	void* blocksub = zmq_socket(ctx, ZMQ_SUB);
	void* worksub = zmq_socket(ctx, ZMQ_SUB);
	void* statspush = zmq_socket(ctx, ZMQ_PUSH);
	void* sharepush = zmq_socket(ctx, ZMQ_PUSH);  
	
	zmq_connect(blocksub, "inproc://blocks");
	zmq_connect(worksub, "inproc://work");
	zmq_connect(statspush, "inproc://stats");
	zmq_connect(sharepush, "inproc://shares");        

	{
		const char one[2] = {1, 0};
		zmq_setsockopt (blocksub, ZMQ_SUBSCRIBE, one, 1);
		zmq_setsockopt (worksub, ZMQ_SUBSCRIBE, one, 1);
	}
	
	proto::Block block;
	proto::Work work;
	proto::Share share;
	
	block.set_height(1);
	work.set_height(0);
	
	share.set_addr(gAddr);
	share.set_name(gClientName);
	share.set_clientid(gClientID);
	
	stats_t stats;
	stats.id = mID;
	stats.errors = 0;
	stats.fps = 0;
	stats.primeprob = 0;
	stats.cpd = 0;
	
	const unsigned mPrimorial = 13;
	uint64_t fermatCount = 1;
	uint64_t primeCount = 1;
	
	time_t time1 = time(0);
	time_t time2 = time(0);
	uint64_t testCount = 0;

	unsigned iteration = 0;
	mpz_class primorial[maxHashPrimorial];
	block_t blockheader;
	search_t hashmod;

  lifoBuffer<hash_t> hashes(PW);
	clBuffer<cl_uint> hashBuf;
	clBuffer<cl_uint> sieveBuf[2];
	clBuffer<cl_uint> sieveOff[2];
  clBuffer<fermat_t> sieveBuffers[SW][FERMAT_PIPELINES][2];
  clBuffer<cl_uint> candidatesCountBuffers[SW][2];
	pipeline_t fermat320;
  pipeline_t fermat352;
	CPrimalityTestParams testParams;
	std::vector<fermat_t> candis;
  unsigned numHashCoeff = 32768;
	
  cl_mem primeBuf[maxHashPrimorial];
  cl_mem primeBuf2[maxHashPrimorial];
  
  for (unsigned i = 0; i < maxHashPrimorial - mPrimorial; i++) {
    cl_int error = 0;
    primeBuf[i] = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                 mConfig.PCOUNT*sizeof(cl_uint), &gPrimes[mPrimorial+i+1], &error);
    OCL(error);
    primeBuf2[i] = clCreateBuffer(_context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  mConfig.PCOUNT*2*sizeof(cl_uint), &gPrimes2[2*(mPrimorial+i)+2], &error);
    OCL(error);
    
    mpz_class p = 1;
    for(unsigned j = 0; j <= mPrimorial+i; j++)
      p *= gPrimes[j];
    
    primorial[i] = p;
  }
  
	{
		unsigned primorialbits = mpz_sizeinbase(primorial[0].get_mpz_t(), 2);
		mpz_class sievesize = mConfig.SIZE*32*mConfig.STRIPES;
		unsigned sievebits = mpz_sizeinbase(sievesize.get_mpz_t(), 2);
		printf("GPU %d: primorial = %s (%d bits)\n", mID, primorial[0].get_str(10).c_str(), primorialbits);
		printf("GPU %d: sieve size = %s (%d bits)\n", mID, sievesize.get_str(10).c_str(), sievebits);
	}
	
	hashmod.midstate.init(_context, 8*sizeof(cl_uint), CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY);
	hashmod.found.init(_context, 128, CL_MEM_ALLOC_HOST_PTR);
  hashmod.primorialBitField.init(_context, 128, CL_MEM_ALLOC_HOST_PTR);
	hashmod.count.init(_context, 1, CL_MEM_ALLOC_HOST_PTR);
	hashBuf.init(_context, PW*mConfig.N, CL_MEM_READ_WRITE);
	
	for(int sieveIdx = 0; sieveIdx < SW; ++sieveIdx) {
    for(int instIdx = 0; instIdx < 2; ++instIdx){    
      for (int pipelineIdx = 0; pipelineIdx < FERMAT_PIPELINES; pipelineIdx++)
        sieveBuffers[sieveIdx][pipelineIdx][instIdx].init(_context, MSO, CL_MEM_HOST_NO_ACCESS);
      
      candidatesCountBuffers[sieveIdx][instIdx].init(_context, FERMAT_PIPELINES, CL_MEM_ALLOC_HOST_PTR);
    }
  }
	
	for(int k = 0; k < 2; ++k){
		sieveBuf[k].init(_context, mConfig.SIZE*mConfig.STRIPES/2*mConfig.WIDTH, CL_MEM_HOST_NO_ACCESS);
		sieveOff[k].init(_context, mConfig.PCOUNT*mConfig.WIDTH, CL_MEM_HOST_NO_ACCESS);
	}
	
	final.info.init(_context, MFS/(4*mDepth), CL_MEM_ALLOC_HOST_PTR);
  final.count.init(_context, 1, CL_MEM_ALLOC_HOST_PTR);	
	
	FermatInit(fermat320, MFS);
  FermatInit(fermat352, MFS);  

  clBuffer<cl_uint> modulosBuf[maxHashPrimorial];
  unsigned modulosBufferSize = mConfig.PCOUNT*(mConfig.N-1);   
  for (unsigned bufIdx = 0; bufIdx < maxHashPrimorial-mPrimorial; bufIdx++) {
    clBuffer<cl_uint> &current = modulosBuf[bufIdx];
    current.init(_context, modulosBufferSize, CL_MEM_READ_ONLY);
    for (unsigned i = 0; i < mConfig.PCOUNT; i++) {
      mpz_class X = 1;
      for (unsigned j = 0; j < mConfig.N-1; j++) {
        X <<= 32;
        mpz_class mod = X % gPrimes[i+mPrimorial+bufIdx+1];
        current[mConfig.PCOUNT*j+i] = mod.get_ui();
      }
    }
    
    current.copyToDevice(mSmall);
  }

	OCL(clSetKernelArg(mHashMod, 0, sizeof(cl_mem), &hashmod.found.DeviceData));
	OCL(clSetKernelArg(mHashMod, 1, sizeof(cl_mem), &hashmod.count.DeviceData));
	OCL(clSetKernelArg(mHashMod, 2, sizeof(cl_mem), &hashmod.primorialBitField.DeviceData));
	OCL(clSetKernelArg(mHashMod, 3, sizeof(cl_mem), &hashmod.midstate.DeviceData));
	OCL(clSetKernelArg(mSieveSetup, 0, sizeof(cl_mem), &sieveOff[0].DeviceData));
	OCL(clSetKernelArg(mSieveSetup, 1, sizeof(cl_mem), &sieveOff[1].DeviceData));
	OCL(clSetKernelArg(mSieveSetup, 3, sizeof(cl_mem), &hashBuf.DeviceData));
	OCL(clSetKernelArg(mSieveSearch, 0, sizeof(cl_mem), &sieveBuf[0].DeviceData));
	OCL(clSetKernelArg(mSieveSearch, 1, sizeof(cl_mem), &sieveBuf[1].DeviceData));
  OCL(clSetKernelArg(mSieveSearch, 7, sizeof(cl_uint), &mDepth));  
	OCL(clSetKernelArg(mFermatSetup, 2, sizeof(cl_mem), &hashBuf.DeviceData));
	OCL(clSetKernelArg(mFermatCheck, 2, sizeof(cl_mem), &final.info.DeviceData));
	OCL(clSetKernelArg(mFermatCheck, 3, sizeof(cl_mem), &final.count.DeviceData));
	OCL(clSetKernelArg(mFermatCheck, 6, sizeof(unsigned), &mDepth));

	czmq_signal(pipe);
	czmq_poll(pipe, -1);
	
	bool run = true;
	while(run){
		if(czmq_poll(pipe, 0)){
			czmq_wait(pipe);
			czmq_wait(pipe);
		}
		
		//printf("\n--------- iteration %d -------\n", iteration);
		{
			time_t currtime = time(0);
			time_t elapsed = currtime - time1;
			if(elapsed > 11){
 				zmq_send(statspush, &stats, sizeof(stats), 0);                          
				time1 = currtime;
			}
			
			elapsed = currtime - time2;
			if(elapsed > 15){
				stats.fps = testCount / elapsed;
				time2 = currtime;
				testCount = 0;
			}
		}
		
		stats.primeprob = pow(double(primeCount)/double(fermatCount), 1./mDepth)
				- 0.0003 * (double(mConfig.TARGET-1)/2. - double(mDepth-1)/2.);
		stats.cpd = 24.*3600. * double(stats.fps) * pow(stats.primeprob, mConfig.TARGET);
		
		// get work
		bool reset = false;
		{
			bool getwork = true;
			while(getwork && run){
				if(czmq_poll(worksub, 0) || work.height() < block.height()){
					run = ReceivePub(work, worksub);
					reset = true;
				}
				
				getwork = false;
				if(czmq_poll(blocksub, 0) || work.height() > block.height()){
					run = ReceivePub(block, blocksub);
					getwork = true;
				}
			}
		}
		if(!run)
			break;
		
		// reset if new work
		if(reset){
      hashes.clear();
			hashmod.count[0] = 0;
			fermat320.bsize = 0;
			fermat320.buffer[0].count[0] = 0;
			fermat320.buffer[1].count[0] = 0;
      fermat352.bsize = 0;
      fermat352.buffer[0].count[0] = 0;
      fermat352.buffer[1].count[0] = 0;      
			final.count[0] = 0;
      
      for(int sieveIdx = 0; sieveIdx < mSievePerRound; ++sieveIdx) {
        for(int instIdx = 0; instIdx < 2; ++instIdx) {
          for (int pipelineIdx = 0; pipelineIdx < FERMAT_PIPELINES; pipelineIdx++)
            (candidatesCountBuffers[sieveIdx][instIdx])[pipelineIdx] = 0;
        }
      }      
			
			blockheader.version = block_t::CURRENT_VERSION;
			blockheader.hashPrevBlock.SetHex(block.hash());
			blockheader.hashMerkleRoot.SetHex(work.merkle());
			blockheader.time = work.time() + mID;
			blockheader.bits = work.bits();
			blockheader.nonce = 1;
			testParams.nBits = blockheader.bits;
			
			unsigned target = TargetGetLength(blockheader.bits);
// 			if(target > mConfig.TARGET){
// 				printf("ERROR: This miner is compiled with the wrong target: %d (required target %d)\n", mConfig.TARGET, target);
// 				return;
// 			}

			simplePrecalcSHA256(&blockheader, hashmod.midstate, mBig, mHashMod);
		}
		
		// hashmod fetch & dispatch
		{
// 			printf("got %d new hashes\n", hashmod.count[0]);
			for(unsigned i = 0; i < hashmod.count[0]; ++i) {
				hash_t hash;
				hash.iter = iteration;
				hash.time = blockheader.time;
				hash.nonce = hashmod.found[i];
        uint32_t primorialBitField = hashmod.primorialBitField[i];
        uint32_t primorialIdx = primorialBitField >> 16;
        uint64_t realPrimorial = 1;
        for (unsigned j = 0; j < primorialIdx+1; j++) {
          if (primorialBitField & (1 << j))
            realPrimorial *= gPrimes[j];
        }      
        
        mpz_class mpzRealPrimorial;        
        mpz_import(mpzRealPrimorial.get_mpz_t(), 2, -1, 4, 0, 0, &realPrimorial);            
        primorialIdx = std::max(mPrimorial, primorialIdx) - mPrimorial;
        mpz_class mpzHashMultiplier = primorial[primorialIdx] / mpzRealPrimorial;
        unsigned hashMultiplierSize = mpz_sizeinbase(mpzHashMultiplier.get_mpz_t(), 2);      
        mpz_import(mpzRealPrimorial.get_mpz_t(), 2, -1, 4, 0, 0, &realPrimorial);        
				
				block_t b = blockheader;
				b.nonce = hash.nonce;
				
				SHA_256 sha;
				sha.init();
				sha.update((const unsigned char*)&b, sizeof(b));
				sha.final((unsigned char*)&hash.hash);
				sha.init();
				sha.update((const unsigned char*)&hash.hash, sizeof(uint256));
				sha.final((unsigned char*)&hash.hash);
				
				if(hash.hash < (uint256(1) << 255)){
					printf("error: hash does not meet minimum.\n");
					stats.errors++;
					continue;
				}
				
				mpz_class mpzHash;
				mpz_set_uint256(mpzHash.get_mpz_t(), hash.hash);
        if(!mpz_divisible_p(mpzHash.get_mpz_t(), mpzRealPrimorial.get_mpz_t())){
					printf("error: mpz_divisible_ui_p failed.\n");
					stats.errors++;
					continue;
				}
				
				hash.primorialIdx = primorialIdx;
        hash.primorial = mpzHashMultiplier;
        hash.shash = mpzHash * hash.primorial;       

        unsigned hid = hashes.push(hash);
        memset(&hashBuf[hid*mConfig.N], 0, sizeof(uint32_t)*mConfig.N);
        mpz_export(&hashBuf[hid*mConfig.N], 0, -1, 4, 0, 0, hashes.get(hid).shash.get_mpz_t());        
			}
			
			if (hashmod.count[0])
        hashBuf.copyToDevice(mSmall, false);
			
			//printf("hashlist.size() = %d\n", (int)hashlist.size());
			hashmod.count[0] = 0;
			
      int numhash = ((int)(16*mSievePerRound) - (int)hashes.remaining()) * numHashCoeff;

			if(numhash > 0){
        numhash += mLSize - numhash % mLSize;
				if(blockheader.nonce > (1u << 31)){
					blockheader.time += mThreads;
					blockheader.nonce = 1;
          simplePrecalcSHA256(&blockheader, hashmod.midstate, mBig, mHashMod);
				}

				size_t globalOffset[] = { blockheader.nonce, 1u, 1u };
				size_t globalSize[] = { (unsigned)numhash, 1u, 1u };
        size_t localSize[] = { mLSize, 1 };      
				hashmod.count.copyToDevice(mBig, false);
        OCL(clEnqueueNDRangeKernel(mBig, mHashMod, 1, globalOffset, globalSize, localSize, 0, 0, 0));
				hashmod.found.copyToHost(mBig, false);
        hashmod.primorialBitField.copyToHost(mBig, false);
				hashmod.count.copyToHost(mBig, false);
				blockheader.nonce += numhash;
			}
			
		}

		int ridx = iteration % 2;
		int widx = ridx xor 1;
		
		// sieve dispatch    
      for (unsigned i = 0; i < mSievePerRound; i++) {
        if(hashes.empty()){
          if (!reset) {
            numHashCoeff += 32768;
            printf(" * warning: ran out of hashes, increasing sha256 work size coefficient to %u\n", numHashCoeff);
          }
          break;
        }
        
        cl_int hid = hashes.pop();
        unsigned primorialIdx = hashes.get(hid).primorialIdx;
        OCL(clSetKernelArg(mSieveSetup, 2, sizeof(cl_mem), &primeBuf[primorialIdx]));
        OCL(clSetKernelArg(mSieveSetup, 5, sizeof(cl_mem), &modulosBuf[primorialIdx].DeviceData));          
        OCL(clSetKernelArg(mSieve, 2, sizeof(cl_mem), &primeBuf2[primorialIdx]));        

        {
					OCL(clSetKernelArg(mSieveSetup, 4, sizeof(cl_int), &hid));
					size_t globalSize[] = { mConfig.PCOUNT, 1, 1 };
					OCL(clEnqueueNDRangeKernel(mSmall, mSieveSetup, 1, 0, globalSize, 0, 0, 0, 0));
				}

        {
          OCL(clSetKernelArg(mSieve, 0, sizeof(cl_mem), &sieveBuf[0].DeviceData));
          OCL(clSetKernelArg(mSieve, 1, sizeof(cl_mem), &sieveOff[0].DeviceData));
          size_t globalSize[] = { mLSize*mConfig.STRIPES/2, mConfig.WIDTH};
          size_t localSize[] = { mLSize, 1 };          
          OCL(clEnqueueNDRangeKernel(mSmall, mSieve, 2, 0, globalSize, localSize, 0, 0, 0));
        }
        
        {
          OCL(clSetKernelArg(mSieve, 0, sizeof(cl_mem), &sieveBuf[1].DeviceData));
          OCL(clSetKernelArg(mSieve, 1, sizeof(cl_mem), &sieveOff[1].DeviceData));              
          size_t globalSize[] = { mLSize*mConfig.STRIPES/2, mConfig.WIDTH};
          size_t localSize[] = { mLSize, 1 };          
          OCL(clEnqueueNDRangeKernel(mSmall, mSieve, 2, 0, globalSize, localSize, 0, 0, 0));
        }         

				candidatesCountBuffers[i][widx].copyToDevice(mSmall, false);
         
				{
          cl_uint multiplierSize = mpz_sizeinbase(hashes.get(hid).shash.get_mpz_t(), 2);
          OCL(clSetKernelArg(mSieveSearch, 2, sizeof(cl_mem), &sieveBuffers[i][0][widx].DeviceData));
          OCL(clSetKernelArg(mSieveSearch, 3, sizeof(cl_mem), &sieveBuffers[i][1][widx].DeviceData));          
          OCL(clSetKernelArg(mSieveSearch, 4, sizeof(cl_mem), &candidatesCountBuffers[i][widx].DeviceData));
					OCL(clSetKernelArg(mSieveSearch, 5, sizeof(cl_int), &hid));
          OCL(clSetKernelArg(mSieveSearch, 6, sizeof(cl_uint), &multiplierSize));
					size_t globalSize[] = { mConfig.SIZE*mConfig.STRIPES/2, 1, 1 };
					OCL(clEnqueueNDRangeKernel(mSmall, mSieveSearch, 1, 0, globalSize, 0, 0, 0, 0));
          
          candidatesCountBuffers[i][widx].copyToHost(mSmall, false);
				}
			}
		
    
		// get candis
		int numcandis = final.count[0];
		numcandis = std::min(numcandis, final.info.Size);
		numcandis = std::max(numcandis, 0);
// 		printf("got %d new candis\n", numcandis);
		candis.resize(numcandis);
		primeCount += numcandis;
		if(numcandis)
			memcpy(&candis[0], final.info.HostData, numcandis*sizeof(fermat_t));
		
    final.count[0] = 0;
    final.count.copyToDevice(mBig, false);    
    FermatDispatch(fermat320, sieveBuffers, candidatesCountBuffers, 0, ridx, widx, testCount, fermatCount, mFermatKernel320, mSievePerRound);    
    FermatDispatch(fermat352, sieveBuffers, candidatesCountBuffers, 1, ridx, widx, testCount, fermatCount, mFermatKernel352, mSievePerRound);
    final.info.copyToHost(mBig, false);
    final.count.copyToHost(mBig, false);         
		
		clFlush(mBig);
		clFlush(mSmall);
    
    // adjust sieves per round
    if (fermat320.buffer[ridx].count[0] && fermat320.buffer[ridx].count[0] < mBlockSize &&
        fermat352.buffer[ridx].count[0] && fermat352.buffer[ridx].count[0] < mBlockSize) {
      mSievePerRound = std::min((unsigned)SW, mSievePerRound+1);
      printf(" * warning: not enough candidates (%u available, must be more than %u\n",
             std::max(fermat320.buffer[ridx].count[0], fermat352.buffer[ridx].count[0]),
             mBlockSize);
             
      printf("   * increase sieves per round to %u\n", mSievePerRound);
    }
		
		// check candis
		if(candis.size()){
			//printf("checking %d candis\n", (int)candis.size());
			mpz_class chainorg;
			mpz_class multi;
			for(unsigned i = 0; i < candis.size(); ++i){
				
				fermat_t& candi = candis[i];
				hash_t& hash = hashes.get(candi.hashid);
				
				unsigned age = iteration - hash.iter;
				if(age > PW/2)
					printf("WARNING: candidate age > PW/2 with %d\n", age);
				
				multi = candi.index;
				multi <<= candi.origin;
				chainorg = hash.shash;
				chainorg *= multi;
				
				testParams.nCandidateType = candi.type;
        bool isblock = ProbablePrimeChainTestFast(chainorg, testParams, mDepth);
				unsigned chainlength = TargetGetLength(testParams.nChainLength);

				/*printf("candi %d: hashid=%d index=%d origin=%d type=%d length=%d\n",
						i, candi.hashid, candi.index, candi.origin, candi.type, chainlength);*/
				if(chainlength >= block.minshare()){
					
					mpz_class sharemulti = hash.primorial * multi;
					share.set_hash(hash.hash.GetHex());
					share.set_merkle(work.merkle());
					share.set_time(hash.time);
					share.set_bits(work.bits());
					share.set_nonce(hash.nonce);
					share.set_multi(sharemulti.get_str(16));
					share.set_height(block.height());
					share.set_length(chainlength);
					share.set_chaintype(candi.type);
					share.set_isblock(isblock);
					
					printf("GPU %d found share: %d-ch type %d\n", mID, chainlength, candi.type+1);
					if(isblock)
						printf("GPU %d found BLOCK!\n", mID);
					
					Send(share, sharepush);
					
				}else if(chainlength < mDepth){
					printf("error: ProbablePrimeChainTestFast %ubits %d/%d\n", (unsigned)mpz_sizeinbase(chainorg.get_mpz_t(), 2), chainlength, mDepth);
          printf(" * origin: %s\n", chainorg.get_str().c_str());
          printf(" * type: %u\n", (unsigned)candi.type);
          printf(" * multiplier: %u\n", (unsigned)candi.index);
          printf(" * layer: %u\n", (unsigned)candi.origin);
          printf(" * hash primorial: %s\n", hash.primorial.get_str().c_str());
          printf("   * primorial multipliers: ");
          for (unsigned i = 0; i < mPrimorial;) {
            if (hash.primorial % gPrimes[i] == 0) {
              hash.primorial /= gPrimes[i];
              printf("[%u]%u ", i+1, gPrimes[i]);
            } else {
              i++;
            }
          }
          printf("\n");
					stats.errors++;
				}
			}
		}
		
		clFinish(mBig);
		clFinish(mSmall);
		
		if(MakeExit)
			break;
		
		iteration++;
	}
	
	printf("GPU %d stopped.\n", mID);
	
  for (unsigned i = 0; i < maxHashPrimorial-mPrimorial; i++) {
	  clReleaseMemObject(primeBuf[i]);
	  clReleaseMemObject(primeBuf2[i]);
  }
	
	zmq_close(blocksub);
	zmq_close(worksub);
	zmq_close(statspush);
	zmq_close(sharepush);
	czmq_signal(pipe);
}

XPMClient::~XPMClient() {
	
	for(unsigned i = 0; i < mWorkers.size(); ++i)
		if(mWorkers[i].first){
			mWorkers[i].first->MakeExit = true;
			if(czmq_poll(mWorkers[i].second, 1000))
				delete mWorkers[i].first;
		}

	zmq_close(mBlockPub);
	zmq_close(mWorkPub);
	zmq_close(mStatsPull);
  if (platformType == ptAMD)
  	clear_adl(mNumDevices);
	
}

void XPMClient::dumpSieveConstants(unsigned weaveDepth,
                                   unsigned threadsNum,
                                   unsigned windowSize,
                                   unsigned *primes,
                                   std::ostream &file) 
{
  unsigned ranges[3] = {0, 0, 0};
  for (unsigned i = 0; i < weaveDepth/threadsNum; i++) {
    unsigned prime = primes[i*threadsNum];
    if (ranges[0] == 0 && windowSize/prime <= 2)
      ranges[0] = i;
    if (ranges[1] == 0 && windowSize/prime <= 1)
      ranges[1] = i;
    if (ranges[2] == 0 && windowSize/prime == 0)
      ranges[2] = i;
  }
  
  file << "__constant uint sieveRanges[3] = {";
  file << ranges[0] << ", ";
  file << ranges[1] << ", " ;
  file << ranges[2];
  file << "};\n";  
}

bool XPMClient::Initialize(Configuration* cfg, bool benchmarkOnly, unsigned adjustedKernelTarget) {
  cl_context gContext[64] = {0};
  cl_program gProgram[64] = {0};
  
	_cfg = cfg;
	
	{
		int np = sizeof(gPrimes)/sizeof(unsigned);
		gPrimes2.resize(np*2);
		for(int i = 0; i < np; ++i){
			unsigned prime = gPrimes[i];
			cl_float fiprime = 1.f / cl_float(prime);
			gPrimes2[i*2] = prime;
			memcpy(&gPrimes2[i*2+1], &fiprime, sizeof(cl_float));
		}
	}
	
	const char *platformId = cfg->lookupString("", "platform");
	const char *platformName = "";
	unsigned clKernelLSize = 0;
	unsigned clKernelLSizeLog2 = 0;
	DeviceTypeTy deviceType = dtUnknown;
  
  if (strcmp(platformId, "amd") == 0) {
    platformName = "AMD Accelerated Parallel Processing";
    platformType = ptAMD;
    clKernelLSize = 256;
    clKernelLSizeLog2 = 8;
		deviceType = dtAMDGCN;
  } else if (strcmp(platformId, "amd legacy") == 0) {
    platformName = "AMD Accelerated Parallel Processing";
    platformType = ptAMD;
    deviceType = dtAMDLegacy;
    clKernelLSize = 256;
    clKernelLSizeLog2 = 8;
  } else if (strcmp(platformId, "amd vega") == 0) {
    platformName = "AMD Accelerated Parallel Processing";
    platformType = ptAMD;
    deviceType = dtAMDVega;
    clKernelLSize = 256;
    clKernelLSizeLog2 = 8;
  }  else if (strcmp(platformId, "nvidia") == 0) {
    platformName = "NVIDIA CUDA";
    platformType = ptNVidia;    
		deviceType = dtNVIDIA;
    clKernelLSize = 1024;
    clKernelLSizeLog2 = 10;
  }

  std::vector<cl_device_id> allgpus;
  if (!clInitialize(platformName, allgpus))
    return false;
  mNumDevices = allgpus.size();
  
	int cpuload = cfg->lookupInt("", "cpuload", 1);
	int depth = 5 - cpuload;
	depth = std::max(depth, 2);
	depth = std::min(depth, 5);
	
  exitType = cfg->lookupInt("", "onCrash", 0);

	unsigned clKernelTarget = adjustedKernelTarget ? adjustedKernelTarget : 10;
	const char *targetValue = cfg->lookupString("", "target", "auto");
	if (strcmp(targetValue, "auto") != 0) {
		clKernelTarget = atoi(targetValue);
		clKernelTargetAutoAdjust = false;
	}

	bool clKernelWidthAutoAdjust = true;
	unsigned clKernelWidth = clKernelTarget*2;
	const char *widthValue = cfg->lookupString("", "width", "auto");
	if (strcmp(widthValue, "auto") != 0) {
		clKernelWidthAutoAdjust = false;
		clKernelWidth = atoi(widthValue);
	}
	
  unsigned clKernelStripes = cfg->lookupInt("", "sieveSize", 420);
  unsigned clKernelPCount = cfg->lookupInt("", "weaveDepth", 40960);
  unsigned clKernelWindowSize = cfg->lookupInt("", "windowSize", 4096);

	unsigned multiplierSizeLimits[3] = {26, 33, 36};
	std::vector<bool> usegpu(mNumDevices, true);
  std::vector<int> sievePerRound(mNumDevices, 5);
	mCoreFreq = std::vector<int>(mNumDevices, -1);
	mMemFreq = std::vector<int>(mNumDevices, -1);
	mPowertune = std::vector<int>(mNumDevices, 42);
  mFanSpeed = std::vector<int>(mNumDevices, 70);
	
	{
		StringVector cmultiplierlimits;
		StringVector cdevices;
		StringVector csieveperround;
		StringVector ccorespeed;
		StringVector cmemspeed;
		StringVector cpowertune;
		StringVector cfanspeed;
    
		try {
			cfg->lookupList("", "devices", cdevices);
			cfg->lookupList("", "corefreq", ccorespeed);
			cfg->lookupList("", "memfreq", cmemspeed);
			cfg->lookupList("", "powertune", cpowertune);
      cfg->lookupList("", "fanspeed", cfanspeed);
			cfg->lookupList("", "multiplierLimits", cmultiplierlimits);
      cfg->lookupList("", "sievePerRound", csieveperround);
		}catch(const ConfigurationException& ex) {}

		if (cmultiplierlimits.length() == 3) {
			multiplierSizeLimits[0] = atoi(cmultiplierlimits[0]);
			multiplierSizeLimits[1] = atoi(cmultiplierlimits[1]);
			multiplierSizeLimits[2] = atoi(cmultiplierlimits[2]);
		} else {
			printf(" * Warning: invalid multiplierLimits parameter in config, must be list of 3 numbers\n");
		}
		
		for(int i = 0; i < (int)mNumDevices; ++i){
			
			if(i < cdevices.length())
				usegpu[i] = !strcmp(cdevices[i], "1");
			
			if (i < csieveperround.length())
				sievePerRound[i] = atoi(csieveperround[i]);
			
			if(i < ccorespeed.length())
				mCoreFreq[i] = atoi(ccorespeed[i]);
			
			if(i < cmemspeed.length())
				mMemFreq[i] = atoi(cmemspeed[i]);
			
			if(i < cpowertune.length())
				mPowertune[i] = atoi(cpowertune[i]);
			
      if(i < cfanspeed.length())
        mFanSpeed[i] = atoi(cfanspeed[i]);                        
		}
	}
	
	std::vector<cl_device_id> gpus;
	for(unsigned i = 0; i < mNumDevices; ++i)
		if(usegpu[i]){
			printf("Using device %d as GPU %d\n", i, (int)gpus.size());
			mDeviceMap[i] = gpus.size();
			mDeviceMapRev[gpus.size()] = i;
			gpus.push_back(allgpus[i]);
		}else{
			mDeviceMap[i] = -1;
		}
	
	if(!gpus.size()){
		printf("EXIT: config.txt says not to use any devices!?\n");
		return false;
	};
	
	for (size_t i = 0; i < gpus.size(); i++) {
		cl_context_properties props[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)gPlatform, 0 };
		cl_int error;
		gContext[i] = clCreateContext(props, 1, &gpus[i], 0, 0, &error);
		OCLR(error, false);
	}
	
	// generate kernel configuration file
  {
    std::ofstream config("xpm/gpu/config.cl", std::fstream::trunc);
    config << "#define STRIPES " << clKernelStripes << '\n';
    config << "#define WIDTH " << clKernelWidth << '\n';
    config << "#define PCOUNT " << clKernelPCount << '\n';
    config << "#define TARGET " << clKernelTarget << '\n';
    config << "#define SIZE " << clKernelWindowSize << '\n';
    config << "#define LSIZE " << clKernelLSize << '\n';
    config << "#define LSIZELOG2 " << clKernelLSizeLog2 << '\n';
    config << "#define LIMIT13 " << multiplierSizeLimits[0] << '\n';
    config << "#define LIMIT14 " << multiplierSizeLimits[1] << '\n';
    config << "#define LIMIT15 " << multiplierSizeLimits[2] << '\n';    
    dumpSieveConstants(clKernelPCount, clKernelLSize, clKernelWindowSize*32, gPrimes+13, config);
  }

  const char *arguments = "";
  if (platformType == ptAMD) {
		if (deviceType == dtAMDLegacy)
      arguments = "-D__AMDLEGACY";
		else if (deviceType == dtAMDVega)
			arguments = "-D__AMDVEGA";
	} else if (platformType == ptNVidia) {
    arguments = "-D__NVIDIA -cl-nv-verbose";
	}
  
  std::vector<cl_int> binstatus;
  binstatus.resize(gpus.size());	
  for (size_t i = 0; i < gpus.size(); i++) {
    char kernelName[64];
    sprintf(kernelName, "kernelxpm_gpu%u.cl", (unsigned)i);

    // force rebuild kernel if adjusted target received
    if (!clCompileKernel(gContext[i],
                         gpus[i],
                         kernelName,
                         { "xpm/gpu/config.cl", "xpm/gpu/procs.cl", "xpm/gpu/fermat.cl", "xpm/gpu/sieve.cl", "xpm/gpu/sha256.cl", "xpm/gpu/benchmarks.cl"},
                         arguments,
                         &binstatus[i],
                         &gProgram[i],
                         adjustedKernelTarget != 0)) {
      return false;
    }
  }

  if (platformType == ptAMD)
    setup_adl();
  
  if (benchmarkOnly) {
    for (unsigned i = 0; i < gpus.size(); i++) {
      if (binstatus[i] == CL_SUCCESS) {
        runBenchmarks(gContext[0], gProgram[0], gpus[i], depth, clKernelLSize);
      }
    }
    
    return false;
  } else {
    for(unsigned i = 0; i < gpus.size(); ++i) {
      std::pair<PrimeMiner*,void*> worker;
      if(binstatus[i] == CL_SUCCESS){
      
        PrimeMiner* miner = new PrimeMiner(i, gpus.size(), sievePerRound[i], depth, clKernelLSize);
        miner->Initialize(gContext[i], gProgram[i], gpus[i]);
        config_t config = miner->getConfig();
				if ((!clKernelTargetAutoAdjust && config.TARGET != clKernelTarget) ||
						(!clKernelWidthAutoAdjust && config.WIDTH != clKernelWidth) ||
						config.PCOUNT != clKernelPCount ||
						config.STRIPES != clKernelStripes ||
						config.SIZE != clKernelWindowSize ||
						config.LIMIT13 != multiplierSizeLimits[0] ||
						config.LIMIT14 != multiplierSizeLimits[1] ||
						config.LIMIT15 != multiplierSizeLimits[2]) {
          printf("Existing OpenCL kernel (kernel.bin) incompatible with configuration\n");
          printf("Please remove kernel.bin file and restart miner\n");
          exit(1);
        }

        void *pipe = czmq_thread_fork(mCtx, &PrimeMiner::InvokeMining, miner);
        czmq_wait(pipe);
        czmq_signal(pipe);
        worker.first = miner;
        worker.second = pipe;
      } else {
        printf("GPU %d: failed to load kernel\n", i);
        worker.first = 0;
        worker.second = 0;
      
      }
    
      mWorkers.push_back(worker);
    }
  }
	  
	return true;
}


void XPMClient::NotifyBlock(const proto::Block& block) {
	
	SendPub(block, mBlockPub);
	
}


bool XPMClient::TakeWork(const proto::Work& work) {
	
	const double TargetIncrease = 0.994;
	const double TargetDecrease = 0.0061;
	
	SendPub(work, mWorkPub);
	
	if (!clKernelTargetAutoAdjust || work.bits() == 0)
		return true;
	double difficulty = GetPrimeDifficulty(work.bits());

	bool needReset = false;
	for(unsigned i = 0; i < mWorkers.size(); ++i) {
		PrimeMiner *miner = mWorkers[i].first;
		double target = miner->getConfig().TARGET;
		if (difficulty > target && difficulty-target >= TargetIncrease) {
			printf("Target with high difficulty detected, need increase miner target\n");
			needReset = true;
			break;
		} else if (difficulty < target && target-difficulty >= TargetDecrease) {
			printf("Target with low difficulty detected, need decrease miner target\n");
			needReset = true;
			break;
		}
	}
	
	if (needReset) {
		unsigned newTarget = TargetGetLength(work.bits());
		if (difficulty - newTarget >= TargetIncrease)
			newTarget++;
		printf("Rebuild miner kernels, adjust target to %u..\n", newTarget);
		// Stop and destroy all workers
		for(unsigned i = 0; i < mWorkers.size(); ++i) {
			printf("attempt to stop GPU %u ...\n", i);
			if(mWorkers[i].first){
				mWorkers[i].first->MakeExit = true;
				if(czmq_poll(mWorkers[i].second, 8000))
					delete mWorkers[i].first;
			}
		}
		
		mWorkers.clear();
		
		// Build new kernels with adjusted target
		mPaused = true;
		Initialize(_cfg, false, newTarget);
		Toggle();
		return false;
	} else {
		return true;
	}
}


int XPMClient::GetStats(proto::ClientStats& stats) {
	
	unsigned nw = mWorkers.size();
	std::vector<bool> running(nw);
	std::vector<stats_t> wstats(nw);
	
	while (czmq_poll(mStatsPull, 0)) {
		zmq_msg_t msg;
		zmq_msg_init(&msg);
		zmq_recvmsg(mStatsPull, &msg, 0);          
		size_t fsize = zmq_msg_size(&msg);
		uint8_t *fbytes = (uint8_t*)zmq_msg_data(&msg);
		if(fsize >= sizeof(stats_t)) {
			stats_t* tmp = (stats_t*)fbytes;
			if(tmp->id < nw){
				running[tmp->id] = true;
				wstats[tmp->id] = *tmp;
			}
		}
		
		zmq_msg_close(&msg); 
	}

	double cpd = 0;
	unsigned errors = 0;
	int maxtemp = 0;
	unsigned ngpus = 0;
  int crashed = 0;
  
	for(unsigned i = 0; i < nw; ++i){
		
		int devid = mDeviceMapRev[i];
		int temp = gpu_temp(devid);
		int activity = gpu_activity(devid);
		
		if(temp > maxtemp)
			maxtemp = temp;
		
		cpd += wstats[i].cpd;
		errors += wstats[i].errors;
		
		if(running[i]){
			ngpus++;
			printf("[GPU %d] T=%dC A=%d%% E=%d primes=%f fermat=%d/sec cpd=%.2f/day\n",
					i, temp, activity, wstats[i].errors, wstats[i].primeprob, wstats[i].fps, wstats[i].cpd);
		}else if(!mWorkers[i].first)
			printf("[GPU %d] failed to start!\n", i);
		else if(mPaused) {
			printf("[GPU %d] paused\n", i);
    } else {
      crashed++;
			printf("[GPU %d] crashed!\n", i);
    }
		
	}
	
  if (crashed) {
    if (exitType == 1) {
      exit(1);
    } else if (exitType == 2) {
#ifdef WIN32
      ExitWindowsEx(EWX_REBOOT, 0);
#else
      system("/sbin/reboot");
#endif
    }
  }
	
	if(mStatCounter % 10 == 0)
		for(unsigned i = 0; i < mNumDevices; ++i){
			int gpuid = mDeviceMap[i];
			if(gpuid >= 0)
				printf("GPU %d: core=%dMHz mem=%dMHz powertune=%d fanspeed=%d\n",
						gpuid, gpu_engineclock(i), gpu_memclock(i), gpu_powertune(i), gpu_fanspeed(i));
		}
	
	stats.set_cpd(cpd);
	stats.set_errors(errors);
	stats.set_temp(maxtemp);
	
	mStatCounter++;
	
	return ngpus;
	
}


void XPMClient::Toggle() {
	
	for(unsigned i = 0; i < mWorkers.size(); ++i) {
		if(mWorkers[i].first)
			czmq_signal(mWorkers[i].second);
	}
	
	mPaused = !mPaused;
	
}


void XPMClient::setup_adl(){
	
	init_adl(mNumDevices);
	
	for(unsigned i = 0; i < mNumDevices; ++i){
		
		if(mCoreFreq[i] > 0)
			if(set_engineclock(i, mCoreFreq[i]))
				printf("set_engineclock(%d, %d) failed.\n", i, mCoreFreq[i]);
		if(mMemFreq[i] > 0)
			if(set_memoryclock(i, mMemFreq[i]))
				printf("set_memoryclock(%d, %d) failed.\n", i, mMemFreq[i]);
		if(mPowertune[i] >= -20 && mPowertune[i] <= 20)
			if(set_powertune(i, mPowertune[i]))
				printf("set_powertune(%d, %d) failed.\n", i, mPowertune[i]);
    if (mFanSpeed[i] > 0)
      if(set_fanspeed(i, mFanSpeed[i]))
        printf("set_fanspeed(%d, %d) failed.\n", i, mFanSpeed[i]);
	}
}
