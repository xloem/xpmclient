#include "benchmarks.h"

#include "gmpxx.h"

#include <time.h>
#include <chrono>
#include <memory>
#if defined(__GXX_EXPERIMENTAL_CXX0X__) && (__cplusplus < 201103L)
#define steady_clock monotonic_clock
#endif  

#include "xpmclient.h"
 
enum OpenCLKernels {
  CLKernelSquareBenchmark320 = 0,
  CLKernelSquareBenchmark352,  
  CLKernelMultiplyBenchmark320,
  CLKernelMultiplyBenchmark352,
  CLKernelFermatTestBenchmark320,
  CLKernelFermatTestBenchmark352,
  CLKernelHashMod,
  CLKernelsNum
};  
 
static const char *gOpenCLKernelNames[] = {
  "squareBenchmark320",
  "squareBenchmark352",
  "multiplyBenchmark320",
  "multiplyBenchmark352",
  "fermatTestBenchMark320",
  "fermatTestBenchMark352",
  "bhashmod14"
};

const unsigned GroupSize = 256;
const unsigned MulOpsNum = 512; 

uint32_t rand32()
{
  uint32_t result = rand();
  result = (result << 16) | rand();
  return result;
}

uint64_t rand64()
{
  uint64_t result = rand();
  result = (result << 16) | rand();
  result = (result << 16) | rand();
  result = (result << 16);
  return result;
} 
 
void multiplyBenchmark(cl_command_queue queue,
                       cl_kernel *kernels,
                       unsigned groupsNum,                       
                       unsigned mulOperandSize,
                       uint32_t elementsNum,
                       bool isSquaring)
{
  unsigned gmpOpSize = mulOperandSize + (mulOperandSize%2);
  unsigned limbsNum = elementsNum*gmpOpSize;
  clBuffer<uint32_t> m1;
  clBuffer<uint32_t> m2;
  clBuffer<uint32_t> mR;
  clBuffer<uint32_t> cpuR;
  
  m1.init(limbsNum, CL_MEM_READ_WRITE);
  m2.init(limbsNum, CL_MEM_READ_WRITE);
  mR.init(limbsNum*2, CL_MEM_READ_WRITE);
  cpuR.init(limbsNum*2, CL_MEM_READ_WRITE);

  memset(&m1.get(0), 0, limbsNum*sizeof(uint32_t));
  memset(&m2.get(0), 0, limbsNum*sizeof(uint32_t));
  memset(&mR.get(0), 0, 2*limbsNum*sizeof(uint32_t));
  memset(&cpuR.get(0), 0, 2*limbsNum*sizeof(uint32_t));  
  for (unsigned i = 0; i < elementsNum; i++) {
    for (unsigned j = 0; j < mulOperandSize; j++) {
      m1[i*gmpOpSize + j] = rand32();
      m2[i*gmpOpSize + j] = rand32();
    }
  }

  m1.copyToDevice(queue);
  m2.copyToDevice(queue);

  cl_kernel kernel;
  if (isSquaring) {
    if (mulOperandSize == 320/32) {
      kernel = kernels[CLKernelSquareBenchmark320];
    } else if (mulOperandSize == 352/32) {
      kernel = kernels[CLKernelSquareBenchmark352];
    } else {
      fprintf(stderr, "Can't multiply %u-size operands on OpenCL device\n", mulOperandSize*32);
      return;
    }
  } else {
    if (mulOperandSize == 320/32) {
      kernel = kernels[CLKernelMultiplyBenchmark320];
    } else if (mulOperandSize == 352/32) {
      kernel = kernels[CLKernelMultiplyBenchmark352];
    } else {
      fprintf(stderr, "Can't multiply %u-size operands on OpenCL device\n", mulOperandSize*32);
      return;
    }
  }

  
  if (isSquaring) {
    clSetKernelArg(kernel, 0, sizeof(cl_mem), &m1.DeviceData);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), &mR.DeviceData);
    clSetKernelArg(kernel, 2, sizeof(elementsNum), &elementsNum);
  } else {
    clSetKernelArg(kernel, 0, sizeof(cl_mem), &m1.DeviceData);
    clSetKernelArg(kernel, 1, sizeof(cl_mem), &m2.DeviceData);
    clSetKernelArg(kernel, 2, sizeof(cl_mem), &mR.DeviceData);
    clSetKernelArg(kernel, 3, sizeof(elementsNum), &elementsNum);
  }
  
  std::unique_ptr<mpz_class[]> cpuM1(new mpz_class[elementsNum]);
  std::unique_ptr<mpz_class[]> cpuM2(new mpz_class[elementsNum]);
  std::unique_ptr<mpz_class[]> cpuResult(new mpz_class[elementsNum]);
  
  for (unsigned i = 0; i < elementsNum; i++) {
    mpz_import(cpuM1[i].get_mpz_t(), mulOperandSize, -1, 4, 0, 0, &m1[i*gmpOpSize]);
    mpz_import(cpuM2[i].get_mpz_t(), mulOperandSize, -1, 4, 0, 0, &m2[i*gmpOpSize]);
    mpz_import(cpuResult[i].get_mpz_t(), mulOperandSize*2, -1, 4, 0, 0, &mR[i*mulOperandSize*2]);
  }

  clFinish(queue);
  auto gpuBegin = std::chrono::steady_clock::now();  
  
  {
    size_t globalThreads[1] = { groupsNum*GroupSize };
    size_t localThreads[1] = { GroupSize };
    cl_event event;
    cl_int result;
    if ((result = clEnqueueNDRangeKernel(queue,
                                         kernel,
                                         1,
                                         0,
                                         globalThreads,
                                         localThreads,
                                         0, 0, &event)) != CL_SUCCESS) {
      fprintf(stderr, "clEnqueueNDRangeKernel error!\n");
      return;
    }

    if (clWaitForEvents(1, &event) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error!\n");
      return;
    }
    
    clReleaseEvent(event);
  }
  
  auto gpuEnd = std::chrono::steady_clock::now();  
  
  if (isSquaring) {
    for (unsigned i = 0; i < elementsNum; i++) {
      unsigned gmpLimbsNum = cpuM1[i].get_mpz_t()->_mp_size;
      mp_limb_t *Operand1 = cpuM1[i].get_mpz_t()->_mp_d;
      mp_limb_t *target = (mp_limb_t*)&cpuR[i*mulOperandSize*2];
      for (unsigned j = 0; j < MulOpsNum; j++) {
        mpn_sqr(target, Operand1, gmpLimbsNum);
        memcpy(Operand1, target, mulOperandSize*sizeof(uint32_t));
      }
    }
  } else {
    for (unsigned i = 0; i < elementsNum; i++) {
      unsigned gmpLimbsNum = cpuM1[i].get_mpz_t()->_mp_size;
      mp_limb_t *Operand1 = cpuM1[i].get_mpz_t()->_mp_d;
      mp_limb_t *Operand2 = cpuM2[i].get_mpz_t()->_mp_d;
      mp_limb_t *target = (mp_limb_t*)&cpuR[i*mulOperandSize*2];
      for (unsigned j = 0; j < MulOpsNum; j++) {
        mpn_mul_n(target, Operand1, Operand2, gmpLimbsNum);
        memcpy(Operand1, target, mulOperandSize*sizeof(uint32_t));
      }
    }
  }

  mR.copyToHost(queue);
  clFinish(queue);

  for (unsigned i = 0; i < elementsNum; i++) {
    if (memcmp(&mR[i*mulOperandSize*2], &cpuR[i*mulOperandSize*2], 4*mulOperandSize*2) != 0) {
      fprintf(stderr, "element index: %u\n", i);
      fprintf(stderr, "gmp: ");
      for (unsigned j = 0; j < mulOperandSize*2; j++)
        fprintf(stderr, "%08X ", cpuR[i*mulOperandSize*2 + j]);
      fprintf(stderr, "\ngpu: ");
      for (unsigned j = 0; j < mulOperandSize*2; j++)
        fprintf(stderr, "%08X ", mR[i*mulOperandSize*2 + j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "results differ!\n");
      break;
    }
  }

  double gpuTime = std::chrono::duration_cast<std::chrono::microseconds>(gpuEnd-gpuBegin).count() / 1000.0;  
  double opsNum = ((elementsNum*MulOpsNum) / 1000000.0) / gpuTime * 1000.0;
  
  printf("%s %u bits: %.3lfms (%.3lfM ops/sec)\n", (isSquaring ? "square" : "multiply"), mulOperandSize*32, gpuTime, opsNum);
}


void fermatTestBenchmark(cl_command_queue queue,
                         cl_kernel *kernels,
                         unsigned groupsNum, 
                         unsigned operandSize,
                         unsigned elementsNum)
{ 
  unsigned numberLimbsNum = elementsNum*operandSize;
  
  clBuffer<uint32_t> numbers;
  clBuffer<uint32_t> gpuResults;
  clBuffer<uint32_t> cpuResults;
  
  numbers.init(numberLimbsNum, CL_MEM_READ_WRITE);
  gpuResults.init(numberLimbsNum, CL_MEM_READ_WRITE);
  cpuResults.init(numberLimbsNum, CL_MEM_READ_WRITE);
  
  for (unsigned i = 0; i < elementsNum; i++) {
    for (unsigned j = 0; j < operandSize; j++)
      numbers[i*operandSize + j] = (j == operandSize-1) ? (1 << (i % 32)) : rand32();
    numbers[i*operandSize] |= 0x1; 
  }

  numbers.copyToDevice(queue);
  gpuResults.copyToDevice(queue);

  cl_kernel kernel;
  if (operandSize == 320/32) {
    kernel = kernels[CLKernelFermatTestBenchmark320];
  } else if (operandSize == 352/32) {
    kernel = kernels[CLKernelFermatTestBenchmark352];
  } else {
    fprintf(stderr, "Can't do Fermat test on %ubit operand\n", operandSize*32);
    return;
  }
  
  clSetKernelArg(kernel, 0, sizeof(cl_mem), &numbers.DeviceData);
  clSetKernelArg(kernel, 1, sizeof(cl_mem), &gpuResults.DeviceData);
  clSetKernelArg(kernel, 2, sizeof(elementsNum), &elementsNum);
  
  std::unique_ptr<mpz_t[]> cpuNumbersBuffer(new mpz_t[elementsNum]);
  std::unique_ptr<mpz_t[]> cpuResultsBuffer(new mpz_t[elementsNum]);
  mpz_class mpzTwo = 2;
  mpz_class mpzE;
  mpz_import(mpzE.get_mpz_t(), operandSize, -1, 4, 0, 0, &numbers[0]);
  for (unsigned i = 0; i < elementsNum; i++) {
    mpz_init(cpuNumbersBuffer[i]);
    mpz_init(cpuResultsBuffer[i]);
    mpz_import(cpuNumbersBuffer[i], operandSize, -1, 4, 0, 0, &numbers[i*operandSize]);
    mpz_import(cpuResultsBuffer[i], operandSize, -1, 4, 0, 0, &cpuResults[i*operandSize]);
  }
  
  clFinish(queue);
  auto gpuBegin = std::chrono::steady_clock::now();  

  {
    size_t globalThreads[1] = { groupsNum*GroupSize };
    size_t localThreads[1] = { GroupSize };
    cl_event event;
    cl_int result;
    if ((result = clEnqueueNDRangeKernel(queue,
                                         kernel,
                                         1,
                                         0,
                                         globalThreads,
                                         localThreads,
                                         0, 0, &event)) != CL_SUCCESS) {
      fprintf(stderr, "clEnqueueNDRangeKernel error!\n");
      return;
    }
      
    cl_int error;
    if ((error = clWaitForEvents(1, &event)) != CL_SUCCESS) {
      fprintf(stderr, "clWaitForEvents error %i!\n", error);
      return;
    }
      
    clReleaseEvent(event);
  }
  
  auto gpuEnd = std::chrono::steady_clock::now();  
  
  
  for (unsigned i = 0; i < elementsNum; i++) {
    mpz_sub_ui(mpzE.get_mpz_t(), cpuNumbersBuffer[i], 1);
    mpz_powm(cpuResultsBuffer[i], mpzTwo.get_mpz_t(), mpzE.get_mpz_t(), cpuNumbersBuffer[i]);
  }

  gpuResults.copyToHost(queue);
  clFinish(queue);
  
  memset(&cpuResults[0], 0, 4*operandSize*elementsNum);
  for (unsigned i = 0; i < elementsNum; i++) {
    size_t exportedLimbs;
    mpz_export(&cpuResults[i*operandSize], &exportedLimbs, -1, 4, 0, 0, cpuResultsBuffer[i]);
    if (memcmp(&gpuResults[i*operandSize], &cpuResults[i*operandSize], 4*operandSize) != 0) {
      fprintf(stderr, "element index: %u\n", i);
      fprintf(stderr, "gmp: ");
      for (unsigned j = 0; j < operandSize; j++)
        fprintf(stderr, "%08X ", cpuResults[i*operandSize + j]);
      fprintf(stderr, "\ngpu: ");
      for (unsigned j = 0; j < operandSize; j++)
        fprintf(stderr, "%08X ", gpuResults[i*operandSize + j]);
      fprintf(stderr, "\n");
      fprintf(stderr, "results differ!\n");
      break;
    }
  }
  
  double gpuTime = std::chrono::duration_cast<std::chrono::microseconds>(gpuEnd-gpuBegin).count() / 1000.0;  
  double opsNum = ((elementsNum) / 1000000.0) / gpuTime * 1000.0;
  
  printf("%s %u bits: %.3lfms (%.3lfM ops/sec)\n", "Fermat tests", operandSize*32, gpuTime, opsNum);
}


void hashmodBenchmark(cl_command_queue queue, cl_kernel *kernels, unsigned groupsNum)
{
  const unsigned iterationsNum = 64;
  cl_kernel mHashMod = kernels[CLKernelHashMod];
  
  PrimeMiner::search_t hashmod;
  PrimeMiner::block_t blockheader;
  
  hashmod.midstate.init(8*sizeof(cl_uint), CL_MEM_READ_ONLY);
  hashmod.found.init(2048, CL_MEM_READ_WRITE);
  hashmod.primorialBitField.init(2048, CL_MEM_READ_WRITE);
  hashmod.count.init(1, CL_MEM_READ_WRITE);

  clSetKernelArg(mHashMod, 0, sizeof(cl_mem), &hashmod.found.DeviceData);
  clSetKernelArg(mHashMod, 1, sizeof(cl_mem), &hashmod.count.DeviceData);
  clSetKernelArg(mHashMod, 2, sizeof(cl_mem), &hashmod.primorialBitField.DeviceData);
  clSetKernelArg(mHashMod, 3, sizeof(cl_mem), &hashmod.midstate.DeviceData);

  uint64_t totalTime = 0;
  unsigned totalHashes = 0;
  int numhash = 64 * 131072;

  unsigned hashm[32];
  memset(hashm, 0, sizeof(hashm));
  
  for (unsigned i = 0; i < iterationsNum; i++) {
    {
      uint8_t *pHeader = (uint8_t*)&blockheader;
      for (unsigned i = 0; i < sizeof(blockheader); i++)
        pHeader[i] = rand();
      blockheader.version = PrimeMiner::block_t::CURRENT_VERSION;
      blockheader.nonce = 1;    
    }    
    
    cl_uint msg_merkle = blockheader.hashMerkleRoot.Get64(3) >> 32;
    cl_uint msg_time = blockheader.time;
    cl_uint msg_bits = blockheader.bits;
 
    clSetKernelArg(mHashMod, 4, sizeof(cl_uint), &msg_merkle);
    clSetKernelArg(mHashMod, 5, sizeof(cl_uint), &msg_time);
    clSetKernelArg(mHashMod, 6, sizeof(cl_uint), &msg_bits);    
    hashmod.count.copyToDevice(queue, false);
    
    size_t globalSize[] = { numhash, 1, 1 };
 
    hashmod.count[0] = 0;
    hashmod.count.copyToDevice(queue);
    clFinish(queue);
    auto gpuBegin = std::chrono::steady_clock::now();  
    
    {
      cl_event event;
      cl_int result;
      if ((result = clEnqueueNDRangeKernel(queue,
                                           mHashMod,
                                           1,
                                           0,
                                           globalSize,
                                           0,
                                           0,
                                           0, &event)) != CL_SUCCESS) {
        fprintf(stderr, "clEnqueueNDRangeKernel error!\n");
        return;
      }
        
      cl_int error;
      if ((error = clWaitForEvents(1, &event)) != CL_SUCCESS) {
        fprintf(stderr, "clWaitForEvents error %i!\n", error);
        return;
      }
        
      clReleaseEvent(event);
    } 
    
    auto gpuEnd = std::chrono::steady_clock::now();  
    
    hashmod.found.copyToHost(queue, false);
    hashmod.primorialBitField.copyToHost(queue, false);
    hashmod.count.copyToHost(queue, false);
    clFinish(queue);
    
    totalTime += std::chrono::duration_cast<std::chrono::microseconds>(gpuEnd-gpuBegin).count();
    totalHashes += hashmod.count[0];
    
    for (unsigned i = 0; i < hashmod.count[0]; i++) {
      uint32_t multiplierBitField = hashmod.primorialBitField[i];
      
      unsigned multiplierCount = 0;
      for (unsigned j = 0; j < 32; j++)
        multiplierCount += ((multiplierBitField & (1 << j)) != 0);
      
      hashm[multiplierCount]++;
    }
    
//     printf("found %u hashes\n", hashmod.count[0]);
  }
  
  double averageHashes = (double)totalHashes / iterationsNum;
  printf("\n *** hashmod benchmark ***\n");
  printf(" MHash per second: %.3lf\n", iterationsNum*numhash / (double)totalTime);
  printf(" Hash per iteration: %.3lf (%.6lf %%)\n", averageHashes, averageHashes*100/numhash);
  printf(" Hashes by multipliers count:\n");
  for (unsigned i = 0; i < 32; i++) {
    if (hashm[i])
      printf("   * [%u] %.3lf (%.3lf%%)\n", i, (double)hashm[i] / iterationsNum, (double)hashm[i] / (double)totalHashes * 100.0);
  }
  
}

void sieveTestBenchmark(cl_command_queue queue, cl_kernel *kernels, unsigned groupsNum)
{
  clBuffer<cl_uint> midstate;
  clBuffer<cl_uint> found;
  clBuffer<cl_uint> primorialBitField;
  clBuffer<cl_uint> count;
  
  midstate.init(8*sizeof(cl_uint), CL_MEM_READ_ONLY);
  found.init(2048, CL_MEM_READ_WRITE);
  primorialBitField.init(2048, CL_MEM_READ_WRITE);
  count.init(1, CL_MEM_READ_WRITE);
//   hashBuf.init(PW*mConfig.N, CL_MEM_ALLOC_HOST_PTR | CL_MEM_READ_ONLY);  
}

void runBenchmarks(cl_context context, cl_program program, cl_device_id deviceId)
{
  char deviceName[128] = {0};
  cl_uint computeUnits;

  clGetDeviceInfo(deviceId, CL_DEVICE_NAME, sizeof(deviceName), deviceName, 0);
  clGetDeviceInfo(deviceId, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(computeUnits), &computeUnits, 0);
  printf("%s; %u compute units\n", deviceName, computeUnits);  
  
  std::unique_ptr<cl_kernel[]> kernels(new cl_kernel[CLKernelsNum]);
  for (unsigned i = 0; i < CLKernelsNum; i++) {
    cl_int clResult;
    kernels[i] = clCreateKernel(program, gOpenCLKernelNames[i], &clResult);
    if (clResult != CL_SUCCESS) {
      fprintf(stderr, " * Error: can't found kernel %s\n", gOpenCLKernelNames[i]);
      return;
    }
  }
  
  cl_int error;  
  cl_command_queue queue = clCreateCommandQueue(context, deviceId, 0, &error);
  if (!queue || error != CL_SUCCESS) {
    fprintf(stderr, " * Error: can't create command queue\n");
    return;
  }
    
  multiplyBenchmark(queue, kernels.get(), computeUnits*4, 320/32, 262144, true);  
   
  multiplyBenchmark(queue, kernels.get(), computeUnits*4, 320/32, 262144, true);
  multiplyBenchmark(queue, kernels.get(), computeUnits*4, 320/32, 262144, false);
  multiplyBenchmark(queue, kernels.get(), computeUnits*4, 352/32, 262144, true);    
  multiplyBenchmark(queue, kernels.get(), computeUnits*4, 352/32, 262144, false);  

  fermatTestBenchmark(queue, kernels.get(), computeUnits*4, 320/32, 131072);
  fermatTestBenchmark(queue, kernels.get(), computeUnits*4, 352/32, 131072);   
  
  hashmodBenchmark(queue, kernels.get(), 0);
}
