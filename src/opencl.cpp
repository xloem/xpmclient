#include "opencl.h"
#include <fstream>
#include <vector>
#include <memory>

extern cl_platform_id gPlatform;

bool clInitialize(const char *requiredPlatform, std::vector<cl_device_id> &gpus)
{
  cl_platform_id platforms[64];
  cl_uint numPlatforms;
  OCLR(clGetPlatformIDs(sizeof(platforms)/sizeof(cl_platform_id), platforms, &numPlatforms), false);
  if (!numPlatforms) {
    LOG_F(ERROR, "no OpenCL platforms found");
    return false;
  }
  
  int platformIdx = -1;
  if (requiredPlatform) {
    for (decltype(numPlatforms) i = 0; i < numPlatforms; i++) {
      char name[1024] = {0};
      OCLR(clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, sizeof(name), name, 0), false);
      LOG_F(INFO, "found platform[%i] name = '%s'", (int)i, name);
      if (strcmp(name, requiredPlatform) == 0) {
        platformIdx = i;
        break;
      }
    }
  } else {
    platformIdx = 0;
  }
  
  
  if (platformIdx == -1) {
    LOG_F(ERROR, "platform %s not exists", requiredPlatform);
    return false;
  }
  
  gPlatform = platforms[platformIdx];
  
  cl_uint numDevices = 0;
  cl_device_id devices[64];
  clGetDeviceIDs(gPlatform, CL_DEVICE_TYPE_GPU, sizeof(devices)/sizeof(cl_device_id), devices, &numDevices);
  if (numDevices) {
    LOG_F(INFO, "found %d devices", numDevices);
  } else {
    LOG_F(ERROR, "no OpenCL GPU devices found.");
    return false;
  }

  for (decltype(numDevices) i = 0; i < numDevices; i++) {
    gpus.push_back(devices[i]);
  }
  
  return true;
}

bool fetchOfflineDevicesAMD(cl_platform_id pid,
                            const char *required,
                            cl_context *offlineCtx,
                            cl_device_id *offlineDevice)
{
#define CL_CONTEXT_OFFLINE_DEVICES_AMD 0x403F

  cl_device_id *offlineDevices = nullptr;
  cl_context ctx = nullptr;
  cl_uint devicesNum = 0;

  cl_context_properties ctxpft[] = {
    CL_CONTEXT_PLATFORM, (cl_context_properties)pid,
    CL_CONTEXT_OFFLINE_DEVICES_AMD, (cl_context_properties)CL_TRUE,
    0
  };

  bool result = false;
  cl_int error;
  if (ctx = clCreateContextFromType(ctxpft, CL_DEVICE_TYPE_ALL, NULL, NULL, &error)) {
    OCL(clGetContextInfo(ctx, CL_CONTEXT_NUM_DEVICES, sizeof(devicesNum), &devicesNum, NULL));

    std::unique_ptr<cl_device_id[]> offlineDevices(new cl_device_id[devicesNum]);
    OCL(clGetContextInfo(ctx, CL_CONTEXT_DEVICES, devicesNum*sizeof(cl_device_id), offlineDevices.get(), NULL));
    LOG_F(INFO, "List available virtual devices:");
    for (cl_uint i = 0; i < devicesNum; i++) {
      char deviceName[128] = {0};
      clGetDeviceInfo(offlineDevices[i], CL_DEVICE_NAME, sizeof(deviceName), deviceName, 0);
      LOG_F(INFO, " * %s", deviceName);
      if (strcmp(deviceName, required) == 0) {
        *offlineCtx = ctx;
        *offlineDevice = offlineDevices[i];
        result = true;
      }
    }
  }

  if (!result)
    clReleaseContext(ctx);
  return result;
}

bool clCompileKernel(cl_context buildContext,
                     cl_device_id buildDevice,
                     cl_context runningContext,
                     cl_device_id runningDevice,
                     const char *binaryName,
                     const std::vector<const char*> &sources,
                     const char *arguments,
                     cl_int *binstatus,
                     cl_program *gProgram,
                     bool needRebuild)
{
  std::ifstream testfile(binaryName);
  if(needRebuild || !testfile) {
    LOG_F(INFO, "compiling ...");
    
    std::string sourceFile;
    for (auto &i: sources) {
      std::ifstream stream(i);
      std::string str((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
      sourceFile.append(str);
    }
    
    LOG_F(INFO, "source: %u bytes", (unsigned)sourceFile.size());
    if(sourceFile.size() < 1){
      LOG_F(ERROR, "source files not found or empty");
      return false;
    }
    
    cl_int error;
    const char *sources[] = { sourceFile.c_str(), 0 };
    *gProgram = clCreateProgramWithSource(buildContext, 1, sources, 0, &error);
    OCLR(error, false);
    
    if (clBuildProgram(*gProgram, 1, &buildDevice, arguments, 0, 0) != CL_SUCCESS) {
      size_t logSize;
      clGetProgramBuildInfo(*gProgram, buildDevice, CL_PROGRAM_BUILD_LOG, 0, 0, &logSize);
      
      std::unique_ptr<char[]> log(new char[logSize]);
      clGetProgramBuildInfo(*gProgram, buildDevice, CL_PROGRAM_BUILD_LOG, logSize, log.get(), 0);
      LOG_F(INFO, "Error: %s\n", log.get());

      return false;
    }
    
    size_t binsize;
    OCLR(clGetProgramInfo(*gProgram, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &binsize, 0), false);
      if(!binsize) {
        LOG_F(ERROR, "no binary available!");
        return false;
      }
    
    LOG_F(INFO, "binsize = %u bytes", (unsigned)binsize);
    std::unique_ptr<unsigned char[]> binary(new unsigned char[binsize+1]);
    OCLR(clGetProgramInfo(*gProgram, CL_PROGRAM_BINARIES, sizeof(void*), &binary, 0), false);
    
    {
      std::ofstream bin(binaryName, std::ofstream::binary | std::ofstream::trunc);
      bin.write((const char*)binary.get(), binsize);
      bin.close();      
    }
   
    OCLR(clReleaseProgram(*gProgram), false);
  }
  
  std::ifstream bfile(binaryName, std::ifstream::binary);
  if(!bfile) {
    LOG_F(ERROR, "%s not found", binaryName);
    return false;
  }  
  
  bfile.seekg(0, bfile.end);
  size_t binsize = bfile.tellg();
  bfile.seekg(0, bfile.beg);
  if(!binsize){
    LOG_F(ERROR, "<error> %s empty", binaryName);
    return false;
  }
  
  std::vector<char> binary(binsize+1);
  bfile.read(&binary[0], binsize);
  bfile.close();
  
  cl_int error;
  const unsigned char *binaryPtr = (const unsigned char*)&binary[0];
  
  *gProgram = clCreateProgramWithBinary(runningContext, 1, &runningDevice, &binsize, &binaryPtr, binstatus, &error);
  OCLR(error, false);
  OCLR(clBuildProgram(*gProgram, 1, &runningDevice, 0, 0, 0), false);
  return true;
}
