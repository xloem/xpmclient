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

bool clCompileKernel(cl_context gContext,
                     cl_device_id gpu,
                     const char *binaryName,
                     const char **sources,
                     unsigned sourcesNum,
                     const char *arguments,
                     cl_int *binstatus,
                     cl_program *gProgram,
                     bool needRebuild)
{
  std::ifstream testfile(binaryName);
  if(needRebuild || !testfile) {
    LOG_F(INFO, "compiling ...");

    std::unique_ptr<std::string[]> sourceData(new std::string[sourcesNum]);
    std::unique_ptr<const char*[]> sourcesPtr(new const char*[sourcesNum+1]);
    for (unsigned i = 0; i < sourcesNum; i++) {
      std::ifstream stream(sources[i]);
      std::string data((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
      sourceData[i] = data;
      sourcesPtr[i] = sourceData[i].c_str();
    }

    sourcesPtr[sourcesNum] = nullptr;

    cl_int error;
    *gProgram = clCreateProgramWithSource(gContext, sourcesNum, sourcesPtr.get(), nullptr, &error);
    OCLR(error, false);
    
    if (clBuildProgram(*gProgram, 1, &gpu, arguments, 0, 0) != CL_SUCCESS) {
      size_t logSize;
      clGetProgramBuildInfo(*gProgram, gpu, CL_PROGRAM_BUILD_LOG, 0, 0, &logSize);
      
      std::unique_ptr<char[]> log(new char[logSize]);
      clGetProgramBuildInfo(*gProgram, gpu, CL_PROGRAM_BUILD_LOG, logSize, log.get(), 0);
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
  
  *gProgram = clCreateProgramWithBinary(gContext, 1, &gpu, &binsize, &binaryPtr, binstatus, &error);
  OCLR(error, false);
  OCLR(clBuildProgram(*gProgram, 1, &gpu, 0, 0, 0), false);  
  return true;
}
