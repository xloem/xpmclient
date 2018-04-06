#include "cudautil.h"
#include <fstream>
#include <iostream>
#include <memory>

bool cudaCompileKernel(const char *kernelName,
                       const std::vector<const char*> &sources,
                       const char **arguments,
                       int argumentsNum,
                       CUmodule *module,
                       bool needRebuild) 
{
  std::ifstream testfile(kernelName);
  if(needRebuild || !testfile) {
    printf("<info> compiling ...\n");
    
    std::string sourceFile;
    for (auto &i: sources) {
      std::ifstream stream(i);
      std::string str((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());
      sourceFile.append(str);
    }
    
    printf("<info> source: %u bytes\n", (unsigned)sourceFile.size());
    if(sourceFile.size() < 1){
      fprintf(stderr, "<error> source files not found or empty\n");
      return false;
    }
    
    nvrtcProgram prog;
    NVRTC_SAFE_CALL(
      nvrtcCreateProgram(&prog,
                         sourceFile.c_str(),
                         "xpm.cu",
                         0,
                         NULL,
                         NULL));

    nvrtcResult compileResult = nvrtcCompileProgram(prog, argumentsNum, arguments);
    
    // Obtain compilation log from the program.
    size_t logSize;
    NVRTC_SAFE_CALL(nvrtcGetProgramLogSize(prog, &logSize));
    char *log = new char[logSize];
    NVRTC_SAFE_CALL(nvrtcGetProgramLog(prog, log));
    std::cout << log << '\n';
    delete[] log;
    if (compileResult != NVRTC_SUCCESS) {
      return false;
    }
    
    // Obtain PTX from the program.
    size_t ptxSize;
    NVRTC_SAFE_CALL(nvrtcGetPTXSize(prog, &ptxSize));
    char *ptx = new char[ptxSize];
    NVRTC_SAFE_CALL(nvrtcGetPTX(prog, ptx));
    
    // Destroy the program.
    NVRTC_SAFE_CALL(nvrtcDestroyProgram(&prog));
    
    {
      std::ofstream bin(kernelName, std::ofstream::binary | std::ofstream::trunc);
      bin.write(ptx, ptxSize);
      bin.close();      
    }
    
    delete[] ptx;
  }
  
  std::ifstream bfile(kernelName, std::ifstream::binary);
  if(!bfile) {
    printf("<error> %s not found\n", kernelName);
    return false;
  }  
  
  bfile.seekg(0, bfile.end);
  size_t binsize = bfile.tellg();
  bfile.seekg(0, bfile.beg);
  if(!binsize){
    printf("<error> %s empty\n", kernelName);
    return false;
  }
  
  std::unique_ptr<char[]> ptx(new char[binsize+1]);
  bfile.read(ptx.get(), binsize);
  bfile.close();
  
  CUDA_SAFE_CALL(cuModuleLoadDataEx(module, ptx.get(), 0, 0, 0));
  return true;
}
