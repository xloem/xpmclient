if (InstallPrefix)
  # GNU MP
  find_path(GMP_INCLUDE_DIRECTORY gmp.h gmpxx.h PATHS ${InstallPrefix}/include NO_DEFAULT_PATH)
  find_library(GMP_LIBRARY gmp PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
  find_library(GMPXX_LIBRARY gmpxx PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)

  # ZeroMQ
  find_path(ZMQ_INCLUDE_DIRECTORY zmq.h PATHS ${InstallPrefix}/include NO_DEFAULT_PATH)
  if (WIN32)
    find_library(ZMQ_LIBRARY zmq.dll PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
    find_library(SODIUM_LIBRARY sodium PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
  else()
    find_library(ZMQ_LIBRARY zmq PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
    find_library(SODIUM_LIBRARY sodium PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
  endif()

  # CLRadeonExtender
  find_path(CLRX_INCLUDE_DIRECTORY CLRX/amdasm/Assembler.h PATHS ${InstallPrefix}/include NO_DEFAULT_PATH)
  find_path(CLRX_CONFIG_INCLUDE_DIRECTORY CLRX/Config.h PATHS ${InstallPrefix}/include NO_DEFAULT_PATH)
  find_library(CLRX_AMDASM_LIBRARY CLRXAmdAsm PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
  find_library(CLRX_AMDBIN_LIBRARY CLRXAmdBin PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
  find_library(CLRX_UTILS_LIBRARY CLRXUtils PATHS ${InstallPrefix}/lib NO_DEFAULT_PATH)
else()
  # GNU MP
  find_path(GMP_INCLUDE_DIRECTORY gmp.h gmpxx.h)
  find_library(GMP_LIBRARY gmp)
  find_library(GMPXX_LIBRARY gmpxx)

  # ZeroMQ
  find_path(ZMQ_INCLUDE_DIRECTORY zmq.h)
  if (WIN32)
    find_library(ZMQ_LIBRARY zmq.dll)
	find_library(SODIUM_LIBRARY sodium)
  else()
    find_library(ZMQ_LIBRARY zmq)
    find_library(SODIUM_LIBRARY sodium)
  endif()

  # CLRadeonExtender
  find_path(CLRX_INCLUDE_DIRECTORY CLRX/amdasm/Assembler.h)
  find_path(CLRX_CONFIG_INCLUDE_DIRECTORY CLRX/Config.h)
  find_library(CLRX_AMDASM_LIBRARY CLRXAmdAsm)
  find_library(CLRX_AMDBIN_LIBRARY CLRXAmdBin)
  find_library(CLRX_UTILS_LIBRARY CLRXUtils)
endif()

set(CLRX_INCLUDE_DIRECTORIES ${CLRX_INCLUDE_DIRECTORY} ${CLRX_CONFIG_INCLUDE_DIRECTORY})
set(CLRX_LIBRARIES ${CLRX_AMDASM_LIBRARY} ${CLRX_AMDBIN_LIBRARY} ${CLRX_UTILS_LIBRARY})
