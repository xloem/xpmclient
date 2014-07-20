cmake_minimum_required(VERSION 2.8)

option(STATIC_BUILD "Build with static libraries on Linux")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  ${CMAKE_SOURCE_DIR}/cmake
)

include(FindDependencies)

include_directories(
  ${GMP_INCLUDE_DIRECTORY}
  ${ZMQ_INCLUDE_DIRECTORY}
  ${CZMQ_INCLUDE_DIRECTORY}
  ${PROTOBUF_INCLUDE_DIRECTORY}
  ${OPENCL_INCLUDE_DIRECTORY}
  ${CMAKE_SOURCE_DIR}/config4cpp/include
)

if (WIN32)
  message("Building for Win32")
  add_definitions(-D_WIN32 -D__WINDOWS__ -DLIBCZMQ_EXPORTS)
  if (InstallPrefix)
    include_directories(${InstallPrefix}/include)
  endif()
else()
  message("Building for Linux")
  add_definitions(-DLINUX)
  if (STATIC_BUILD)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libgcc -static-libstdc++")
  endif()
endif()

file(GLOB config4cppFiles "${CMAKE_SOURCE_DIR}/config4cpp/src/*.cpp")
add_library(config4cpp ${config4cppFiles})

add_executable(xpmclient
  adl.cpp
  baseclient.cpp
  prime.cpp
  primes.cpp
  protocol.pb.cpp
  sha256.cpp
  xpmclient.cpp
)

SET_TARGET_PROPERTIES(xpmclient PROPERTIES
  COMPILE_FLAGS "-std=c++0x"
)

target_link_libraries(xpmclient
  ${CZMQ_LIBRARY}
  ${ZMQ_LIBRARY}
  ${SODIUM_LIBRARY}
  ${GMP_LIBRARY}
  ${PROTOBUF_LIBRARY}
  ${OPENCL_LIBRARY}
  config4cpp
)

if (WIN32)
  target_link_libraries(xpmclient ws2_32)
else()
  target_link_libraries(xpmclient dl rt)
  if (STATIC_BUILD)
    set_target_properties(xpmclient PROPERTIES LINK_SEARCH_END_STATIC 1)
  endif()
endif()
