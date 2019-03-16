#!/bin/bash

export DIST_LANG="cpp"

# Linux static build

# gmp
cd /home/user/build/deps-linux
tar --lzip -xvf ../gmp-6.1.2.tar.lz
cd gmp-6.1.2
./configure --build corei7 --prefix=/home/user/install/x86_64-Linux --enable-cxx --enable-static --disable-shared 
make -j4
make install

# sodium
cd /home/user/build/deps-linux
tar -xzf ../libsodium-1.0.17.tar.gz
cd libsodium-1.0.17
./configure --prefix=/home/user/install/x86_64-Linux --enable-static --disable-shared 
make -j4
make install

# zmq
cd /home/user/build/deps-linux
tar -xzf ../zeromq-4.3.1.tar.gz
cd zeromq-4.3.1
# disable glibc 2.25 'getrandom' usage for compatibility
sed -i 's/libzmq_cv_getrandom=\"yes\"/libzmq_cv_getrandom=\"no\"/' configure
./configure --prefix=/home/user/install/x86_64-Linux --enable-static --disable-shared 
make -j4
make install

# protobuf
cd /home/user/build/deps-linux
tar -xzf ../protobuf-cpp-3.6.1.tar.gz
cd protobuf-3.6.1
./configure --prefix=/home/user/install/x86_64-Linux --enable-static --disable-shared
make -j4
make install

# xpmclient
mkdir /home/user/build/xpmclient/x86_64-Linux
cd /home/user/build/xpmclient/x86_64-Linux
cmake ../src -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home/user/install/x86_64-Linux -DSTATIC_BUILD=ON -DOPENCL_INCLUDE_DIRECTORY=/usr/local/cuda-10.0/include -DOPENCL_LIBRARY=/usr/local/cuda-10.0/lib64/libOpenCL.so -DCUDA_driver_LIBRARY=/usr/local/cuda-10.0/compat/libcuda.so
make -j4 xpmclient xpmclientnv
strip xpmclient
strip xpmclientnv

# make AMD(OpenCL) distr
mkdir xpmclient-opencl-10.3-linux
cd xpmclient-opencl-10.3-linux
cp ../xpmclient .
cp ../../src/xpm/opencl/config.txt .
mkdir -p xpm/opencl
cp ../../src/xpm/opencl/*.cl xpm/opencl
cd ..
tar -czf xpmclient-opencl-10.3-linux.tar.gz xpmclient-opencl-10.3-linux

# make NVidia distr
mkdir xpmclient-cuda-10.3-linux
cd xpmclient-cuda-10.3-linux
cp ../xpmclientnv ./miner
echo "#/bin/bash" > xpmclientnv
echo "DIR=\$(dirname \"\$0\")" >> xpmclientnv
echo "LD_LIBRARY_PATH=\$DIR/. ./miner \$@" >> xpmclientnv
chmod +x xpmclientnv
cp ../../src/xpm/cuda/config.txt .
mkdir -p xpm/cuda
cp ../../src/xpm/cuda/*.cu xpm/cuda
cp /usr/local/cuda-10.0/lib64/libnvrtc.so.10.0 .
cp /usr/local/cuda-10.0/lib64/libnvrtc-builtins.so .
cd ..
tar -czf xpmclient-cuda-10.3-linux.tar.gz xpmclient-cuda-10.3-linux

# Windows (MingGW) static build

# gmp
cd /home/user/build/deps-win32
tar --lzip -xvf ../gmp-6.1.2.tar.lz
cd gmp-6.1.2
./configure --host=x86_64-w64-mingw32 --build corei7 --prefix=/home/user/install/x86_64-w64-mingw32 --enable-cxx --enable-static --disable-shared
make -j4
make install

# sodium
cd /home/user/build/deps-win32
tar -xzf ../libsodium-1.0.17.tar.gz
cd libsodium-1.0.17
./configure --host=x86_64-w64-mingw32 --prefix=/home/user/install/x86_64-w64-mingw32 --enable-static --disable-shared 
make -j4
make install

# zmq (dll)
cd /home/user/build/deps-win32
tar -xzf ../zeromq-4.3.1.tar.gz
cd zeromq-4.3.1
./configure --host=x86_64-w64-mingw32 --prefix=/home/user/install/x86_64-w64-mingw32
make -j4
make install

# protobuf
cd /home/user/build/deps-win32
tar -xzf ../protobuf-cpp-3.6.1.tar.gz
cd protobuf-3.6.1
./configure --host=x86_64-w64-mingw32 --prefix=/home/user/install/x86_64-w64-mingw32 --enable-static --disable-shared
make -j4
make install

# xpmclient
mkdir /home/user/build/xpmclient/x86_64-w64-mingw32
cd /home/user/build/xpmclient/x86_64-w64-mingw32
cmake ../src -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=../src/cmake/Toolchain-x86_64-w64-mingw32.cmake -DCMAKE_INSTALL_PREFIX=/home/user/install/x86_64-w64-mingw32 -DSTATIC_BUILD=ON -DProtobuf_PROTOC_EXECUTABLE=/home/user/install/x86_64-Linux/bin/protoc -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-10.0-win32 -DOPENCL_INCLUDE_DIRECTORY=/usr/local/cuda-10.0-win32/include -DOPENCL_LIBRARY=/usr/local/cuda-10.0-win32/lib/x64/OpenCL.lib
make -j4 xpmclient xpmclientnv
x86_64-w64-mingw32-strip xpmclient.exe
x86_64-w64-mingw32-strip xpmclientnv.exe

# make AMD(OpenCL) distr
mkdir xpmclient-opencl-10.3-win64
cd xpmclient-opencl-10.3-win64
cp ../xpmclient.exe .
cp ../../src/xpm/opencl/config.txt .
mkdir -p xpm/opencl
cp ../../src/xpm/opencl/*.cl xpm/opencl
cp /usr/lib/gcc/x86_64-w64-mingw32/7.3-posix/libgcc_s_seh-1.dll .
cp /usr/lib/gcc/x86_64-w64-mingw32/7.3-posix/libstdc++-6.dll .
cp /usr/x86_64-w64-mingw32/lib/libwinpthread-1.dll .
cp /home/user/install/x86_64-w64-mingw32/bin/libzmq.dll .
cd ..
zip -9 -r xpmclient-opencl-10.3-win64.zip xpmclient-opencl-10.3-win64

# make Nvidia distr
mkdir xpmclient-cuda-10.3-win64
cd xpmclient-cuda-10.3-win64
cp ../xpmclientnv.exe .
cp ../../src/xpm/cuda/config.txt .
mkdir -p xpm/cuda
cp ../../src/xpm/cuda/*.cu xpm/cuda
cp /usr/lib/gcc/x86_64-w64-mingw32/7.3-posix/libgcc_s_seh-1.dll .
cp /usr/lib/gcc/x86_64-w64-mingw32/7.3-posix/libstdc++-6.dll .
cp /usr/x86_64-w64-mingw32/lib/libwinpthread-1.dll .
cp /home/user/install/x86_64-w64-mingw32/bin/libzmq.dll .
cp /usr/local/cuda-10.0-win32/bin/nvrtc64_100_0.dll .
cp /usr/local/cuda-10.0-win32/bin/nvrtc-builtins64_100.dll .
cd ..
zip -9 -r xpmclient-cuda-10.3-win64.zip xpmclient-cuda-10.3-win64

# Calculate SHA256 checksum
cd /home/user/build/xpmclient
sha256sum /home/user/build/xpmclient/x86_64-Linux/xpmclient-opencl-10.3-linux.tar.gz > xpmclient-10.3-sha256.txt
sha256sum /home/user/build/xpmclient/x86_64-Linux/xpmclient-cuda-10.3-linux.tar.gz >> xpmclient-10.3-sha256.txt
sha256sum /home/user/build/xpmclient/x86_64-w64-mingw32/xpmclient-opencl-10.3-win64.zip >> xpmclient-10.3-sha256.txt
sha256sum /home/user/build/xpmclient/x86_64-w64-mingw32/xpmclient-cuda-10.3-win64.zip >> xpmclient-10.3-sha256.txt
