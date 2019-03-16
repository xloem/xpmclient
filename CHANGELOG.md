xpmclient-10.3

1. Speedup:

  - AMD Radeon RX Vega: +15% or more!
  - All NVidia GPUs: 5-10%
  - Other GPUs: 0-3% depends on OS, drivers, etc.
  
2. NVidia version uses CUDA 10.0, RTX2xxx series support added
3. Deploy script added, everybody can build xpmclient binaries using docker (look contrib directory)
4. Protocol changed.

  Previous versions of xpmclient receives signal on every new block incoming, ask pool about work and receive answer from pool (Three ZeroMQ messages). New version protocol include work into "signal" message.
  For connect to old servers, use -c argument:
  
  - xpmclient -c
  
5. Some minor bug fixed
