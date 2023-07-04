# Compiling mSWEEP
## Compiler flags
### CPU instructions
If you intend to run mSWEEP on the machine used in compiling the
source code, you might want to add the '-march=native -mtune=native'
flags if compiling with GCC by running
```
> cmake -DCMAKE_CXX_FLAGS="-march=native -mtune=native" -DCMAKE_C_FLAGS="-march=native -mtune=native" ..
```
Using these options significantly reduces the runtime of mSWEEP in
some environments (e.g. most HPC setups).

### Link time optimization
mSWEEP can be compiled with link time optimization. This can improve performance on large problems and can be enabled with
```
> cmake -DCMAKE_BUILD_WITH_FLTO=1 ..
```

## Using the Intel C compiler
The [Intel C++ compiler](https://software.intel.com/en-us/c-compilers) can be enabled by compiling mSWEEP with
```
> cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ..
```

## MPI support (experimental)
mSWEEP can be compiled with MPI support, distributing the mixture
component estimation part of the program to several processes. To
compile with MPI support, set your environment appropriately and build
mSWEEP with the following commands:
```
> mkdir build
> cd build
> module load mpi/openmpi
> cmake -DCMAKE_ENABLE_MPI_SUPPORT=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
> make
```

The `CMAKE_ENABLE_MPI_SUPPORT` flag configures the project accordingly.

### Running with MPI
Using mSWEEP through MPI will result in increased memory usage.

Distribute computation to 4 processes by calling mSWEEP with:
```
mpirun -np 4 mSWEEP --themisto-1 fwd.txt --themisto-2 rev.txt -i cluster_indicators.txt
```

Enable hybrid parallellization with multiple threads per process through use of the `-t` flag:
```
mpirun -np 2 mSWEEP --themisto-1 forward_aln.gz --themisto-2 reverse_aln.gz -i cluster_indicators.txt -t 2
```
Hybrid parallelization might require binding the processes (refer to your HPC documentation on how to do this).
The optimal configuration between ranks and threads will depend on the size and
structure of your data.

When mSWEEP is called through MPI, the root process will
handle all read and write operations and only the estimation part is
distributed.
