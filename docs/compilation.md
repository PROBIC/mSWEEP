# Compilation tips for mSWEEP
## Improving performance
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
mSWEEP can be compiled with link time optimization. This can slightly improve the performance on large problems and can be enabled with
```
> cmake -DCMAKE_BUILD_WITH_FLTO=1 ..
```

### Intel C compiler
If the [Intel C++ compiler](https://software.intel.com/en-us/c-compilers) is available in your environment, you might want to use that to compile mSWEEP — especially if running on Intel hardware. The compiler can be specified by running
```
> cmake -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc ..
```

## MPI support
mSWEEP can be compiled with MPI support, distributing the mixture
component estimation part of the program to several processes. To
compile with MPI support, set your environment appropriately and build
mSWEEP with the following commands (example case for OpenMPI):
```
> mkdir build
> cd build
> module load mpi/openmpi
> cmake -DCMAKE_ENABLE_MPI_SUPPORT=1 -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx ..
> make
```

The project should configure itself appropriately. To distribute the
computation to 4 processes after compiling, call mSWEEP with:
```
mpirun -np 4 mSWEEP --themisto-1 forward_aln.gz --themisto-2 reverse_aln.gz -i cluster_indicators.txt
```

In some cases it might be useful to use hybrid parallellization with
multiple threads per process. This can be accomplished through use of
the `-t` flag:
```
mpirun -np 2 mSWEEP --themisto-1 forward_aln.gz --themisto-2 reverse_aln.gz -i cluster_indicators.txt -t 2
```

which will distribute computation to two processes with two
threads. The optimal configuration will depend on the size and
structure of your data.

Note that when mSWEEP is called through MPI, the root process will
handle all read and write operations and only the estimation part is
distributed.

WARNING: Using mSWEEP through MPI will result in increased memory usage.
