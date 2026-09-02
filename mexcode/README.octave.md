# BUILDING THE MEX FILES OF THE TOOLBOX WITH LINKS TO BLAS AND LAPACK

Almost all mex functions to be compiled rely in some way on matrix manipulations and linear algebra, so that a proper implementation for them is crucial. This is solved by systematically calling BLAS/LAPACK routines. As opposed to Matlab, GNU Linux's Octave does not embed itself a particular implementation of BLAS, but instead it depends on that provided by the system (typically, Netlib's BLAS/LAPACK through arpack). For this reason, dmrimatlab offers several options for the mex code to be linked against. You can choose the one that best fits your needs by creating a file named 'config.octave' in this same sub-folder ('mexcode'). This file will be created from scratch the first time 'makefile_mexcode' is called. It should include variables definition like:

   >> BLAS_CONFIG=blis-local

(or BLAS_CONFIG=netlib, BLAS_CONFIG=openblas, BLAS_CONFIG=openblas-local or BLAS_CONFIG=mkl, see description below) to tell the actual version of BLAS to be used, one of:

- [netlib] So that the system-provided implementation (Netlib's BLAS/LAPACK) will be used. Note that Netlib's implementations are the baseline for performance, since they are very little optimized. Moreover, this option will properly work only if GNU Octave itself is using Netlib's implementations (see note below).

- [openblas] So that a system-wide implementation of OpenBLAS will be used. NOTE: it is assumed that you or your system admin has installed the "openblas" package in your GNU Linux distribution (with pacman, apt, ...). It will provide:

   - /usr/include/openblas/cblas.h, /usr/include/openblas/lapacke.h, /usr/lib/libopenblas.so(or alike).

   OpenBLAS is more efficient than Netlib's implementation, with dedicated kernels for each processor type. Besides, it has multi-threading capabilities. Unfortunately, this is a problem for our mex code, which is also multi-threaded. Though our code prevents this situation by calling "blas_set_num_threads()" as needed, problems can still arise in certain systems with a large number of cores. They can be fixed by running Octave as:

   >> $ OMP_NUM_THREADS=1 octave

   so that OpenBLAS will always run single-thread. Note, however, that any other Octave-related software using Open MP will also run single-thread with this fix, which might impact its performance. Moreover, this option will properly work only if GNU Octave itself is using OpenBLAS as its BLAS/LAPACK implementation (see note below).

- [openblas-local] This will automatically download the source code and compile a local OpenBLAS library without Open MP, so that it will be single-threaded in nature, without conflicting with any other existing libraries. It does not require manually installing any additional software, and will properly work regardless on the BLAS implementation GNU Octave actually uses. However, it is slightly more complex since it has to compile OpenBLAS (and it will take some more time); note that you will need gcc, g++ and gfortran.

- [blis-local] This is a more efficient implementation than OpenBLAS, and it is the default. It will use SHPC's optimized implementation of BLAS, namely BLIS, together with SHPC's FLAME for most of the LAPACK calls. However, since FLAME does not implement each and every existing LAPACK routine, it still relies on Netlib's LAPACK for certain tasks. Both BLIS and FLAME use dedicated kernels for the particular architecture used. The script will automatically download and install these three software pieces without the need of manually installing any additional software, and this option will also work regardless on the BLAS implementation GNU Octave actually uses. Once again, this option is also more complex and will take some time to run (besides, you will need CMake).

- [mkl] In case you are using a x86_64 based processor, the by-far more efficient BLAS/LAPACK implementation is Intel's Oneapi MKL. Using this option is the only way to attain a similar performance as in the Matlab version of the toolbox. Note this requires using your preferred software manager (pacman, apt, ...) to install the package intel-oneapi-mkl (or alike). The installation root for this package will be something like "/opt/intel/mkl", where you should be able to find folders named "lib" and "include". This path must be provided in the config file by setting the variable:

   >> MKL_ROOT=/opt/intel/mkl

   Since MKL depends on libiomp5.so (among others), you must also provide the location of these redistributed libraries, typically something like:

   >> MKL_REDIST=/opt/intel/oneapi/compiler/2025.0/lib

## IMPORTANT NOTE:

GNU Octave is designed in a way that, when mex files are loaded, they preferably use symbols already loaded over those present in the dynamic libraries they are linked against:

- Assume your GNU Octave is using Netlib's BLAS but you have linked your mex files against system-wide OpenBLAS.
- If a mex file calls, for example, dgemv_, you have two implementations for it: one from Netlib's BLAS, one from OpenBLAS.
- Due to GNU Octave's design, the version from Netlib's BLAS will be used because it was loaded earlier (when GNU Octave was started).
- The same thing will happen the other way around, i.e. if you link mex files against Netlib's but GNU Octave is using OpenBLAS.

In the deployed version of the toolbox [https://github.com/atriveg/dmrimatlab-octave-deploy] this is solved by patching the source code of Octave, but this solution does not apply if you have installed Octave with a package manager. Hence, in the non-deployed version:

- Using either [netlib] or [openblas] options is discouraged (unless you ensure consistency).
- If you use [openblas-local] or [blis-local], then the "BUILD_DYNAMIC" variable within config.octave MUST be set "no" (BUILD_DYNAMIC=no). This means that mex files will be statically linked with the respective libraries. It will translate in far larger mex files, but they can be safely used regardless of GNU Octave's configuration.
- If you use [mkl], then the mex files will always be dynamically linked without further issues. This is because Intel Oneapi MKL provides two different interfaces for each routine, e.g. dgemv and dgemv_. By using the former, we avoid "shadowing" symbols with the library GNU Octave is using.
