# DMRIMatlab
This is a complete and self-contained toolbox for basic-to-advanced diffusion MRI reconstruction, processing, and representation (see complete description in: http://www.lpi.tel.uva.es/dmrilab). It works with both Matlab (multi-platform) and Octave (Linux, by now).

It includes:

- The classical Diffusion Tensor Imaging (DTI) formalism with linear and non-linear signal fitting.
- High Angular Resolution Diffusion Imaging (HARDI):
        - Spherical Harmonics (SH).
        - Orientation Probability Density Transform (OPDT).
        - Higher Order Tensors (HOT).
- (Constrained, regularized) Diffusion Spectrum Imaging (DSI).
- Spherical Means and Spherical Deconvolution (MiSFIT).
- (Constrained, regularized) Mean Apparent Propagator MRI (MAP-MRI/MAPL).
- Free-water elimination.
- Denoising of Diffusion Weighted Images.
- Representation of scalar maps, color-by-orientation, DTI glyphs, or ODF fields.
- Tractography based on DTI and global geodesics (Finsler tractography).

The toolbox is designed pursuing computational performance, so that large databases can be processed within a reasonable time to perform group studies, connectomics, atlasing or IA training afterwards. This is attained by coding the key parts of the algorithms as (efficiently) multi-threaded C/C++ mex code, and exploiting advanced BLAS/LAPACK optimizations when possible.

Authors:

        Antonio Tristán Vega, Santiago Aja-Fernández, Guillem París, Tomasz Pieciak.
        Laboratorio de Procesado de Imagen. Universidad de Valladolid
        Spain

In case you use the package for your own research, we ask you to kindly cite it as:

        Antonio Tristán-Vega, Santiago Aja-Fernández, Guillem París and Tomasz Pieciak. "dMRI-Lab: advanced diffusion MRI with Matlab" [Online resource] https://www.lpi.tel.uva.es/dmrilab September 2025. Universidad de Valladolid. Spain.

## Getting started:

1. Download the code from https://github.com/atriveg/dmrimatlab.
2. From the Matlab/Octave command window, cd to the home folder, i.e. that containing the setup script "setup__DMRIMatlab_toolbox.m".
3. Run the setup script as either (the 'useparallel' option has no effect in Octave by now):

        >> setup__DMRIMatlab_toolbox('useparallel',true);
        >> setup__DMRIMatlab_toolbox('useparallel',false);

   or simply:

        >> setup__DMRIMatlab_toolbox;

   for using/avoid using the Parallel Computing Toolbox. In case you don't have a working license for it, the script will not throw an error. This will setup your Matlab/Octave path for the present session (it won't make any permanent changes).
4. Download the test data used for demo files and tests by running [YOU ONLY NEED TO (SUCCESSFULLY) RUN THIS COMMAND ONCE]:

        >> download_dmritestdata

5. Make sure you have a suitable C/C++ compiler installed in your computer, preferably the GCC suite (in Windows, we recommend installing the "MinGW" Add-On from the Matlab interface). In Octave, you might need a Fortran compiler (gfortran) as well. In Matlab, type:

        >> mex-setup
        >> mex-setup C++

In Octave:

        >> mex --print CC
        >> mex --print CXX

6. Run the script [YOU ONLY NEED TO RUN THIS COMMAND EACH TIME YOU PULL A NEW VERSION FROM THE REPO]:

        >> makefile_mexcode

    which will build all the necessary mex files within the toolbox [FOR OCTAVE BUILDS, PLEASE CHECK THE README.octave FILE IN THE 'mexcode' SUBFOLDER].
7. Open/run some of the demo files from the launch menu by typing:

        >> dmrimatlab_demo

    In Matlab, .mlx notebook files are used. In Octave, they are .ipynb noteboks [IN OCTAVE YOU WILL ADDITIONALLY NEED JupyterLab WITH THE KERNEL FOR OCTAVE].

## A note on Octave's computational performance

Unlike Matlab, GNU Octave does not pack any particular implementations of BLAS/LAPACK, but instead it uses those provided by the system (usually Netlib's libblas and liblapack through arpack). Note this may dramatically decrease its computational performance compared to Matlab. To avoid this issue, you may install OpenBLAS or, even better for compatible hardware, Intel's MKL and force octave to use them by doing something like:

        $ LD_PRELOAD="/usr/lib/libopenblas.so:${LD_PRELOAD}" octave

or:

        $ LD_PRELOAD="/opt/intel/mkl/lib/intel64/libmkl_rt.so:${LD_PRELOAD}" octave

You may even write your own wrapper for GNU Octave, something like:

        /usr/bin/octave-mkl:
                #!/bin/bash
                LD_PRELOAD="/opt/intel/mkl/lib/intel64/libmkl_rt.so:${LD_PRELOAD}" /usr/bin/octave "$@"

## Additional software required

The toolbox is designed to be self-contained, so that core methods do not rely in any external software not directly provided. It does not depend either on any particular Matlab's or Octave's toolbox/package. Note, however, that developer tools (GCC suite or alike) are required to compile the mex functions of the toolbox. Depending on your system configuration, in Octave builds this might imply the need of installing certain packages (mainly: cblas, lapacke, openblas or others, see the README.octave file in the 'mexcode' subfolder).

Yet, test programs make use of certain specific functions from toolboxes/packages that you will need to install in case you want to actually run these tests. These are:

- For Matlab: Optimization, Signal Processing, Statistics, and Symbolic Maths.
- For Octave: gsl, optim, signal, statistics.

NOTE: in some platforms it is likely that you get a compile error when trying to install the optim package. In this case you will have to:
        1. Download the source code of the package, optim-1.6.2.tar.gz, from: https://gnu-octave.github.io/packages/optim/
        2. Untar and look for the file named: src/__max_nargin_optim__.cc
        3. Patch this file according to: https://sourceforge.net/p/octave/optim/ci/d8c28ab3f3f37d439bc2d961d79f4a8caa9830d1/ i.e. change:
               (fcn.user_function_value ()->parameter_list ()->length ());
           to:
               (fcn.user_function_value ()->parameter_list ()->size ());
        4. Pack again the folder to optim-1.6.2.tar.gz
        5. Within Octave, browse to the folder where your resulting optim-1.6.2.tar.gz is, then run <pkg install optim-1.6.2.tar.gz>

Finally, Octave's demos are written as Jupyter notebooks, so that you will need JupyterLab and Jupyter's kernel for Octave to run them. NOTE: since JupyterLab uses by default the octave-cli program, it is very convenient that you run it like this:

        $ OCTAVE_EXECUTABLE=/usr/bin/octave jupyter-lab

so that it actually uses the full-featured Octave program and you can enjoy Qt-based graphics instead of the more limited gnuplot (this will be automatically done for you in case you run the demos from the launch menu by using the dmrimatlab_demo command).
