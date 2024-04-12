# dmrimatlab
This is a complete Matlab toolbox (including multi-threaded C/C++-mex code) for diffusion MRI processing (see complete description in: http://www.lpi.tel.uva.es/dmrilab).

It includes:

- The classical Diffusion Tensor Imaging (DTI) formalism with linear and non-linear signal fitting.
- High Angular Resolution Diffusion Imaging (HARDI):
        - Spherical Harmonics (SH).
        - Orientation Probability Density Transform (OPDT).
        - Higher Order Tensors (HOT).
- Diffusion Spectrum Imaging (DSI).
- Spherical Means and Spherical Deconvolution (MiSFIT).
- Mean Apparent Propagator MRI (MAP-MRI/MAPL).
- Free-water elimination.
- Denoising of Diffusion Weighted Images.
- Representation of scalar maps, color-by-orientation, DTI glyphs, or ODF fields.

Authors:

&nbsp;&nbsp;&nbsp;&nbsp; Antonio Tristán Vega, Santiago Aja-Fernández, Guillem París.

&nbsp;&nbsp;&nbsp;&nbsp; Laboratorio de Procesado de Imagen. Universidad de Valladolid

&nbsp;&nbsp;&nbsp;&nbsp; Spain

In case you use the package for your own research, we ask you to kindly cite it as:

&nbsp;&nbsp;&nbsp;&nbsp; Antonio Tristán-Vega, Santiago Aja-Fernández and Guillem París. "dMRI-Lab: advanced diffusion MRI with Matlab" [Online resource] https://www.lpi.tel.uva.es/dmrilab January 2022. Universidad de Valladolid. Spain

## Getting started:

1. Download the code from https://github.com/atriveg/dmrimatlab.
2. From the Matlab command window, cd to the home folder, i.e. that containing the setup script "setup__DMRIMatlab_toolbox.m".
3. Run the setup script as either:

        >> setup__DMRIMatlab_toolbox('useparallel',true);
        >> setup__DMRIMatlab_toolbox('useparallel',false);

   or simply:

        >> setup__DMRIMatlab_toolbox;

   for using/avoid using the Parallel Computing Toolbox. In case you don't have a working license for it, the script will not throw an error. This will setup your Matlab path for the present session (it won't make any permanent changes).
4. Download the test data used for demo files and tests by running [YOU ONLY NEED TO (SUCCESSFULLY) RUN THIS COMMAND ONCE]:

        >> download_dmritestdata

5. Make sure you have a suitable C/C++ compiler installed in your computer (in Windows, we recommend installing the "MinGW" Add-On from the Matlab interface).
6. Run the script [YOU ONLY NEED TO RUN THIS COMMAND EACH TIME YOU PULL A NEW VERSION FROM THE REPO]:

        >> makefile_mexcode

    which will build all the necessary mex files within the toolbox.
7. Open/run some of the demo files in the "examples" folder to get started.



