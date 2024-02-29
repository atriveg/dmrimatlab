# dmrimatlab
This is complete Matlab toolbox (including multi-threaded C/C++-mex code) for diffusion MRI processing. See complete description in: http://www.lpi.tel.uva.es/dmrilab

Authors:

&nbsp;&nbsp;&nbsp;&nbsp; Antonio Tristán Vega, Santiago Aja-Fernández, Guillem París.

&nbsp;&nbsp;&nbsp;&nbsp; Laboratorio de Procesado de Imagen. Universidad de Valladolid

&nbsp;&nbsp;&nbsp;&nbsp; Spain

In case you use the package for your own research, we ask you to kindly cite it as:

&nbsp;&nbsp;&nbsp;&nbsp; Antonio Tristán-Vega, Santiago Aja-Fernández and Guillem París. "dMRI-Lab: advanced diffusion MRI with Matlab" [Online resource] https://www.lpi.tel.uva.es/dmrilab January 2022. Universidad de Valladolid. Spain

## Getting started:

1. Download the code.
2. From the Matlab command window, cd to the home folder, i.e. that containing the setup script "setup__DMRIMatlab_toolbox.m".
3. Run the setup script as either:

        >> setup\_\_DMRIMatlab\_toolbox('useparallel',true);
        >> setup\_\_DMRIMatlab\_toolbox('useparallel',false); or simply: >> setup\_\_DMRIMatlab\_toolbox;

    for using/avoid using the Parallel Computing Toolbox. In case you don't have a working license for it, the script will not throw an error. This will setup your Matlab path for the present session (it won't make any permanent changes).
4. Make sure you have a suitable C/C++ compiler installed in your computer (in Windows, we recommend installing the "MinGW" Add-On from the Matlab interface).
5. Run the script:

        >> makefile\_mexcode

    which will build all the necessary mex files within the toolbox.
6. Open/run some of the demo files in the "examples" folder to get started.



