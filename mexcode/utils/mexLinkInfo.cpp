/*==========================================================
 * mexLinkInfo.cpp
 *
 * This function simply prints on screen the name of the
 * BLAS/LAPACK libraries mex files are linked against.
 * Only for GNU Octave builds under GNU Linux
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "../mathsmex/matrixCalculus.h"
#include <iostream>
#include <dlfcn.h>

#ifndef OCTAVE_BUILD
#error "This mex file should be compiled only for GNU Octave builds"
#endif

void check_symbol_library( const void*, const char*, const char* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    const void* test_symbol1 = reinterpret_cast<const void*>(BLASCALLFCN(dgemm));
    const void* test_symbol2 = reinterpret_cast<const void*>(LAPACKCALLFCN(dposv));
    const char* library1 = "BLAS";
    const char* library2 = "LAPACK";
    const char* basename = "mexLinkInfo.mex";

    check_symbol_library( test_symbol1, library1, basename );
    check_symbol_library( test_symbol2, library2, basename );
    
    return;
}

void check_symbol_library( const void* symbol, const char* library, const char* basename )
{
    Dl_info info;

    if( dladdr( symbol, &info ) ){
        if( info.dli_fname != NULL ){
            if( info.dli_sname == NULL )
                std::cout << "   " << library << ": symbols statically linked (check the value of BLAS_CONFIG in config.octave)" << std::endl;
            else{
                if(   (strlen(info.dli_fname) >= strlen(basename)) && (strcmp(basename,&(info.dli_fname)[strlen(info.dli_fname)-strlen(basename)]) == 0)   )
                    std::cout << "   " << library << ": symbols statically linked (check the value of BLAS_CONFIG in config.octave)" << std::endl;
                else
                    std::cout << "   " << library << ": " << info.dli_fname << std::endl;
            }
        }
        else{
            std::cerr << "   Could not determine " << library << " implementation for mex files" << std::endl;
        }
    }
    else{
        std::cerr << "    Could not resolve the origin of " << library << " calls in mex files" << std::endl;
    }

    return;
}
