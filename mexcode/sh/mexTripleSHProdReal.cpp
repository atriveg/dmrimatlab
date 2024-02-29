/*==========================================================
 * mexTripleSHProdReal.c
 *
 * Implements triple SH real products
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "../mathsmex/sphericalHarmonics.h"
#include "../mathsmex/mexToMathsTypes.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double indicesD[6];
    int    indices[6];
    if(nrhs==1){
        if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ){
            mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:notDouble","Inputs must be type double.");
        }    
        if( mxGetNumberOfElements(prhs[0]) != 6 ){
            mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:nrhs","If one argument is provided, it must be a 6 element vector.");
        }
        BufferType args = mxGetDoubles(prhs[0]);
        for(unsigned int n=0; n<6; ++n )
            indicesD[n] = args[n];
    }
    else if(nrhs==6){
        for( unsigned int n=0; n<6; ++n ){
            if( !mxIsDouble(prhs[n]) || mxIsComplex(prhs[n]) ){
                mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:notDouble","Inputs must be type double.");
            }
            if( mxGetNumberOfElements(prhs[n]) != 1 ){
                mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:nrhs","If six arguments are provided, each one must be a scalar.");
            }
            indicesD[n] = mxGetScalar(prhs[n]);
        }
    }
    else{
        mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:nrhs","Either one or six arguments allowed.");
    }
    for( unsigned int n=0; n<3; ++n ){
        if( indicesD[n]<0 ){
            mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:lpos","l-indices must be positive.");
        }
        indices[n] = (int)indicesD[n];
        double err = indicesD[n] - (double)indices[n];
        err = ( err>0 ? err : -err );
        if(err>10*mxGetEps()){
            mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:notInteger","l-indices must be integer.");
        }
        if( !shmaths::isEven(indices[n]) ){
            mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:nonEven","l-indices must be even.");
        }
    }
    for( unsigned int n=3; n<6; ++n ){
        indices[n] = (int)indicesD[n];
        double err = indicesD[n] - (double)indices[n];
        err = ( err>0 ? err : -err );
        if(err>10*mxGetEps()){
            mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:notInteger","m-indices must be integer.");
        }
    }    
    bool nonnull;
    bool resnan;
    unsigned long NF = indices[0]+indices[1]+indices[2]+1;
    double* factorials = new double[NF+1];
    shmaths::cumprod( NF, factorials );
    double res = shmaths::computeTripleSHProd(
        (unsigned int)indices[0],  (unsigned int)indices[1], (unsigned int)indices[2],
        indices[3], indices[4], indices[5], nonnull, resnan, factorials );
    delete[] factorials;
    if(resnan){
        mexErrMsgIdAndTxt("MyToolbox:mexTripleSHProdReal:nanvalue","The requested combination of l and m is undefined.");
    }
    
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL );
        BufferType output = mxGetDoubles(plhs[0]);
        output[0] = res;
    }
}
