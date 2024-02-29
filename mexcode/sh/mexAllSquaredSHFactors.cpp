/*==========================================================
 * mexAllSquaredSHFactors.c
 *
 * Implements triple SH real products
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "../mathsmex/sphericalHarmonics.h"
#include <string.h>
#include "../mathsmex/mexToMathsTypes.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if(nrhs!=1){
        mexErrMsgIdAndTxt("MyToolbox:mexAllSquaredSHFactors:nrhs","Only one input (L) admitted.");
    }
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ){
        mexErrMsgIdAndTxt("MyToolbox:mexAllSquaredSHFactors:notDouble","Input must be type double.");
    }
    if( mxGetNumberOfElements(prhs[0]) != 1 ){
        mexErrMsgIdAndTxt("MyToolbox:mexAllSquaredSHFactors:scalar","Input must be a scalar.");
    }
    double Ld = mxGetScalar(prhs[0]);
    if( Ld<0 ){
        mexErrMsgIdAndTxt("MyToolbox:mexAllSquaredSHFactors:sign","Input must be positive.");
    }
    unsigned int L = (unsigned int)Ld;
    if( (L-(double)Ld) > 10*mxGetEps() ){
        mexErrMsgIdAndTxt("MyToolbox:mexAllSquaredSHFactors:int","Input must be an integer.");
    }
    if( L != 2*(L/2) ){
        mexErrMsgIdAndTxt("MyToolbox:mexAllSquaredSHFactors:even","Input must be an even integer.");
    }
    
    SizeType N = shmaths::computeNumberOfSquaredSHFactors( L );
    
    BufferType factors = new ElementType[N];
    unsigned int*  l1 = new unsigned int[N];
    unsigned int*  l2 = new unsigned int[N];
    unsigned int*  l3 = new unsigned int[N];
    int*           m1 = new int[N];
    int*           m2 = new int[N];
    int*           m3 = new int[N];
    
    shmaths::computeSquaredSHFactors( L, factors, l1, l2, l3, m1, m2, m3 );
    
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( 1, N, mxREAL );
        memcpy( mxGetDoubles(plhs[0]), factors, N*sizeof(ElementType) );
    }
    if(nlhs>1){
        plhs[1] = mxCreateDoubleMatrix( 1, N, mxREAL );
        for( SizeType n=0; n<N; ++n )
            mxGetDoubles(plhs[1])[n] = (double)(l1[n]);
    }
    if(nlhs>2){
        plhs[2] = mxCreateDoubleMatrix( 1, N, mxREAL );
        for( SizeType n=0; n<N; ++n )
            mxGetDoubles(plhs[2])[n] = (double)(l2[n]);
    }
    if(nlhs>3){
        plhs[3] = mxCreateDoubleMatrix( 1, N, mxREAL );
        for( SizeType n=0; n<N; ++n )
            mxGetDoubles(plhs[3])[n] = (double)(l3[n]);
    }
    if(nlhs>4){
        plhs[4] = mxCreateDoubleMatrix( 1, N, mxREAL );
        for( SizeType n=0; n<N; ++n )
            mxGetDoubles(plhs[4])[n] = (double)(m1[n]);
    }
    if(nlhs>5){
        plhs[5] = mxCreateDoubleMatrix( 1, N, mxREAL );
        for( SizeType n=0; n<N; ++n )
            mxGetDoubles(plhs[5])[n] = (double)(m2[n]);
    }
    if(nlhs>6){
        plhs[6] = mxCreateDoubleMatrix( 1, N, mxREAL );
        for( SizeType n=0; n<N; ++n )
            mxGetDoubles(plhs[6])[n] = (double)(m3[n]);
    }
    delete[] factors;
    delete[] l1;
    delete[] l2;
    delete[] l3;
    delete[] m1;
    delete[] m2;
    delete[] m3;
}
