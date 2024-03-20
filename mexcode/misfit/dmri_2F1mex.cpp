/*==========================================================
 * dmri_2F1mex.c
 *
 * Implements dmri_2F1mex as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "../mathsmex/hypergeom2F1.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

class ThArgs : public DMRIThreader
{
public:
    ElementType g;
    BufferType z;
    BufferType f;
};

THFCNRET dmri_2F1_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    if( (nrhs<2) || (nrhs>3) )
        mexErrMsgIdAndTxt("MyToolbox:dmri_2F1:nrhs","Either 2 or 3 input arguments allowed");
    //=======================================================================================
    if( mxGetNumberOfElements(prhs[0]) != 1 )
         mexErrMsgIdAndTxt("MyToolbox:dmri_2F1:dim","Input 1 must be a scalar");
    ElementType g = mxGetDoubles(prhs[0])[0];
    //=======================================================================================
    SizeType N = mxGetNumberOfElements( prhs[1] );
    //=======================================================================================
    unsigned int maxthreads = 1000000;
    if( nrhs>2 ){
        if( mxGetNumberOfElements(prhs[2]) != 1 )
            mexErrMsgIdAndTxt("MyToolbox:dmri_2F1:dim","Input 3 must be a scalar");
        maxthreads = (unsigned int)( mxGetDoubles(prhs[2])[0] );
    }
    maxthreads = get_number_of_threads( maxthreads );
    //=======================================================================================
    if(nlhs==0)
        return;
    else if(nlhs==1){
        plhs[0] = mxCreateNumericArray(
                mxGetNumberOfDimensions(prhs[1]),
                mxGetDimensions(prhs[1]),
                mxDOUBLE_CLASS,
                mxREAL );
    }
    else
        mexErrMsgIdAndTxt("MyToolbox:dmri_2F1:nlhs","At most 1 output can be returned");
    //=======================================================================================
    unsigned int chunksz = (N/maxthreads)/20;
    if(chunksz<1)
        chunksz = 1;
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, chunksz );
    // Own values:
    threader.g = g;        //
    threader.z = mxGetDoubles(prhs[1]); //
    threader.f = mxGetDoubles(plhs[0]); //
    //=======================================================================================
    threader.threadedProcess( maxthreads, dmri_2F1_process_fcn );
    //=======================================================================================
}

THFCNRET dmri_2F1_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    ElementType g = args->g;
    
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( IndexType i=start; i<end; ++i ){
            ElementType z = args->z[i];
            ElementType val;
            int result = hypegeo::hyperGeom2F1( g, z, val );
            args->f[i] = val;
        }
    }
    while( start < args->getN() );
    
    return (THFCNRET)NULL;
}
