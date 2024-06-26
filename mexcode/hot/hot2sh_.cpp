/*==========================================================
 * hot2sh_.c
 *
 * Implements the core operations of hot2sh as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "../mathsmex/sphericalHarmonics.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/sh2hot.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

class ThArgs : public DMRIThreader
{
public:
    unsigned int K;     // (L+1)(L+2)/2;
    BufferType A;       // The values of the conversion matrix
    BufferType hotin;   // Size M X (L+1)(L+2)/2, where M is the number of voxels. INPUT
    BufferType shout;   // Size M X (L+1)(L+2)/2, where M is the number of voxels. OUTPUT
};

THFCNRET hot2sh_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:hot2sh_:callstack","This function should only be called from hot2sh");
    else if(strcmp(callerFunc,"hot2sh"))
        mexErrMsgIdAndTxt("MyToolbox:hot2sh_:callstack","This function should only be called from hot2sh");
    //=======================================================================================
    if( nrhs!=2 )
        mexErrMsgIdAndTxt("MyToolbox:hot2sh_:nrhs","Only 2 input arguments allowed");
    SizeType ND1  = mxGetNumberOfDimensions( prhs[0] );
    if( ND1!=2 )
        mexErrMsgIdAndTxt("MyToolbox:hot2sh_:dim","Input 2 must be 2-D");
    const SizeType* dims = mxGetDimensions( prhs[0] );
    //=======================================================================================
    double Kd = (double)dims[1];
    double Ld = (::sqrt(8*(Kd-1)+9)-3)/2;
    unsigned int L = (unsigned int)Ld;
    if( (Ld-(double)L)>100*mxGetEps() )
        mexErrMsgIdAndTxt("MyToolbox:hot2sh_:dim","Weird size for the last dimension of Input 1");
    //=======================================================================================
    SizeType K = (L+1)*(L+2)/2;
    unsigned long dim = sh2hot::sh2hotHardcodesDim( L );
    if(dim==0)
        mexErrMsgIdAndTxt("MyToolbox:hot2sh_:L","The maximum order L of the SH is too large. You should re-generate sh2hothardcodes.cxx");
    //=======================================================================================
    if( nlhs>2 )
        mexErrMsgIdAndTxt("MyToolbox:hot2sh_:nlhs","Only 1 or 2 outputs can be returned");
    if( nlhs==0 )
        return;
    //=======================================================================================
    plhs[0] = mxCreateNumericArray( ND1, dims, mxDOUBLE_CLASS, mxREAL );
    //=======================================================================================
    BufferType A = (BufferType)NULL;
    if( nlhs<2 ) // No need to return the conversion matrix
        A = new ElementType[ dim ];
    else{
        SizeType* dims2 = new SizeType[2];
        dims2[0] = K;
        dims2[1] = K;
        plhs[1] = mxCreateNumericArray( 2, dims2, mxDOUBLE_CLASS, mxREAL );
        delete[] dims2;
        A = mxGetDoubles( plhs[1] );
    }
    sh2hot::sh2hotHardcodes(L,A);
    //==========================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[1]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( dims[0], 20 );
    // Own values:
    threader.K = K;
    threader.A = A;
    threader.hotin = mxGetDoubles( prhs[0] );
    threader.shout = mxGetDoubles( plhs[0] );
    //=======================================================================================
    threader.threadedProcess( maxthreads, hot2sh_process_fcn );
    //=======================================================================================
    if( nlhs<2 )
        delete[] A;
    //=======================================================================================
    return;
}

THFCNRET hot2sh_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    SizeType N = args->getN();
    BufferType in  = new ElementType[args->K];
    BufferType out = new ElementType[args->K];

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads(1);
    
    // ---------------------------------------------------------------
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( IndexType i=start; i<end; ++i ){
            // Fill input:
            for( SizeType k1=0; k1<args->K; ++k1 )
                in[k1] = args->hotin[k1*N+i];
            // Matrix product:
            mataux::multiplyMxArrays( args->A, in, out, args->K, args->K, 1 );
            // Fill output:
            for( SizeType k2=0; k2<args->K; ++k2 )
                args->shout[k2*N+i] = out[k2];
        }
    }
    while( start < N );

    blas_num_threads(blas_threads);
    
    delete[] in;
    delete[] out;
    return (THFCNRET)NULL;
}
