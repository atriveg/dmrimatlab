/*==========================================================
 * sh2squaredsh_.c
 *
 * Implements sh2squaredsh_ as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
//#include "matrix.h"
#include "math.h"
#include "../mathsmex/sphericalHarmonics.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

class ThArgs : public DMRIThreader
{
public:
    SizeType P;         // The number of non-null factors to compute squared SH from SH
    IndexType* k1;      // Size P, positions of the first input SH coefficient
    IndexType* k2;      // Size P, positions of the second input SH coefficient
    IndexType* k3;      // Size P, positions of the output SH coefficient
    BufferType factors; // Size P, components of the quadratic form relating input and output SH coeffs.
    BufferType shin;    // Size M X (L+1)(L+2)/2, where M is the number of voxels. INPUT
    BufferType shout;   // Size M X (2L+1)(2L+2)/2, where M is the number of voxels. OUTPUT
    SizeType fov;       // M, the total number of voxels
    SizeType procs;     // The number of voxels processed by this thread
    SizeType K1;        // (L+1)(L+2)/2, the number of input SH coefficients
    SizeType K2;        // (2L+1)(2L+2)/2, the number of output SH coefficients
};

THFCNRET sh2sqsh_process_fcn( void* );

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
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh_:callstack","This function should only be called from sh2squaredsh");
    else if( strcmp(callerFunc,"sh2squaredsh") )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh_:callstack","This function should only be called from sh2squaredsh");
    //=======================================================================================
    if( (nrhs<1) || (nrhs>2) )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh_:nrhs","Either 1 or 2 input arguments allowed");
    SizeType ND1  = mxGetNumberOfDimensions( prhs[0] );
    if( ND1<2 )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh_:dim","Input 1 must be at least 2-D");
    //=======================================================================================
    SizeType fov  = 1;
    const SizeType* dims = mxGetDimensions( prhs[0] );
    for( SizeType d=0; d<ND1-1; ++d )
        fov *= dims[d];
    //=======================================================================================
    double K = (double)dims[ND1-1];
    double Ld = (::sqrt(8*(K-1)+9)-3)/2;
    unsigned int L = (unsigned int)Ld;
    if( (Ld-(double)L)>100*mxGetEps() )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh_:dim","Weird size for the last dimension of Input 1");
    //=======================================================================================
    SizeType Kin  = (L+1)*(L+2)/2;
    SizeType Kout = (2*L+1)*(2*L+2)/2;
    //=======================================================================================
    // Check the number of threads to use
    unsigned int maxthreads = 1000000;
    if( nrhs>1 ){
        if( mxGetNumberOfElements(prhs[1]) != 1 ){
            mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh_:dim","Input 2 should be a scalar integer");
        }
        maxthreads = (unsigned int)( mxGetDoubles( prhs[1] )[0] );
    }
    maxthreads = get_number_of_threads( maxthreads );
    //=======================================================================================
    if( nlhs>1 )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh_:nlhs","Only 1 output can be returned");
    if( nlhs==0 )
        return;
    SizeType* dims2 = new SizeType[ND1];
    for( unsigned int d=0; d<ND1-1; ++d )
        dims2[d] = dims[d];
    dims2[ND1-1] = Kout;
    plhs[0] = mxCreateNumericArray( ND1, dims2, mxDOUBLE_CLASS, mxREAL );
    delete[] dims2;
    //=======================================================================================
    SizeType P = shmaths::computeNumberOfSquaredSHFactors( L );
    BufferType factors = new ElementType[P];
    unsigned int*  l1 = new unsigned int[P];
    unsigned int*  l2 = new unsigned int[P];
    unsigned int*  l3 = new unsigned int[P];
    int*           m1 = new int[P];
    int*           m2 = new int[P];
    int*           m3 = new int[P];
    shmaths::computeSquaredSHFactors( L, factors, l1, l2, l3, m1, m2, m3 );
    IndexType* k1 = new IndexType[P];
    IndexType* k2 = new IndexType[P];
    IndexType* k3 = new IndexType[P];
    shmaths::unrollEvenSHIndices( P, l1, m1, k1 );
    shmaths::unrollEvenSHIndices( P, l2, m2, k2 );
    shmaths::unrollEvenSHIndices( P, l3, m3, k3 );
    delete[] l1;
    delete[] l2;
    delete[] l3;
    delete[] m1;
    delete[] m2;
    delete[] m3;
    //=======================================================================================
    BufferType shin  = mxGetDoubles( prhs[0] );
    BufferType shout = mxGetDoubles( plhs[0] );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( fov, 20 );
    // Own values:
    threader.P       = P;        //
    threader.k1      = k1;       //
    threader.k2      = k2;       //
    threader.k3      = k3;       //
    threader.factors = factors;  //
    threader.shin    = shin;     //
    threader.shout   = shout;    //
    threader.fov     = fov;      //
    threader.procs   = 0;        //
    threader.K1      = Kin;      //
    threader.K2      = Kout;     //
    //=======================================================================================
    threader.threadedProcess( maxthreads, sh2sqsh_process_fcn );
    //=======================================================================================
    delete[] factors;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    //=======================================================================================
}

THFCNRET sh2sqsh_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    BufferType in  = new ElementType[args->K1];
    BufferType out = new ElementType[args->K2];
    SizeType fov = args->fov;

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads_thread(1);
    
    // ---------------------------------------------------------------
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( IndexType i=start; i<end; ++i ){
            for( SizeType k1=0; k1<args->K1; ++k1 )
                in[k1] = args->shin[k1*fov+i];
            for( SizeType k3=0; k3<args->K2; ++k3 )
                out[k3] = 0.0f;
            for( SizeType p=0; p<args->P; ++p )
                out[args->k3[p]] += in[args->k1[p]]*in[args->k2[p]]*(args->factors[p]);
            for( SizeType k3=0; k3<args->K2; ++k3 )
                args->shout[k3*fov+i] = out[k3];
            (args->procs)++;
        }
    }
    while( start < args->getN() );

    blas_num_threads_thread(blas_threads);
    
    delete[] in;
    delete[] out;
    return (THFCNRET)NULL;
}
