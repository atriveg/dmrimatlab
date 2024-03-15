/*==========================================================
 * hot2signal_.c
 *
 * Implements the core operations of hot2signal as a mex function
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
    unsigned int L;
    unsigned int G;
    unsigned int* nx;
    unsigned int* ny;
    unsigned int* nz;
    unsigned long* mu;
    double* x;
    double* y;
    double* z;
    BufferType hotin;   // Size N X (L+1)(L+2)/2, where N is the number of voxels. INPUT
    BufferType evals;   // Size N X 1, where N is the number of voxels. OUTPUT
};

THFCNRET hot2signal_process_fcn( void* );

/** The gateway function:
 * INPUTS:
 *     prhs[0]: the input tensor, NxK, K=(L+1)*(L+2)/2, L even
 *     prhs[1]: the gradients table, Gx3, assumed (but not checked) to be unit-norm
 *     prhs[2]: the maximum number of threads to be used, 1x1
 * * OUTPUTS:
 *     plhs[0]: the evaluation of the tensor, NxG
 *     plhs[1]: the multiplicities of the tensor components, Kx1
 *     plhs[2]: the powers of (x,y,z) for each component of the tensor, Kx3
 * */

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:callstack","This function should only be called from hot2signal");
    else if(strcmp(callerFunc,"hot2signal"))
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:callstack","This function should only be called from hot2signal");
    //=======================================================================================
    //=======================================================================================
    if( nrhs!=3 )
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:nrhs","Only 3 input arguments allowed");
    //=======================================================================================
    SizeType ND1  = mxGetNumberOfDimensions( prhs[0] );
    if( ND1!=2 )
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:dim","Input 1 must be 2-D");
    const SizeType* dims = mxGetDimensions( prhs[0] );
    //=======================================================================================
    double Kd = (double)dims[1];
    double Ld = (::sqrt(8*(Kd-1)+9)-3)/2;
    unsigned int L = (unsigned int)Ld;
    if( (Ld-(double)L)>100*mxGetEps() )
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:dim","Weird size for the last dimension of Input 1");
    SizeType K = (L+1)*(L+2)/2;
    //=======================================================================================
    //=======================================================================================
    SizeType ND2  = mxGetNumberOfDimensions( prhs[1] );
    if( ND2!=2 )
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:dim","Input 2 must be 2-D");
    const SizeType* dims2 = mxGetDimensions( prhs[1] );
    if(dims2[1]!=3)
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:dim","Input 2 must have size Gx3");
    unsigned int G = dims2[0];
    //=======================================================================================
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[2]) );
    //=======================================================================================
    //=======================================================================================
    if( nlhs>3 )
        mexErrMsgIdAndTxt("MyToolbox:hot2signal_:nlhs","At most 3 outputs can be returned");
    if( nlhs==0 )
        return;
    //=======================================================================================
    //=======================================================================================
    unsigned int* nx = new unsigned int[K];
    unsigned int* ny = new unsigned int[K];
    unsigned int* nz = new unsigned int[K];
    sh2hot::computeHOTPowers( L, nx, ny, nz );
    unsigned long* mu = new unsigned long[K];
    sh2hot::computeHOTMultiplicity( L, mu, nx, ny, nz );
    //=======================================================================================
    SizeType* dims_ = new SizeType[2];
    dims_[0] = dims[0];
    dims_[1] = G;
    plhs[0] = mxCreateNumericArray( ND1, dims_, mxDOUBLE_CLASS, mxREAL );
    delete[] dims_;
    //=======================================================================================
    if(nlhs>1){
        SizeType* dims2 = new SizeType[2];
        dims2[0] = K;
        dims2[1] = 1;
        plhs[1] = mxCreateNumericArray( 2, dims2, mxDOUBLE_CLASS, mxREAL );
        delete[] dims2;
        for( unsigned int k=0; k<K; ++k )
            mxGetDoubles(plhs[1])[k] = (ElementType)(mu[k]);
    }
    //=======================================================================================
    if(nlhs>2){
        SizeType* dims2 = new SizeType[2];
        dims2[0] = K;
        dims2[1] = 3;
        plhs[2] = mxCreateNumericArray( 2, dims2, mxDOUBLE_CLASS, mxREAL );
        delete[] dims2;
        for( unsigned int k=0; k<K; ++k ){
            mxGetDoubles(plhs[2])[k]     = (ElementType)(nx[k]);
            mxGetDoubles(plhs[2])[k+K]   = (ElementType)(ny[k]);
            mxGetDoubles(plhs[2])[k+2*K] = (ElementType)(nz[k]);
        }
    }
    //=======================================================================================
    double* x = new double[G];
    double* y = new double[G];
    double* z = new double[G];
    for( unsigned int g=0; g<G; ++g ){
        x[g] = mxGetDoubles(prhs[1])[g];
        y[g] = mxGetDoubles(prhs[1])[g+G];
        z[g] = mxGetDoubles(prhs[1])[g+2*G];
    }
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( dims[0], 20 );
    // Own values:
    threader.K = K;
    threader.L = L;
    threader.G = G;
    threader.nx = nx;
    threader.ny = ny;
    threader.nz = nz;
    threader.mu = mu;
    threader.x = x;
    threader.y = y;
    threader.z = z;
    threader.hotin = mxGetDoubles( prhs[0] );
    threader.evals = mxGetDoubles( plhs[0] );
    //=======================================================================================
    threader.threadedProcess( maxthreads, hot2signal_process_fcn );
    //=======================================================================================
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] nx;
    delete[] ny;
    delete[] nz;
    delete[] mu;
    //=======================================================================================
    return;
}

THFCNRET hot2signal_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    SizeType N = args->getN();
    BufferType in = new ElementType[args->K];
    BufferType out = new ElementType[args->G];

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
            for( SizeType k1=0; k1<(SizeType)(args->K); ++k1 )
                in[k1] = args->hotin[k1*N+i];
            // Evaluate:
            sh2hot::evaluateHOT( args->G, args->x, args->y, args->z,
                                out, args->L, args->nx, args->ny, args->nz, args->mu, in );
            // Fill the output:
            for( SizeType k2=0; k2<(SizeType)(args->G); ++k2 )
                args->evals[k2*N+i] = out[k2];
        }
    }
    while( start < N );

    blas_num_threads(blas_threads);
    
    delete[] in;
    delete[] out;
    return (THFCNRET)NULL;
}
