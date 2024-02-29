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

#if defined (_HAS_POSIX_THREADS)
#include <pthread.h>
#include <unistd.h>
#include "omp.h" // Threads management with Lapack/BLAS
#endif

typedef struct{
    unsigned int K;     // (L+1)(L+2)/2;
    unsigned int L;
    unsigned int G;
    unsigned long N;    // The number of voxels to process
    unsigned int* nx;
    unsigned int* ny;
    unsigned int* nz;
    unsigned long* mu;
    double* x;
    double* y;
    double* z;
    BufferType hotin;   // Size N X (L+1)(L+2)/2, where N is the number of voxels. INPUT
    BufferType evals;   // Size N X 1, where N is the number of voxels. OUTPUT
    unsigned int tid;   // Which thread is being run?
    unsigned int nth;   // What is the total number of threads?
} ThArgs;

void* hot2signal_process_fcn( void* );

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
#if defined(_HAS_POSIX_THREADS)
    unsigned int maxthreads = sysconf(_SC_NPROCESSORS_CONF);
    maxthreads = ( (unsigned int)mxGetScalar(prhs[2])<maxthreads ? (unsigned int)mxGetScalar(prhs[2]) : maxthreads );
#else
    unsigned int maxthreads = 1;
#endif
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
    // Put all the information together in the threaded structure
    ThArgs* args = new ThArgs[maxthreads];
    for( unsigned int tid=0; tid<maxthreads; ++tid ){
        args[tid].K = K;
        args[tid].L = L;
        args[tid].G = G;
        args[tid].N = dims[0];
        args[tid].nx = nx;
        args[tid].ny = ny;
        args[tid].nz = nz;
        args[tid].mu = mu;
        args[tid].x = x;
        args[tid].y = y;
        args[tid].z = z;
        args[tid].hotin = mxGetDoubles( prhs[0] );
        args[tid].evals = mxGetDoubles( plhs[0] );
        args[tid].tid = tid;
        args[tid].nth = maxthreads;
    }
    //=======================================================================================
#if defined(_HAS_POSIX_THREADS)
    //=======================================================================================
    pthread_t* threads = new pthread_t[maxthreads];
    int*       rets    = new int[maxthreads];
    for( unsigned int tid=0; tid<maxthreads; ++tid ){
        rets[tid] = pthread_create(
            &(threads[tid]),
            NULL,
            hot2signal_process_fcn, 
            (void*)(&args[tid])    );
    }
    for( unsigned int tid=0; tid<maxthreads; ++tid ){
        pthread_join( threads[tid], NULL);
    }
    delete[] threads;
    delete[] rets;
#else
    hot2signal_process_fcn( (void*)(&args[0]) );
#endif
    //=======================================================================================
    delete[] x;
    delete[] y;
    delete[] z;
    delete[] nx;
    delete[] ny;
    delete[] nz;
    delete[] mu;
    delete[] args;
    //=======================================================================================
    return;
}

void* hot2signal_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    BufferType in = new ElementType[args->K];
    BufferType out = new ElementType[args->G];
#if defined(_HAS_POSIX_THREADS)
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we cab let Open MP do its
    // stuff.
    int omp_mthreads = omp_get_max_threads();
    omp_set_num_threads(1);
#endif
    for( SizeType i=(args->tid); i<(args->N); i+=(args->nth) ){
        // Fill input:
        for( SizeType k1=0; k1<(SizeType)(args->K); ++k1 )
            in[k1] = args->hotin[k1*(args->N)+i];
        // Evaluate:
        sh2hot::evaluateHOT( args->G, args->x, args->y, args->z,
                             out, args->L, args->nx, args->ny, args->nz, args->mu, in );
        // Fill the output:
        for( SizeType k2=0; k2<(SizeType)(args->G); ++k2 )
            args->evals[k2*(args->N)+i] = out[k2];
    }
#if defined(_HAS_POSIX_THREADS)
    omp_set_num_threads(omp_mthreads);
#endif
    delete[] in;
    delete[] out;
    return (void*)NULL;
}
