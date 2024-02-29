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

#if defined (_HAS_POSIX_THREADS)
#include <pthread.h>
#include <unistd.h>
#include "omp.h" // Threads management with Lapack/BLAS
#endif

typedef struct{
    unsigned int K;     // (L+1)(L+2)/2;
    unsigned long N;    // The number of voxels to process
    BufferType A;       // The values of the conversion matrix
    BufferType hotin;   // Size M X (L+1)(L+2)/2, where M is the number of voxels. INPUT
    BufferType shout;   // Size M X (L+1)(L+2)/2, where M is the number of voxels. OUTPUT
    unsigned int tid;   // Which thread is being run?
    unsigned int nth;   // What is the total number of threads?
} ThArgs;

void* hot2sh_process_fcn( void* );

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
#if defined(_HAS_POSIX_THREADS)
    unsigned int maxthreads = sysconf(_SC_NPROCESSORS_CONF);
    maxthreads = ( (unsigned int)mxGetScalar(prhs[1])<maxthreads ? (unsigned int)mxGetScalar(prhs[1]) : maxthreads );
#else
    unsigned int maxthreads = 1;
#endif
    //=======================================================================================
    // Put all the information together in the threaded structure
    ThArgs* args = new ThArgs[maxthreads];
    for( unsigned int tid=0; tid<maxthreads; ++tid ){
        args[tid].K = K;
        args[tid].N = dims[0];
        args[tid].A = A;
        args[tid].hotin = mxGetDoubles( prhs[0] );
        args[tid].shout = mxGetDoubles( plhs[0] );
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
            hot2sh_process_fcn, 
            (void*)(&args[tid])    );
    }
    for( unsigned int tid=0; tid<maxthreads; ++tid ){
        pthread_join( threads[tid], NULL);
    }
    delete[] threads;
    delete[] rets;
#else
    hot2sh_process_fcn( (void*)(&args[0]) );
#endif
    //=======================================================================================
    if( nlhs<2 )
        delete[] A;
    delete[] args;
    //=======================================================================================
    return;
}

void* hot2sh_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    BufferType in  = new ElementType[args->K];
    BufferType out = new ElementType[args->K];
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
        for( SizeType k1=0; k1<args->K; ++k1 )
            in[k1] = args->hotin[k1*(args->N)+i];
        // Matrix product:
        mataux::multiplyMxArrays( args->A, in, out, args->K, args->K, 1 );
        // Fill output:
        for( SizeType k2=0; k2<args->K; ++k2 )
            args->shout[k2*(args->N)+i] = out[k2];
    }
#if defined(_HAS_POSIX_THREADS)
    omp_set_num_threads(omp_mthreads);
#endif
    delete[] in;
    delete[] out;
    return (void*)NULL;
}
