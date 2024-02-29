/*==========================================================
 * sh2squaredsh.c
 *
 * Implements sh2squaredsh as a mex function
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
#include "../mathsmex/mexToMathsTypes.h"

#if defined(_HAS_POSIX_THREADS)
#include <pthread.h>
#include <unistd.h>
#include "omp.h" // Threads management with Lapack/BLAS
#endif

typedef struct{
    SizeType N;         // The number of non-null factors to compute squared SH from SH
    IndexType* k1;      // Size N, positions of the first input SH coefficient
    IndexType* k2;      // Size N, positions of the second input SH coefficient
    IndexType* k3;      // Size N, positions of the output SH coefficient
    BufferType factors; // Size N, components of the quadratic form relating input and output SH coeffs.
    BufferType shin;    // Size M X (L+1)(L+2)/2, where M is the number of voxels. INPUT
    BufferType shout;   // Size M X (2L+1)(2L+2)/2, where M is the number of voxels. OUTPUT
    SizeType fov;       // M, the total number of masked voxels
    unsigned int tid;   // Which thread is being run?
    unsigned int nth;   // What is the total number of threads?
    SizeType masked;    // The number of masked voxels (which are those actually processed)
    SizeType procs;     // The number of voxels processed by this thread
    IndexType* idx;     // Size masked, the masked voxels (which will be actually processed)
    SizeType K1;        // (L+1)(L+2)/2, the number of input SH coefficients
    SizeType K2;        // (2L+1)(2L+2)/2, the number of output SH coefficients
} ThArgs;

void* sh2sqsh_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    if( (nrhs<1) || (nrhs>2) )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh:nrhs","Only 1 or 2 input arguments allowed");
    SizeType ND1  = mxGetNumberOfDimensions( prhs[0] );
    if( ND1<2 )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh:dim","Input 1 must be at least 2-D");
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
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh:dim","Weird size for the last dimension of Input 1");
    //=======================================================================================
    SizeType Kin  = (L+1)*(L+2)/2;
    SizeType Kout = (2*L+1)*(2*L+2)/2;
    //=======================================================================================
    if( nrhs>1 ){
        SizeType ND2  = mxGetNumberOfDimensions( prhs[1] );
        if( ND2==ND1 ){
            if( mxGetDimensions(prhs[1])[ND2-1]!=1 )
                mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh:dim","Input 2 must match the dimensions of Input 1");
            ND2--;
        }
        else if( ND2!=(ND1-1) )
            mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh:dim","Input 2 must match the dimensions of Input 1");
        for( SizeType d=0; d<ND2; ++d ){
            if( mxGetDimensions(prhs[1])[d] != dims[d] ){
                mexPrintf("size(prhs[0],%i)=%i; size(prhs[1],%i)=%i\n",d,dims[d],d,mxGetDimensions(prhs[1])[d]);
                mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh:fov","Inputs 1 and 2 must have the same FOVs");
            }
        }            
    }
    //=======================================================================================
    if( nlhs>1 )
        mexErrMsgIdAndTxt("MyToolbox:sh2squaredsh:nlhs","Only 1 output can be returned");
    if( nlhs==0 )
        return;
    SizeType* dims2 = new SizeType[ND1];
    for( unsigned int d=0; d<ND1-1; ++d )
        dims2[d] = dims[d];
    dims2[ND1-1] = Kout;
    plhs[0] = mxCreateNumericArray( ND1, dims2, mxDOUBLE_CLASS, mxREAL );
    delete[] dims2;
    //=======================================================================================
    SizeType N = shmaths::computeNumberOfSquaredSHFactors( L );
    BufferType factors = new ElementType[N];
    unsigned int*  l1 = new unsigned int[N];
    unsigned int*  l2 = new unsigned int[N];
    unsigned int*  l3 = new unsigned int[N];
    int*           m1 = new int[N];
    int*           m2 = new int[N];
    int*           m3 = new int[N];
    shmaths::computeSquaredSHFactors( L, factors, l1, l2, l3, m1, m2, m3 );
    IndexType* k1 = new IndexType[N];
    IndexType* k2 = new IndexType[N];
    IndexType* k3 = new IndexType[N];
    shmaths::unrollEvenSHIndices( N, l1, m1, k1 );
    shmaths::unrollEvenSHIndices( N, l2, m2, k2 );
    shmaths::unrollEvenSHIndices( N, l3, m3, k3 );
    delete[] l1;
    delete[] l2;
    delete[] l3;
    delete[] m1;
    delete[] m2;
    delete[] m3;
    //=======================================================================================
    BufferType shin  = mxGetDoubles( prhs[0] );
    BufferType mask;
    if(nrhs>1)
        mask = mxGetDoubles( prhs[1] );
    else{
        mask = new ElementType[fov];
        for( IndexType j=0; j<(IndexType)fov; ++j )
            mask[j] = 1.0;
    }
    BufferType shout = mxGetDoubles( plhs[0] );
    SizeType masked  = mataux::nonnullsMxNDArray( mask, fov );
    IndexType* idx   = new IndexType[masked];
    for( SizeType n=0, i=0; n<fov; ++n ){
        if( mataux::notNull(mask[n]) )
            idx[i++] = (IndexType)n;
    }
    //=======================================================================================
#if defined(_HAS_POSIX_THREADS)
    //=======================================================================================
    const unsigned int NTHREADS = sysconf(_SC_NPROCESSORS_CONF);
    pthread_t* threads = new pthread_t[NTHREADS];
    int*       rets    = new int[NTHREADS];
    ThArgs*    args    = new ThArgs[NTHREADS];
    for( unsigned int tid=0; tid<NTHREADS; ++tid ){
        //===================================================================================
        args[tid].N       = N;        //
        args[tid].k1      = k1;       //
        args[tid].k2      = k2;       //
        args[tid].k3      = k3;       //
        args[tid].factors = factors;  //
        args[tid].shin    = shin;     //
        args[tid].shout   = shout;    //
        args[tid].fov     = fov;      //
        args[tid].tid     = tid;      //
        args[tid].nth     = NTHREADS; //
        args[tid].masked  = masked;   //
        args[tid].procs   = 0;        //
        args[tid].idx     = idx;      //
        args[tid].K1      = Kin;      //
        args[tid].K2      = Kout;     //
        //===================================================================================
        rets[tid] = pthread_create(
            &(threads[tid]),
            NULL,
            sh2sqsh_process_fcn, 
            (void*)(&args[tid])    );
        //===================================================================================
    }
    for( unsigned int tid=0; tid<NTHREADS; ++tid ){
        //===================================================================================
        pthread_join( threads[tid], NULL);
        //===================================================================================
    }
    delete[] threads;
    delete[] rets;
    delete[] args;
#else
    ThArgs args;
    args.N       = N;        //
    args.k1      = k1;       //
    args.k2      = k2;       //
    args.k3      = k3;       //
    args.factors = factors;  //
    args.shin    = shin;     //
    args.shout   = shout;    //
    args.fov     = fov;      //
    args.tid     = 0;        //
    args.nth     = 1;        //
    args.masked  = masked;   //
    args.procs   = 0;        //
    args.idx     = idx;      //
    args.K1      = Kin;      //
    args.K2      = Kout;     //
    sh2sqsh_process_fcn( (void*)&args );
#endif
    //=======================================================================================
    delete[] idx;
    delete[] factors;
    delete[] k1;
    delete[] k2;
    delete[] k3;
    if(nrhs<2)
        delete[] mask;
    //=======================================================================================
}

void* sh2sqsh_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    BufferType in  = new ElementType[args->K1];
    BufferType out = new ElementType[args->K2];
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
    for( SizeType i=(args->tid); i<(args->masked); i+=(args->nth) ){
        IndexType voxel = args->idx[i];
        for( SizeType k1=0; k1<args->K1; ++k1 )
            in[k1] = args->shin[k1*(args->fov)+voxel];
        for( SizeType k3=0; k3<args->K2; ++k3 )
            out[k3] = 0.0f;
        for( SizeType n=0; n<args->N; ++n )
            out[args->k3[n]] += in[args->k1[n]]*in[args->k2[n]]*(args->factors[n]);
        for( SizeType k3=0; k3<args->K2; ++k3 )
            args->shout[k3*(args->fov)+voxel] = out[k3];
        (args->procs)++;
    }
#if defined(_HAS_POSIX_THREADS)
    omp_set_num_threads(omp_mthreads);
#endif
    delete[] in;
    delete[] out;
    return (void*)NULL;
}
