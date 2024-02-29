/*==========================================================
 * shodf2samples.c
 *
 * Implements shodf2samples as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "../mathsmex/sphericalHarmonics.h"
#include "../mathsmex/mexToMathsTypes.h"
#include <stdlib.h>
#include <chrono>

#if defined(_HAS_POSIX_THREADS)
#include <pthread.h>
#include <unistd.h>
typedef struct drand48_data drand48_data;
#else
typedef int drand48_data;
#include <cstdlib>
#endif

typedef struct{
    SizeType N;                // The number of sampels to generate
    SizeType Nth;              // The number of samples generated at this thread
    BufferType theta;          // The buffer (size N) where the theta coordinate is stored
    BufferType phi;            // The buffer (size N) where the phi coordinate is stored
    unsigned int L;            // The (even) order of the SH expansion that represents the ODF
    BufferType SH;             // The buffer (size K) with the SH coefficients expanding the ODF
    ElementType c;             // The *magic* constant used by the rejection method
    unsigned short seed16v[3]; // Random generator seeding
    unsigned int tid;          // Which thread is being run?
    unsigned int nth;          // What is the total number of threads?
} ThArgs;

ElementType shodf2samples_compute_c_from_coeffs( const BufferType, const BufferType, const SizeType );
ElementType shodf2samples_compute_c_brute_force( const BufferType, const unsigned int );
void shodf2samples_seed( unsigned short* );
void shodf2samples_randn3d( drand48_data*, double*, double*, double* );
double shodf2samples_evalshodf( const BufferType, const unsigned int, const ElementType,
                                       const ElementType, BufferType, double*, double* );
void* shodf2samples_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    // INPUTS:
    //      prhs[0] (MANDATORY): The SH coefficients, 1xK
    //      prhs[1] (MANDATORY): The number of samples to generate, 1x1
    //      prhs[2] (OPTIONAL): Either c itself (1x1) or the factors to compute it (1xK)
    //      prhs[3] (OPTIONAL): A 1x3 vector to initialize the random number generator
    // ---------------------
    if( (nrhs<2) || (nrhs>4) )
        mexErrMsgIdAndTxt("MyToolbox:shodf2samples:nrhs","Only 2, 3, or 4 input arguments allowed");
    // ---------------------
    if( mxGetNumberOfDimensions(prhs[0]) != 2 )
        mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 1 must be 2-D");
    if( mxGetM(prhs[0]) != 1 )
        mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 1 must be a row vector");
    SizeType K = mxGetN(prhs[0]);
    unsigned int L = (unsigned int)round( (::sqrt(8*(K-1)+9)-3)/2 );
    if( ((SizeType)L+1)*((SizeType)L+2)/2 != K )
        mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 1 has a weird size for an SH expansion");
    BufferType SH = mxGetDoubles( prhs[0] );
    // ---------------------
    if( mxGetNumberOfElements(prhs[1]) != 1 )
        mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 2 must be scalar");
    SizeType N = (SizeType)( mxGetDoubles( prhs[1] )[0] );
    // ---------------------
    ElementType c;
    if( nrhs>2 ){
        if( mxGetNumberOfElements(prhs[2]) == 0 ){
            // Empty array passed
            c = shodf2samples_compute_c_brute_force( SH, L );
        }
        else if( mxGetNumberOfElements(prhs[2]) == 1 ){
            // c itself was passed
            c = mxGetDoubles( prhs[2] )[0];
        }
        else if( mxGetNumberOfElements(prhs[2]) < K ){
            mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","The length of input 3 is shorter that the number of SH coefficients ");
        }
        else{
            if( mxGetNumberOfDimensions(prhs[2]) != 2 )
                mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 3 must be either a scalar or a row vector");
            if( mxGetM(prhs[2]) != 1 )
                mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 3 must be either a scalar or a row vector");
            c = shodf2samples_compute_c_from_coeffs( SH, mxGetDoubles(prhs[2]), K );
        }
    }
    else{
        c = shodf2samples_compute_c_brute_force( SH, L );
    }
    // ---------------------
    unsigned short seed16v[3];
    if( nrhs>3 ){
        if( mxGetNumberOfElements(prhs[3]) != 3 ){
            mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 4 must be a 3-element vector ");
        }
        seed16v[0] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[3])[0]) % 0x1000000000000 );
        seed16v[1] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[3])[1]) % 0x1000000000000 );
        seed16v[2] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[3])[2]) % 0x1000000000000 );
    }
    else{
        shodf2samples_seed( (unsigned short*)seed16v );
    }
    //=======================================================================================
    // OUTPUTS:
    //      plhs[0]: The theta coordiantes (Nx1)
    //      plhs[1]: The phi coordiantes (Nx1)
    //      plhs[2]: The constant c (1x1)
    // ---------------------
    if(nlhs>3)
        mexErrMsgIdAndTxt("MyToolbox:shodf2samples:nlhs","Only up to 3 outputs can be returned");
    if(nlhs==0)
        return;
    SizeType* dims = new SizeType[2];
    dims[1] = 1;
    if(nlhs==1){
        dims[0] = 1;
        plhs[0] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
        delete[] dims;
        mxGetDoubles(plhs[0])[0] = c;
        return;
    }
    dims[0] = N;
    plhs[0] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    plhs[1] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    BufferType theta = mxGetDoubles(plhs[0]);
    BufferType phi   = mxGetDoubles(plhs[1]);
    if(nlhs==3){
        dims[0] = 1;
        plhs[2] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
        mxGetDoubles(plhs[2])[0] = c;
    }
    delete[] dims;
    //=======================================================================================
#if defined(_HAS_POSIX_THREADS)
    //=======================================================================================
    unsigned int NTHREADS = sysconf(_SC_NPROCESSORS_CONF);
    if(NTHREADS>N)
        NTHREADS = N;
    pthread_t* threads = new pthread_t[NTHREADS];
    int*       rets    = new int[NTHREADS];
    ThArgs*    args    = new ThArgs[NTHREADS];
    for( unsigned int tid=0; tid<NTHREADS; ++tid ){
        //===================================================================================
        args[tid].N       = N;        //
        args[tid].Nth     = 0;        //
        args[tid].theta   = theta;    //
        args[tid].phi     = phi;      //
        args[tid].L       = L;        //
        args[tid].SH      = SH;       //
        args[tid].c       = c;        //
        args[tid].seed16v[0] = (unsigned short)(   ( (unsigned long)seed16v[0] * (tid+1) * 7919 ) % 0x1000000000000 );
        args[tid].seed16v[1] = (unsigned short)(   ( (unsigned long)seed16v[1] * (tid+1) * 6983 ) % 0x1000000000000 );
        args[tid].seed16v[2] = (unsigned short)(   ( (unsigned long)seed16v[2] * (tid+1) * 3911 ) % 0x1000000000000 );
        args[tid].tid     = tid;      //
        args[tid].nth     = NTHREADS; //
        //===================================================================================
        rets[tid] = pthread_create(
            &(threads[tid]),
            NULL,
            shodf2samples_process_fcn,
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
    args.Nth     = 0;        //
    args.theta   = theta;    //
    args.phi     = phi;      //
    args.L       = L;        //
    args.SH      = SH;       //
    args.c       = c;        //
    args.seed16v[0] = seed16v[0];
    args.seed16v[1] = seed16v[1];
    args.seed16v[2] = seed16v[2];
    args.tid     = 0;        //
    args.nth     = 1;        //
    shodf2samples_process_fcn( (void*)&args );
#endif
}

void* shodf2samples_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    drand48_data rnd;
#if defined(_HAS_POSIX_THREADS)
    // Re-entrant random number generator
    seed48_r( args->seed16v, &rnd );
#else
    srand( (unsigned int)(args->seed16v[0])*(unsigned int)(args->seed16v[1]) + (unsigned int)(args->seed16v[2]) );
#endif
    double x,y,z,r;
    double theta, phi;
    double sample,target;
    IndexType pos;
    BufferType Ylm  = new ElementType[shmaths::getNumberOfEvenAssociatedLegendrePolynomials(args->L)];
    double* buffer  = shmaths::allocateBufferForEvenAssociatedLegendrePolynomials(args->L);
    double* buffer2 = shmaths::allocateBufferForAssociatedLegendrePolynomials(args->L);
    for( IndexType i=(IndexType)(args->tid); i<(IndexType)(args->N); i+=(args->nth) ){
        SizeType iters = 0;
        SizeType MAXITERS = (SizeType)::ceil(100*(args->c));
        bool rejected = true;
        while( rejected & (iters<MAXITERS) ){
            // ----------------------------------------------------------
            shodf2samples_randn3d( &rnd, &x, &y, &z );
            r = ::sqrt(x*x+y*y+z*z);
            if(r==0)
                continue;
            else{
                theta = ::acos(z/r);
                phi   = ::atan2( y, x );
                if(phi<0)
                    phi += 2*PI;
            }
            // ----------------------------------------------------------
#if defined(_HAS_POSIX_THREADS)
            drand48_r( &rnd, &sample );
#else
            sample = (double)(rand())/RAND_MAX;
#endif
            sample *= (args->c) / (4*PI);
            // ----------------------------------------------------------
            pos    = 0;
            target = shodf2samples_evalshodf( args->SH, args->L,
                                              phi, theta, Ylm, buffer, buffer2 );
            rejected = (sample>target);
            (args->theta)[i] = theta;
            (args->phi)[i] = phi;
            // ----------------------------------------------------------
            ++iters;
        }
        (args->Nth)++;
    }
    delete[] Ylm;
    delete[] buffer;
    delete[] buffer2;
    return (void*)NULL;
}

ElementType shodf2samples_compute_c_from_coeffs( const BufferType SH, const BufferType coeffs, const SizeType K )
{
    ElementType c = 0.0;
    for( IndexType k=0; k<(IndexType)K; ++k )
        c += coeffs[k] * ::abs(SH[k]);
    return c*4*PI;
}

ElementType shodf2samples_compute_c_brute_force( const BufferType SH, const unsigned int L )
{
    BufferType Ylm  = new ElementType[shmaths::getNumberOfEvenAssociatedLegendrePolynomials(L)];
    double* buffer  = shmaths::allocateBufferForEvenAssociatedLegendrePolynomials(L);
    double* buffer2 = shmaths::allocateBufferForAssociatedLegendrePolynomials(L);
    ElementType c   = 0.0;
    unsigned int RES = 100;
    for( unsigned int t=0; t<RES; ++t ){
        ElementType theta = ( (ElementType)t/(ElementType)RES )*PI;
        unsigned int NP = (SizeType)::floor( ::sin(theta) * RES ) + 1;
        for( unsigned int p=0; p<NP; ++p ){
            ElementType phi = ( (ElementType)p/(ElementType)NP )*PI;
            ElementType eval = shodf2samples_evalshodf(SH, L, phi, theta,
                                                       Ylm, buffer, buffer2 );
            if(eval>c)
                c = eval;
        }
    }
    delete[] Ylm;
    delete[] buffer;
    delete[] buffer2;
    return c*4*PI;
}

void shodf2samples_seed( unsigned short* seed16v )
{
    std::chrono::time_point<std::chrono::high_resolution_clock> now;
    unsigned long long us;
    double res;
    // Get current CPU time:
    now = std::chrono::high_resolution_clock::now();
    us = std::chrono::duration_cast<std::chrono::microseconds>
            (now.time_since_epoch()).count();
    seed16v[0] = (unsigned short)( us % 0x1000000000000 );
    // ----
    for( unsigned short k=0, res=1.0f; k<seed16v[0]; ++k ){ res *= PI; }
    now = std::chrono::high_resolution_clock::now();
    us = std::chrono::duration_cast<std::chrono::microseconds>
            (now.time_since_epoch()).count();
    seed16v[1] = (unsigned short)( (us*7907) % 0x1000000000000 );
    // ----
    for( unsigned short k=0, res=1.0f; k<seed16v[1]; ++k ){ res *= PI; }
    now = std::chrono::high_resolution_clock::now();
    us = std::chrono::duration_cast<std::chrono::microseconds>
            (now.time_since_epoch()).count();
    seed16v[2] = (unsigned short)( (us*6991) % 0x1000000000000 );
}

void shodf2samples_randn3d( drand48_data* rnd, double* x, double* y, double* z )
{

    // Use Box-Muller transform
    double u1, u2;
#if defined(_HAS_POSIX_THREADS)
    drand48_r( rnd, &u1 );
    drand48_r( rnd, &u2 );
#else
    u1 = (double)(rand())/RAND_MAX;
    u2 = (double)(rand())/RAND_MAX;
#endif
    *x = ::sqrt( -2.0 * ::log(u1) ) * ::cos(2*PI*u2);
    *y = ::sqrt( -2.0 * ::log(u1) ) * ::sin(2*PI*u2);
#if defined(_HAS_POSIX_THREADS)
    drand48_r( rnd, &u1 );
    drand48_r( rnd, &u2 );
#else
    u1 = (double)(rand())/RAND_MAX;
    u2 = (double)(rand())/RAND_MAX;
#endif
    *z = ::sqrt( -2.0 * ::log(u1) ) * ::cos(2*PI*u2);
    return;
}

double shodf2samples_evalshodf(
    const BufferType SH,
    const unsigned int L,
    const ElementType phi,
    const ElementType theta,
    BufferType Ylm,
    double* buffer,
    double* buffer2
)
{
    double eval = 0.0f;
    shmaths::computeSHMatrixSymmetric( 1, &theta, &phi, L, Ylm, buffer, buffer2 );
    for( unsigned int p=0; p<shmaths::getNumberOfEvenAssociatedLegendrePolynomials(L); ++p )
        eval += ( Ylm[p] * SH[p] );
    return eval;
}
