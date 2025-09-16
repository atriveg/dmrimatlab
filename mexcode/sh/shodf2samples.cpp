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
//#include "matrix.h"
#include "math.h"
#include "../mathsmex/sphericalHarmonics.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

#include <stdlib.h>
#include <chrono>


#if ( defined (__APPLE__) || defined (__MACH__) || defined (_WIN32) )
#include <limits.h>
#include <stdint.h>
#include <sys/types.h>
/**
 * The re-entrant routine drand48_r seems to be unavailable for
 * both MAC and Windows environments. We will mimick its implementation
 * (taken from GNU Gnulib) to be able to generate "independent" random
 * numbers from concurrent threads
*/
typedef struct _drand48_data_ {
    unsigned short int __x[3] = {0,0,0};
    unsigned short int __c = 0xb;     // 11
    unsigned long long int __a = 0x5deece66d; // 25214903917
} drand48_data;
void seed48_r( unsigned short int seed16v[3], drand48_data* );
int drand48_r( drand48_data*, double* );
#else
typedef struct drand48_data drand48_data;
#endif

class ThArgs : public DMRIThreader
{
public:
    SizeType N;                // The number of samples to generate
    SizeType Nth;              // The number of samples generated at this thread
    BufferType theta;          // The buffer (size N) where the theta coordinate is stored
    BufferType phi;            // The buffer (size N) where the phi coordinate is stored
    unsigned int L;            // The (even) order of the SH expansion that represents the ODF
    BufferType SH;             // The buffer (size K) with the SH coefficients expanding the ODF
    ElementType c;             // The *magic* constant used by the rejection method
    unsigned short seed16v[3]; // Random generator seeding
    unsigned int tid;          // Which thread is being run?
    unsigned int nth;          // What is the total number of threads?
};

ElementType shodf2samples_compute_c_from_coeffs( const BufferType, const BufferType, const SizeType );
ElementType shodf2samples_compute_c_brute_force( const BufferType, const unsigned int );
void shodf2samples_seed( unsigned short* );
void shodf2samples_randn3d( drand48_data*, double*, double*, double* );
double shodf2samples_evalshodf( const BufferType, const unsigned int, const ElementType,
                                       const ElementType, BufferType, double*, double* );
THFCNRET shodf2samples_process_fcn( void* );

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
    //      prhs[4] (OPTIONAL): A 1x1 scalar with the number of threads to use
    // ---------------------
    if( (nrhs<2) || (nrhs>5) )
        mexErrMsgIdAndTxt("MyToolbox:shodf2samples:nrhs","Only 2 to 5 input arguments allowed");
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
            mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","The length of input 3 is shorter that the number of SH coefficients");
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
        if( mxGetNumberOfElements(prhs[3]) == 0 ){
            // Empty array passed -> auto-seed
            shodf2samples_seed( (unsigned short*)seed16v );
        }
        else if( mxGetNumberOfElements(prhs[3]) == 3 ){
            // Seed actually passed. Use it
            seed16v[0] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[3])[0]) % 0x1000000000000 );
            seed16v[1] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[3])[1]) % 0x1000000000000 );
            seed16v[2] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[3])[2]) % 0x1000000000000 );
        }
        else{
            // Weird size
            mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 4 must be a 3-element vector or empty ");
        }
    }
    else{
        shodf2samples_seed( (unsigned short*)seed16v );
    }
    // ---------------------
    unsigned int maxthreads = 1000000;
    if( nrhs>4 ){
        if( mxGetNumberOfElements(prhs[4]) != 1 ){
            mexErrMsgIdAndTxt("MyToolbox:shodf2samples:dim","Input 4 should be a scalar integer");
        }
        maxthreads = (unsigned int)( mxGetDoubles( prhs[4] )[0] );
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
    maxthreads = get_number_of_threads( maxthreads );
    if(maxthreads>N)
        maxthreads = N;
    //=======================================================================================
    unsigned int chunksz = (N/maxthreads)/20;
    if(chunksz<1)
        chunksz = 1;
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, chunksz );
    // Own values:
    threader.theta   = theta;    //
    threader.phi     = phi;      //
    threader.L       = L;        //
    threader.SH      = SH;       //
    threader.c       = c;        //
    threader.seed16v[0] = seed16v[0];
    threader.seed16v[1] = seed16v[1];
    threader.seed16v[2] = seed16v[2];
    //=======================================================================================
    threader.threadedProcess( maxthreads, shodf2samples_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET shodf2samples_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;
    
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads_thread(1);

    // Get a unique thread ordinal to create different seeds for each thread
    unsigned short seed16v[3];
    unsigned int thid = args->getThid();
    seed16v[0] = (unsigned short)(   ( (unsigned long)(args->seed16v[0]) * (thid+1) * 7919 ) % 0x1000000000000 );
    seed16v[1] = (unsigned short)(   ( (unsigned long)(args->seed16v[1]) * (thid+1) * 6983 ) % 0x1000000000000 );
    seed16v[2] = (unsigned short)(   ( (unsigned long)(args->seed16v[2]) * (thid+1) * 3911 ) % 0x1000000000000 );
            
    // Init the random number generator
    drand48_data rnd;
    // Re-entrant random number generator
    seed48_r( seed16v, &rnd );
    
    double x,y,z,r;
    double theta, phi;
    double sample,target;
    
    BufferType Ylm  = new ElementType[shmaths::getNumberOfEvenAssociatedLegendrePolynomials(args->L)];
    double* buffer  = shmaths::allocateBufferForEvenAssociatedLegendrePolynomials(args->L);
    double* buffer2 = shmaths::allocateBufferForAssociatedLegendrePolynomials(args->L);
    
    SizeType MAXITERS = (SizeType)::ceil(100*(args->c));
    
    // ---------------------------------------------------------------
    // Loop through the sample positions
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( IndexType i=start; i<end; ++i ){
            SizeType iters = 0;
            bool rejected = true;
            while( rejected && (iters<MAXITERS) ){
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
                drand48_r( &rnd, &sample );
                sample *= (args->c) / (4*PI);
                // ----------------------------------------------------------
                target = shodf2samples_evalshodf( args->SH, args->L,
                        phi, theta, Ylm, buffer, buffer2 );
                rejected = (sample>target);
                (args->theta)[i] = theta;
                (args->phi)[i] = phi;
                // ----------------------------------------------------------
                ++iters;
            }
        }
    }
    while( start < args->getN() );
    
    blas_num_threads_thread(blas_threads);

    delete[] Ylm;
    delete[] buffer;
    delete[] buffer2;
    return (THFCNRET)NULL;
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

    drand48_r( rnd, &u1 );
    drand48_r( rnd, &u2 );
    *x = ::sqrt( -2.0 * ::log(u1) ) * ::cos(2*PI*u2);
    *y = ::sqrt( -2.0 * ::log(u1) ) * ::sin(2*PI*u2);

    drand48_r( rnd, &u1 );
    drand48_r( rnd, &u2 );
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

#if ( defined (__APPLE__) || defined (__MACH__) || defined (_WIN32) )
/**
 * The re-entrant routine drand48_r seems to be unavailable for
 * both MAC and Windows environments. We will mimick its implementation
 * (taken from GNU Gnulib) to be able to generate "independent" random
 * numbers from concurrent threads
*/
void seed48_r( unsigned short int seed16v[3], drand48_data* rnd )
{
    rnd->__x[0] = seed16v[0];
    rnd->__x[1] = seed16v[1];
    rnd->__x[2] = seed16v[2];
    rnd->__c = 0xb; // 11
    rnd->__a = 0x5deece66d; // 25214903917
    return;
}

/**
* Implementation directly taken from the __erand48_r
*/
int drand48_r( drand48_data* rnd, double* x )
{
    /* Compute next state.  */
    uint64_t X;
    uint64_t result;
    X = (uint64_t) rnd->__x[2] << 32 | (uint32_t) rnd->__x[1] << 16 | rnd->__x[0];
    result = X * rnd->__a + rnd->__c;
    rnd->__x[0] = result & 0xffff;
    rnd->__x[1] = (result >> 16) & 0xffff;
    rnd->__x[2] = (result >> 32) & 0xffff;
    /* The GNU implementation is based on IEEE-754 representations
     * of floating point numbers. Since ieee754.h is not present in
     * Windows or MAC, use floating point arithmetic instead:
     */
    *x = (double)( result & 0xffffffffffff ) / (double)(0x1000000000000);
    return 0;
}

#endif
