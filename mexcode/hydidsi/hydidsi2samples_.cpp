/*==========================================================
 * hydidsi2samples_.cpp
 *
 * Implements hydidsi2samples as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"
#include "../mathsmex/DSIInterp.h"
#include <stdlib.h>
#include <chrono>
#include <cmath>

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
    SizeType Nx;               // Lattice size 0
    SizeType Ny;               // Lattice size 1
    SizeType Nz;               // Lattice size 2
    SizeType Nth;              // The number of samples generated at this thread
    BufferType x;              // The buffer (size N) where the x coordinate is stored
    BufferType y;              // The buffer (size N) where the y coordinate is stored
    BufferType z;              // The buffer (size N) where the z coordinate is stored
    BufferType grid;           // The expanded grid (size (2*Nx+5)x(2*Ny+5)x(2*Nz+5)) where the EAP is rendered
    ElementType c;             // The *magic* constant used by the rejection method
    ElementType sigma;         // The variance of the Gaussian candidates
    char interp;               // The interpolation type
    unsigned short seed16v[3]; // Random generator seeding
    unsigned int tid;          // Which thread is being run?
    unsigned int nth;          // What is the total number of threads?
};

#ifndef PI
#define PI 3.14159265358979323846
#endif

THFCNRET hydidsi2samples_process_fcn( void* );

void hydidsi2samples_randn3d( drand48_data*, double&, double&, double& );

void hydidsi2samples_thomas3D( const BufferType, BufferType*, BufferType*, BufferType, const unsigned int*, const unsigned int = 2 );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    // INPUTS:
    //      prhs[0]: The HYDI-DSI coefficients, 1 x K, with K = ((2*Nx+1)(2*Ny+1)(2*Nz+1)+1)/2
    //      prhs[1]: The lattice size, 3x1
    //      prhs[2]: The constant c, 1x1
    //      prhs[3]: The normalized variance sigma, 1x1
    //      prhs[4]: The number of samples to generate, 1x1
    //      prhs[5]: 'n', 'l' or 's', the interpolation type,
    //      prhs[6]: A 1x3 vector to initialize the random number generator
    //      prhs[7]: A 1x1 scalar with the number of threads to use
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2samples_:callstack","This function should only be called from hydidsi2samples");
    else if( strcmp(callerFunc,"hydidsi2samples") )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2samples_:callstack","This function should only be called from hydidsi2samples");
    //=======================================================================================
    if( nrhs!=8 )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2samples_:nrhs","Exactly 8 arguments must be passed");
    //=======================================================================================
    SizeType Nx = (SizeType)(mxGetDoubles(prhs[1])[0]);
    SizeType Ny = (SizeType)(mxGetDoubles(prhs[1])[1]);
    SizeType Nz = (SizeType)(mxGetDoubles(prhs[1])[2]);
    SizeType K = ( (2*Nx+1)*(2*Ny+1)*(2*Nz+1) + 1 )/2; // Number of EAP coefficients
    //=======================================================================================
    ElementType c = mxGetScalar( prhs[2] );
    ElementType sigma = mxGetScalar( prhs[3] );
    SizeType N = (SizeType)(mxGetScalar(prhs[4]));
    char interp = (char)( mxGetChars(prhs[5])[0] );
    unsigned short seed16v[3];
    seed16v[0] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[6])[0]) % 0x1000000000000 );
    seed16v[1] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[6])[1]) % 0x1000000000000 );
    seed16v[2] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[6])[2]) % 0x1000000000000 );
    unsigned int maxthreads = (unsigned int)( mxGetDoubles(prhs[7])[0] );
    //=======================================================================================
    // OUTPUTS:
    //      plhs[0]: The x coordiantes (Nx1)
    //      plhs[1]: The y coordiantes (Nx1)
    //      plhs[2]: The z coordiantes (Nx1)
    if(nlhs!=3)
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2samples_:nlhs","Exactly 3 outputs will be generated");
    plhs[0] = mxCreateDoubleMatrix( N, 1, mxREAL );
    plhs[1] = mxCreateDoubleMatrix( N, 1, mxREAL );
    plhs[2] = mxCreateDoubleMatrix( N, 1, mxREAL );
    //=======================================================================================
    maxthreads = get_number_of_threads( maxthreads );
    if(maxthreads>N)
        maxthreads = N;
    //=======================================================================================
    unsigned int chunksz = (N/maxthreads)/10;
    if(chunksz<1)
        chunksz = 1;
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, chunksz );
    // Own values:
    threader.Nx     = Nx;
    threader.Ny     = Ny;
    threader.Nz     = Nz;
    threader.x      = mxGetDoubles(plhs[0]);
    threader.y      = mxGetDoubles(plhs[1]);
    threader.z      = mxGetDoubles(plhs[2]);
    threader.c      = c;
    threader.sigma  = sigma;
    threader.interp = interp;
    threader.seed16v[0] = seed16v[0];
    threader.seed16v[1] = seed16v[1];
    threader.seed16v[2] = seed16v[2];
    // From the EAP in the "collapsed" form, compute the EAP
    // at the extended grid points:
    BufferType values;
    SizeType M = dsiinterp::allocate_extended_grid( values, Nx, Ny, Nz );
    mataux::setValueMxArray(values,M,1,0.0);
    dsiinterp::fill_extended_grid( mxGetDoubles(prhs[0]), values, Nx, Ny, Nz );
    // The threader object will store a buffer the same size:
    dsiinterp::allocate_extended_grid( threader.grid, Nx, Ny, Nz );
    mataux::setValueMxArray(threader.grid,M,1,0.0);
    if(interp=='s'){
        // If spline interpolation is used, this buffer is computed by
        // solving a block-tridiagonal system with Thomas' algorithm
        SizeType sizes[3];
        sizes[0] = 2*Ny+5;  // The 0 -> y swap is not a bug, but consistent with meshgrid
        sizes[1] = 2*Nx+5;  // The 1 -> x swap is not a bug, but consistent with meshgrid
        sizes[2] = 2*Nz+5;
        BufferType* gamma;
        BufferType* theta;
        dsiinterp::allocate_thomas3D( gamma, theta, (SizeBuffer)sizes );
        dsiinterp::thomas3D( values, gamma, theta, threader.grid, (SizeBuffer)sizes );
        delete[] gamma; // No longer necessary
        delete[] theta; // No longer necessary
    }
    else{
        // Otherwise, this buffer is copied straight to the threader
        // object
        memcpy( threader.grid, values, M*sizeof(ElementType) );
    }
    delete[] values; // No longer necessary
    //=======================================================================================
    threader.threadedProcess( maxthreads, hydidsi2samples_process_fcn );
    //=======================================================================================
    delete[] threader.grid;
    return;
}

THFCNRET hydidsi2samples_process_fcn( void* inargs )
{
    ThArgs* args = (ThArgs*)inargs;

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

    double rx, ry, rz;
    ElementType pdf, sample, target;
    ElementType accepted[3];
    ElementType acceptedr[3];

    SizeType MAXITERS = (SizeType)::ceil(1000*(args->c));
    
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
            (args->x)[i] = NAN;
            (args->y)[i] = NAN;
            (args->z)[i] = NAN;
            bool rejected = true;
            while( rejected && (iters<MAXITERS) ){
                // ----------------------------------------------------------
                // In principle, we should draw a 3-D sample distributed as
                // a Gaussian with zero mean and identity covariance matrix
                // (which is the 0-th order expansion of MAPL). However, if
                // we geneate a Gaussian with a variance slighly greate than
                // 1 we achieve heavier tales in its PDF, so that it
                // becomes easier to fulfill the requirement c*g(x,y,z) >
                // EAP(x,y,z) necessary for the rejection method, vastly
                // improving the value of c:
                hydidsi2samples_randn3d( &rnd, rx, ry, rz );
                // Instead of doing here:
                // rx *= args->sigma
                // ry *= args->sigma
                // rz *= args->sigma
                // and re-normalizing when computing the pdf, we will
                // correct the values later on...
                // ----------------------------------------------------------
                // Evaluate the synthetic 3-D PDF at this sample, and multiply
                // by the *magic* constant to get the uniform sample:
                pdf  = exp(-0.5*rx*rx) * exp(-0.5*ry*ry) * exp(-0.5*rz*rz);
                pdf /= sqrt( 8.0*PI*PI*PI );
                // ----------------------------------------------------------
                // Correct the 3-D Gaussian random sample with a variance
                // slightly greater than 1 to ease the fulfillment of the
                // condition c*g(x,y,z)>EAP(x,y,z) in the rejection method.
                // Correct also with the number of nodes along each dimension
                rx *= (args->sigma) * (args->Nx);
                ry *= (args->sigma) * (args->Ny);
                rz *= (args->sigma) * (args->Nz);
                pdf /= ( (args->sigma) * (args->sigma) * (args->sigma) );
                pdf /= ( (args->Nx) * (args->Ny) * (args->Nz) );
                // ----------------------------------------------------------
                // Generate a uniform sample in [0,c*g(rx,ry,rz)]
                drand48_r( &rnd, &sample );
                sample *= (pdf*(args->c));
                // ----------------------------------------------------------
                // Evaluate the true EAP by cubic-spline interpolation in
                // the lattice:
                switch(args->interp){
                    case 'n':
                        dsiinterp::interpolate_NN( args->grid,
                                  args->Nx, args->Ny, args->Nz,
                                  rx, ry, rz, target );
                        break;
                    case 'l':
                        dsiinterp::interpolate_LN( args->grid,
                                  args->Nx, args->Ny, args->Nz,
                                  rx, ry, rz, target );
                        break;
                    case 's':
                        dsiinterp::interpolate_BS( args->grid,
                                  args->Nx, args->Ny, args->Nz,
                                  rx, ry, rz, target );
                        break;
                }
                // ----------------------------------------------------------
                // Check if we must accept or reject:
                rejected = (sample>target);
                if(!rejected){
                    (args->x)[i] = rx;
                    (args->y)[i] = ry;
                    (args->z)[i] = rz;
                }
                // ----------------------------------------------------------
                // Continue
                ++iters;
            }
        }
    }
    while( start < args->getN() );

    return (THFCNRET)NULL;
}

void hydidsi2samples_randn3d( drand48_data* rnd, double& x, double& y, double& z )
{

    // Use Box-Muller transform
    double u1, u2;

    drand48_r( rnd, &u1 );
    drand48_r( rnd, &u2 );
    x = sqrt( -2.0 * log(u1) ) * cos(2*PI*u2);
    y = sqrt( -2.0 * log(u1) ) * sin(2*PI*u2);

    drand48_r( rnd, &u1 );
    drand48_r( rnd, &u2 );
    z = sqrt( -2.0 * log(u1) ) * cos(2*PI*u2);

    return;
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
