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
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"
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
    BufferType coeffs;         // The buffer (size K) with the HYDI-DSI coefficients expanding the EAP
    ElementType c;             // The *magic* constant used by the rejection method
    ElementType sigma;         // The variance of the Gaussian candidates
    unsigned short seed16v[3]; // Random generator seeding
    unsigned int tid;          // Which thread is being run?
    unsigned int nth;          // What is the total number of threads?
};

#ifndef PI
#define PI 3.14159265358979323846
#endif

THFCNRET hydidsi2samples_process_fcn( void* );

void hydidsi2samples_randn3d( drand48_data*, double&, double&, double& );

bool hydidsi2samples_interpolate( const BufferType,
                                  const SizeType&, const SizeType&, const SizeType&,
                                  const double&, const double&, const double&, double& );

void hydidsi2samples_splines( const double&, double&, double&, double&, double& );

void hydidsi2samples_thomas3D( const BufferType, BufferType, const unsigned int*, const unsigned int = 2 );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    // INPUTS:
    //      prhs[0]: The HYDI-DSI coefficients, (2*Nx+5)x(2*Ny+5)x(2*Nz+5)
    //      prhs[1]: The constant c, 1x1
    //      prhs[2]: The normalized variance sigma, 1x1
    //      prhs[3]: The number of samples to generate, 1x1
    //      prhs[4]: A 1x3 vector to initialize the random number generator
    //      prhs[5]: A 1x1 scalar with the number of threads to use
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
    if( nrhs!=6 )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2samples_:nrhs","Exactly 6 arguments must be passed");
    //=======================================================================================
    if( mxGetNumberOfDimensions(prhs[0]) != 3 )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2samples_:dim","Input 1 must be 3-D");
    SizeType M  = mxGetNumberOfElements(prhs[0]);
    SizeType Nx = mxGetDimensions(prhs[0])[0];
    SizeType Ny = mxGetDimensions(prhs[0])[1];
    SizeType Nz = mxGetDimensions(prhs[0])[2];
    unsigned int sizes[3];
    sizes[0] = Nx;
    sizes[1] = Ny;
    sizes[2] = Nz;
    Nx = (Nx-5)/2;
    Ny = (Ny-5)/2;
    Nz = (Nz-5)/2;
    if( (Nx<=0) || (Ny<=0) || (Nz<=0) )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2samples_:nrhs","Wrong sizes at input 0");
    //=======================================================================================
    ElementType c = mxGetScalar( prhs[1] );
    ElementType sigma = mxGetScalar( prhs[2] );
    SizeType N = (SizeType)(mxGetScalar(prhs[3]));
    unsigned short seed16v[3];
    seed16v[0] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[4])[0]) % 0x1000000000000 );
    seed16v[1] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[4])[1]) % 0x1000000000000 );
    seed16v[2] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[4])[2]) % 0x1000000000000 );
    unsigned int maxthreads = (unsigned int)( mxGetDoubles(prhs[5])[0] );
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
    threader.seed16v[0] = seed16v[0];
    threader.seed16v[1] = seed16v[1];
    threader.seed16v[2] = seed16v[2];
    // From the values of the EAP at grid points, use Thomas'
    // algorithm recursively to obtain the grid coefficients
    // of the B-spline
    threader.coeffs = new ElementType[M];
    hydidsi2samples_thomas3D( mxGetDoubles(prhs[0]), threader.coeffs, (unsigned int*)sizes );
    //=======================================================================================
    threader.threadedProcess( maxthreads, hydidsi2samples_process_fcn );
    //=======================================================================================
    delete[] threader.coeffs;
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
                // a Gaussian with zero mean and Identity covariance matrix
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
                hydidsi2samples_interpolate( args->coeffs,
                                  args->Nx, args->Ny, args->Nz,
                                  rx, ry, rz, target );
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

bool hydidsi2samples_interpolate( const BufferType grid,
                                  const SizeType& Nx, const SizeType& Ny, const SizeType& Nz,
                                  const double& x, const double& y, const double& z, double& result )
{
    // Perform B-spline interpolation within a grid with size (2Nx+5) x (2Ny+5)  x (2Nz+5)
    // The EAP is estimated in a mesh (-Nx..Nx) x (-Ny..Ny) x (-Nz..Nz). We extend each dimension
    // with a 0 value at each end: (-Nx-1..Nx+1) x (-Ny-1..Ny+1) x (-Nz-1..Nz+1). Then, we
    // add yet and additional node at each end because B-splines require it.
    // Beyond the additional 0-valued node, the EAP is assumed to be 0:
    result = 0.0;
    if( (x<=-(double)Nx-1) || (x>=(double)Nx+1) || (y<=-(double)Ny-1) || (y>=(double)Ny+1) || (z<=-(double)Nz-1) || (z>=(double)Nz+1) )
        return false;
    // Otherwise, we are inbounds and we can actually interpolate. Compute
    // the associated grid coordinates for each dimension:
    unsigned int xi = (unsigned int)( floor(x) + (Nx+2.0) ); // 1 < xi < 2*Nx+3
    unsigned int yi = (unsigned int)( floor(y) + (Ny+2.0) ); // 1 < yi < 2*Ny+3
    unsigned int zi = (unsigned int)( floor(z) + (Nz+2.0) ); // 1 < zi < 2*Nz+3
    double xg = x - floor(x); // 0 <= xg <= 1
    double yg = y - floor(y); // 0 <= yg <= 1
    double zg = z - floor(z); // 0 <= zg <= 1

    // Evaluate all the kernel functions for B-spline interpolation:
    double Bsplines[3][4];
    hydidsi2samples_splines( xg, Bsplines[0][0], Bsplines[0][1], Bsplines[0][2], Bsplines[0][3] );
    hydidsi2samples_splines( yg, Bsplines[1][0], Bsplines[1][1], Bsplines[1][2], Bsplines[1][3] );
    hydidsi2samples_splines( zg, Bsplines[2][0], Bsplines[2][1], Bsplines[2][2], Bsplines[2][3] );
    // B-spline interpolations will use the 4 x 4 x 4 nearest neighbors of the grid point
    // with coordinates [xi,yi,zi]. Let's compute:
    for( unsigned short i=0; i<4; ++i ){
        for( unsigned short j=0; j<4; ++j ){
            for( unsigned short k=0; k<4; ++k ){
                // The value at each grid point is multiplied by the
                // tensor product of the 4 kernels:
                // result += (grid[xi+i-1][yi+j-1][zi+k-1])
                //     * Bsplines[0][i] * Bsplines[1][j] * Bsplines[2][k];
                // BUT: the grid is stored in a linear format, so we have to
                // translate the triple indexing [xi,yi,zi] to a unique linear
                // index
                IndexType lidx = ((int)(xi+i)-1) + ((int)(yi+j)-1)*(2*Nx+5) + ((int)(zi+k)-1)*(2*Nx+5)*(2*Ny+5);
                result += (grid[lidx])
                     * Bsplines[0][i] * Bsplines[1][j] * Bsplines[2][k];
            }
        }
    }
    return true;
}

void hydidsi2samples_splines( const double& x, double& b0, double& b1, double& b2, double& b3 )
{
    // Evaluate the 4 cubic B-spline functions using Horner's algotihm
    // Note x is assumed (but not checked) to be in [0,1]
    // Cubic terms:
    b0 = -1.0;
    b1 = 3.0;
    b2 = -3.0;
    b3 = 1.0;
    // Cuadratic terms:
    b0 = b0*x + 3.0;
    b1 = b1*x - 6.0;
    b2 = b2*x + 3.0;
    b3 = b3*x; // Just 0!
    // Linear terms:
    b0 = b0*x - 3.0;
    b1 = b1*x; // Just 0!
    b2 = b2*x + 3.0;
    b3 = b3*x; // Just 0!
    // Independent terms:
    b0 = b0*x + 1.0;
    b1 = b1*x + 4.0;
    b2 = b2*x + 1.0;
    b3 = b3*x; // Just 0!
    // Factor 6.0:
    b0 /= 6;
    b1 /= 6;
    b2 /= 6;
    b3 /= 6;
    return;
}

/** Solve a block-tridiagonal system by means of the Thomas algorithm */
void hydidsi2samples_thomas3D( const BufferType values, BufferType results, const unsigned int* sizes, const unsigned int order )
{
    // ---------------------------
    unsigned int M = 1;
    for( unsigned int k=0; k<order; ++k ){ M *= sizes[k]; }
    // ---------------------------
    BufferType  gamma = new ElementType[ sizes[order]+1 ];
    BufferType* theta = new BufferType[  sizes[order]+1 ];
    for( unsigned int k=0; k<sizes[order]+1; ++k )
        theta[k] = new ElementType[M];
    // ---------------------------
    // Forwards swept to compute the recursive multiplicative factors:
    gamma[0] = 1;
    for( unsigned int d=0; d<M; ++d )
        theta[0][d] = 0.0;
    for( unsigned int k=0; k<sizes[order]; ++k ){
        gamma[k+1] = -1.0/( gamma[k] + 4.0 );
        if( order==0 )
            theta[k+1][0] = (-gamma[k+1])*( 6.0*values[k] - theta[k][0] );
        else{
            hydidsi2samples_thomas3D( &values[k*M], theta[k+1], (unsigned int*)sizes, order-1 );
            for( unsigned int d=0; d<M; ++d )
                theta[k+1][d] = (-gamma[k+1])*( 6.0*theta[k+1][d] - theta[k][d] );
        }
    }
    // ---------------------------
    // Backwards swept
    for( unsigned int d=0; d<M; ++d )
        results[(sizes[order]-1)*M+d] = (  (theta[sizes[order]][d])/(1.0-gamma[sizes[order]])  )*gamma[sizes[order]] + theta[sizes[order]][d];
    for( int k=sizes[order]-2; k>=0; --k ){
        for( unsigned int d=0; d<M; ++d )
            results[k*M+d] = (gamma[k+1])*(results[(k+1)*M+d]) + (theta[k+1][d]);
    }
    // ---------------------------
    // Memory cleanup
    for( unsigned int k=0; k<sizes[order]+1; ++k )
        delete[] theta[k];
    delete[] gamma;
    delete[] theta;
    // ---------------------------
    return;
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
