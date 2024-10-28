/*==========================================================
 * mapl2samples.cpp
 *
 * Implements mapl2samples as a mex function
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "hermitePols.h"
#include "maplMaths.h"
#include "../mathsmex/sanityCheckDTI.h"
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
    SizeType K;                // The number of MAPL coefficients
    SizeType Nth;              // The number of samples generated at this thread
    BufferType x;              // The buffer (size N) where the x coordinate is stored
    BufferType y;              // The buffer (size N) where the y coordinate is stored
    BufferType z;              // The buffer (size N) where the z coordinate is stored
    unsigned int Nmax;         // The maximum degree of the Hermite polynomials
    BufferType coeffs;         // The buffer (size K) with the MAPL coefficients expanding the EAP
    BufferType dti;            // The doffusion tensor representation
    ElementType tau;           // The effective diffusion time
    ElementType c;             // The *magic* constant used by the rejection method
    ElementType sigma;         // The variance of the Gaussian candidates
    unsigned short seed16v[3]; // Random generator seeding
    unsigned int tid;          // Which thread is being run?
    unsigned int nth;          // What is the total number of threads?
};

ElementType mapl2samples_compute_c_brute_force_core( const BufferType, const unsigned int, const ElementType );
ElementType mapl2samples_compute_c_brute_force( const BufferType, const unsigned int, ElementType& );

void mapl2samples_seed( unsigned short* );
void mapl2samples_randn3d( drand48_data*, double*, double*, double* );

THFCNRET mapl2samples_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    // INPUTS:
    //      prhs[0] (MANDATORY): The MAPL coefficients, 1xK
    //      prhs[1] (MANDATORY): The Diffusion Tensor for the voxel
    //      prhs[2] (MANDATORY): The effective diffusion time of the acqusition
    //      prhs[3] (MANDATORY): The number of samples to generate, 1x1
    //      prhs[4] (OPTIONAL): Either c itself (1x1) or the factors to compute it (1xK)
    //      prhs[5] (OPTIONAL): A 1x3 vector to initialize the random number generator
    //      prhs[6] (OPTIONAL): A 1x1 scalar with the number of threads to use
    // ---------------------
    if( (nrhs<4) || (nrhs>7) )
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:nrhs","Only 4 to 7 input arguments allowed");
    // ---------------------
    if( mxGetNumberOfDimensions(prhs[0]) != 2 )
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 1 must be 2-D");
    if( mxGetM(prhs[0]) != 1 )
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 1 must be a row vector");
    // The number of coefficients should be (L+2)(L+4)(2*L+3)/24 for some even L
    SizeType K = mxGetN(prhs[0]);
    unsigned int Nmax = 0;
    while( (Nmax+2)*(Nmax+4)*(2*Nmax+3)/24 < K )
        Nmax += 2;
    if( (Nmax+2)*(Nmax+4)*(2*Nmax+3)/24 > K)
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 1 has a weird size for a MAPL expansion");
    BufferType coeffs = mxGetDoubles( prhs[0] );
    // ---------------------
    if( mxGetNumberOfDimensions(prhs[1]) != 2 )
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 2 must be 2-D");
    if( mxGetNumberOfElements(prhs[1]) != 6 )
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 2 must have exactly 6 elements");
    BufferType dti = mxGetDoubles( prhs[1] );
    // ---------------------
    if( mxGetNumberOfElements(prhs[2]) != 1 )
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 3 must be scalar");
    ElementType tau = mxGetDoubles( prhs[2] )[0];
    // ---------------------
    if( mxGetNumberOfElements(prhs[3]) != 1 )
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 4 must be scalar");
    SizeType N = (SizeType)( mxGetDoubles( prhs[3] )[0] );
    // ---------------------
    ElementType c;
    ElementType sigma;
    if( nrhs>4 ){
        if( mxGetNumberOfElements(prhs[4]) == 0 ){
            // Empty array passed
            c = mapl2samples_compute_c_brute_force( coeffs, Nmax, sigma );
        }
        else if( mxGetNumberOfElements(prhs[4]) == 1 ){
            // c itself was passed
            c = mxGetDoubles( prhs[4] )[0];
            sigma = 1.5;
        }
        else if( mxGetNumberOfElements(prhs[4]) < K ){
            mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","The length of input 5 is shorter that the number of MAPL coefficients");
        }
        else{
            if( mxGetNumberOfDimensions(prhs[4]) != 2 )
                mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 5 must be either a scalar or a row vector");
            if( mxGetM(prhs[4]) != 1 )
                mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 5 must be either a scalar or a row vector");
            c = mapl2samples_compute_c_brute_force( coeffs, Nmax, sigma );
        }
    }
    else{
        c = mapl2samples_compute_c_brute_force( coeffs, Nmax, sigma );
    }
    // ---------------------
    unsigned short seed16v[3];
    if( nrhs>5 ){
        if( mxGetNumberOfElements(prhs[5]) == 0 ){
            // Empty array passed -> auto-seed
            mapl2samples_seed( (unsigned short*)seed16v );
        }
        else if( mxGetNumberOfElements(prhs[5]) == 3 ){
            // Seed actually passed. Use it
            seed16v[0] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[5])[0]) % 0x1000000000000 );
            seed16v[1] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[5])[1]) % 0x1000000000000 );
            seed16v[2] = (unsigned short)( (unsigned long long)(mxGetDoubles(prhs[5])[2]) % 0x1000000000000 );
        }
        else{
            // Weird size
            mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 6 must be a 3-element vector or empty ");
        }
    }
    else{
        mapl2samples_seed( (unsigned short*)seed16v );
    }
    // ---------------------
    unsigned int maxthreads = 1000000;
    if( nrhs>6 ){
        if( mxGetNumberOfElements(prhs[6]) != 1 ){
            mexErrMsgIdAndTxt("MyToolbox:mapl2samples:dim","Input 7 should be a scalar integer");
        }
        maxthreads = (unsigned int)( mxGetDoubles( prhs[6] )[0] );
    }
    //=======================================================================================
    // OUTPUTS:
    //      plhs[0]: The x coordiantes (Nx1)
    //      plhs[1]: The y coordiantes (Nx1)
    //      plhs[2]: The z coordiantes (Nx1)
    //      plhs[3]: The constant c (1x1)
    // ---------------------
    if(nlhs>4)
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:nlhs","Only up to 4 outputs can be returned");
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
    if(nlhs==2){
        delete[] dims;
        mexErrMsgIdAndTxt("MyToolbox:mapl2samples:nlhs","Calls with 2 outputs are ambiguous. Ask for either 0,1,3, or 4 outputs");
    }
    plhs[0] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    plhs[1] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    plhs[2] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
    BufferType x = mxGetDoubles(plhs[0]);
    BufferType y = mxGetDoubles(plhs[1]);
    BufferType z = mxGetDoubles(plhs[2]);
    if(nlhs==4){
        dims[0] = 1;
        plhs[3] = mxCreateNumericArray( 2, dims, mxDOUBLE_CLASS, mxREAL );
        mxGetDoubles(plhs[3])[0] = c;
    }
    delete[] dims;
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
    threader.K      = K;        //
    threader.x      = x;        //
    threader.y      = y;        //
    threader.z      = z;        //
    threader.Nmax   = Nmax;     //
    threader.coeffs = coeffs;   //
    threader.dti    = dti;      //
    threader.tau    = tau;      //
    threader.c      = c;        //
    threader.sigma  = sigma;    //
    threader.seed16v[0] = seed16v[0];
    threader.seed16v[1] = seed16v[1];
    threader.seed16v[2] = seed16v[2];
    //=======================================================================================
    threader.threadedProcess( maxthreads, mapl2samples_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET mapl2samples_process_fcn( void* inargs )
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

    // Allocate the coefficients of Hermite polynomials:
    BufferType coeffsHn;
    hermpols::allocatePolynomialCoeffs( coeffsHn, args->Nmax ); // (Nmax+1) x (Nmax+1)
    // Compute the coefficients:
    hermpols::computePols( coeffsHn, args->Nmax );

    // Allocate memory for basis functions computations:
    BufferType rx, ry, rz, rxHnVals, ryHnVals, rzHnVals, Psi;
    mapl::allocateDictionaryBuffers( rx, ry, rz, rxHnVals, ryHnVals,
                                     rzHnVals, Psi, 1, args->Nmax );

    // Compute eigenvalues and eigenvectors of the tensor model, which
    // will be used later on to un-normalize the space:
    ElementType dti[6];
    ElementType eigval[3];
    ElementType eigvec[9];
    memcpy( dti, args->dti, 6*sizeof(ElementType) );
    dtisc::sanityCheckDTI( (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, 3.0e-3, 'D' );
    ElementType ux0 = sqrt(2.0*eigval[0]*(args->tau));
    ElementType uy0 = sqrt(2.0*eigval[1]*(args->tau));
    ElementType uz0 = sqrt(2.0*eigval[2]*(args->tau));

    // Use normalized scale factors to evaluate at a normalized space:
    ElementType ux = 1.0f;
    ElementType uy = 1.0f;
    ElementType uz = 1.0f;

    ElementType pdf,sample,target;
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
                mapl2samples_randn3d( &rnd, rx, ry, rz );
                // Instead of doing here:
                // rx[0] *= args->sigma
                // ry[0] *= args->sigma
                // rz[0] *= args->sigma
                // and re-normalizing when computing the pdf, we will
                // correct the values later on...
                // ----------------------------------------------------------
                // Evaluate the synthetic 3-D PDF at this sample, and multiply
                // by the *magic* constant to get the uniform sample:
                pdf  = ::exp(-0.5*rx[0]*rx[0]) * ::exp(-0.5*ry[0]*ry[0]) * ::exp(-0.5*rz[0]*rz[0]);
                pdf /= ::sqrt( 8.0*PI*PI*PI );
                // ----------------------------------------------------------
                // Correct the 3-D Gaussian random sample with a variance
                // slightly greater than 1 to ease the fulfillment of the
                // condition c*g(x,y,z)>EAP(x,y,z) in the rejection method:
                rx[0] *= args->sigma;
                ry[0] *= args->sigma;
                rz[0] *= args->sigma;
                pdf /= ( (args->sigma) * (args->sigma) * (args->sigma) );
                // ----------------------------------------------------------
                // Generate a uniform sample in [0,c*g(rx,ry,rz)]
                drand48_r( &rnd, &sample );
                sample *= (pdf*(args->c));
                // ----------------------------------------------------------
                // Evaluate the true EAP described by the MAPL coefficients
                // at the candidate point rx, ry, rz
                hermpols::evaluateAllPols( coeffsHn, rx, rxHnVals, args->Nmax, 1 );
                hermpols::evaluateAllPols( coeffsHn, ry, ryHnVals, args->Nmax, 1 );
                hermpols::evaluateAllPols( coeffsHn, rz, rzHnVals, args->Nmax, 1 );
                mapl::computePsiDictionary( rx, ry, rz,
                                            rxHnVals, ryHnVals, rzHnVals,
                                            Psi,
                                            ux, uy, uz,
                                            1, args->Nmax ); // Psi is 1 x M
                mataux::multiplyMxArrays( Psi, args->coeffs, &target, 1, args->K, 1 ); // 1 x 1
                // ----------------------------------------------------------
                // Check if we must accept or reject:
                rejected = (sample>target);
                if(!rejected){
                    // Unnormalize to get samples in the actual 3-D space
                    accepted[0] = rx[0]*ux0;
                    accepted[1] = ry[0]*uy0;
                    accepted[2] = rz[0]*uz0;
                    // Finally, rotate [rx,ry,rz] according to the matrix
                    // of eigenvectors
                    mataux::multiplyMxArrays( (BufferType)eigvec,
                                              (BufferType)accepted,
                                              (BufferType)acceptedr, 3, 3, 1 );
                    (args->x)[i] = acceptedr[0];
                    (args->y)[i] = acceptedr[1];
                    (args->z)[i] = acceptedr[2];

                }
                // ----------------------------------------------------------
                // Continue
                ++iters;
            }
        }
    }
    while( start < args->getN() );
    
    // Free memory previously allocated
    hermpols::destroyPolynomialCoeffs( coeffsHn );
    mapl::destroyDictionaryBuffers( rx, ry, rz,
                                    rxHnVals, ryHnVals, rzHnVals, Psi );

    return (THFCNRET)NULL;
}

ElementType mapl2samples_compute_c_brute_force( const BufferType coeffs, const unsigned int Nmax, ElementType& sigma0 )
{
    ElementType sigma = 1.1;
    sigma0 = sigma;
    ElementType c0 = INFINITY;
    ElementType c;
    do{
        sigma = sigma0 + 0.05;
        c = mapl2samples_compute_c_brute_force_core( coeffs, Nmax, sigma );
        if(c>=c0)
            break;
        c0 = c;
        sigma0 = sigma;
    }
    while( true );
    return c0;
}

ElementType mapl2samples_compute_c_brute_force_core( const BufferType coeffs, const unsigned int Nmax, const ElementType sigma )
{

    // Allocate the coefficients of Hermite polynomials:
    BufferType coeffsHn;
    hermpols::allocatePolynomialCoeffs( coeffsHn, Nmax ); // (Nmax+1) x (Nmax+1)
    // Compute the coefficients:
    hermpols::computePols( coeffsHn, Nmax );
    // Arrange a grid where the EAP will be computed:
    SizeType    rgrid   = 5;//61;
    SizeType    ngrid   = ( (2*rgrid+1)*(2*rgrid+1)*(2*rgrid+1) + 1 )/2; // Avoid antipodal symmetries
    ElementType extreme = 5.0f;
    ElementType delta   = extreme/rgrid;
    // Allocate memory for basis functions computations:
    BufferType rx, ry, rz, rxHnVals, ryHnVals, rzHnVals, Psi;
    mapl::allocateDictionaryBuffers( rx, ry, rz, rxHnVals, ryHnVals,
                                     rzHnVals, Psi, ngrid, Nmax );
    // Fill the grid
    IndexType pos = 0;
    for( IndexType k=0; k<(2*rgrid+1); ++k ){
        for( IndexType j=0; j<(2*rgrid+1); ++j ){
            for( IndexType i=0; i<(2*rgrid+1); ++i ){
                rx[pos] = i*delta - extreme;
                ry[pos] = j*delta - extreme;
                rz[pos] = k*delta - extreme;
                ++pos;
                if(pos==ngrid){break;}
            }
            if(pos==ngrid){break;}
        }
        if(pos==ngrid){break;}
    }
    // Evaluate Hermite polynomials at grid points
    hermpols::evaluateAllPols( coeffsHn, rx, rxHnVals, Nmax, ngrid );
    hermpols::evaluateAllPols( coeffsHn, ry, ryHnVals, Nmax, ngrid );
    hermpols::evaluateAllPols( coeffsHn, rz, rzHnVals, Nmax, ngrid );
    // Evaluate basis functions at grid points:
    mapl::computePsiDictionary( rx, ry, rz,
                                rxHnVals, ryHnVals, rzHnVals,
                                Psi,
                                1.0f, 1.0f, 1.0f,
                                ngrid, Nmax ); // Psi is ngrid x M
    // Create a vector to store the evaluations of the EAP at grid points:
    BufferType evals = new ElementType[ngrid];
    // Multiply the dictionary by the expansion coefficients to get
    // the EAP at grid points:
    SizeType M = mapl::numBasisFunctions(Nmax);
    mataux::multiplyMxArrays( Psi, coeffs, evals, ngrid, M, 1 ); // ngrid x 1
    // Note the last position of evals, evals[ngrid-1], corresponds
    // to the origin {rx=0,ry=0,rz=0}, where the EAP is assumed to
    // show its maximum value. We will use this value to establish a
    // threshold for the tails of the EAP:
    ElementType threshold = 1.0e-4 * evals[ngrid-1];
    // Free memory previously allocated and no longer necessary:
    hermpols::destroyPolynomialCoeffs( coeffsHn );
    mapl::destroyDictionaryBuffers( rx, ry, rz,
                                    rxHnVals, ryHnVals, rzHnVals, Psi );
    // Divide these EAP values by the values of a 3-D N(0,1) PDF defined at
    // these same grid points, and keep the maximum value
    ElementType maxval = -1.0e6;
    pos = 0;
    for( IndexType k=0; k<(2*rgrid+1); ++k ){
        for( IndexType j=0; j<(2*rgrid+1); ++j ){
            for( IndexType i=0; i<(2*rgrid+1); ++i ){
                ElementType c = 1.0/::sqrt(8*PI*PI*PI)/(sigma*sigma*sigma);
                c *= ::exp( -0.5*(i*delta - extreme)*(i*delta - extreme) / (sigma*sigma) );
                c *= ::exp( -0.5*(j*delta - extreme)*(j*delta - extreme) / (sigma*sigma) );
                c *= ::exp( -0.5*(k*delta - extreme)*(k*delta - extreme) / (sigma*sigma) );
                c  = (evals[pos]+threshold)/(c+threshold);
                if( c>maxval )
                    maxval = c;
                ++pos;
                if(pos==ngrid){break;}
            }
            if(pos==ngrid){break;}
        }
        if(pos==ngrid){break;}
    }
    // Delete allocated memory
    delete[] evals;

    return maxval;
}

void mapl2samples_seed( unsigned short* seed16v )
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

void mapl2samples_randn3d( drand48_data* rnd, double* x, double* y, double* z )
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
