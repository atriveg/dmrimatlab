/*==========================================================
 * coeffs2localdispersion.cpp
 *
 * This is a core function to several functions in the
 * toolbox, and should only be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2026 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"
#include "../pixelIterators/iterators.h"

#include <iostream>


typedef struct LDispesionIOData{
    // Input:
    BufferType coeffs;
    BufferType mask;
    // Output:
    BufferType disp;
} LDispesionIOData;

typedef struct LDispesionParameters{
    unsigned int M;   // The number of expansion coefficients.
    BufferType wm;    // Weights assigned to compute dot products
    SizeBuffer fov;   // Field of view of the images
    SizeBuffer nhood; // Size of the neighborhood to search

} LDispesionParameters;

class ThArgs : public DMRIThreader
{
public:
    LDispesionParameters* params; // Compute (or not) the eigenvectors
    LDispesionIOData*     io;
};

THFCNRET coeffs2localdispersion_process_fcn( void* );

bool coeffs2localdispersion_callstack( const char* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     *
     * prhs[0]: coeffs, the signal coefficients, M x X x Y x Z
     * prhs[1]: mask, the mask to be used, X x Y x Z
     * prhs[2]: weights, M x 1, the positive weights to define the dot product
     * prhs[3]: neigh, 3x1, the radii of the neighborhood used for computations
     * prhs[4]: maxthreads, 1x1, the maximum number of threads to use
     *
     *  OUTPUTS:
     *
     * plhs[0]: disp, X x Y x Z, the local dispersion index
     * *
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( !coeffs2localdispersion_callstack(callerFunc) )
        mexErrMsgIdAndTxt("MyToolbox:coeffs2localdispersion:callstack","Function called from wrong m-file");

    //=======================================================================================
    if(nrhs!=5)
        mexErrMsgIdAndTxt("MyToolbox:coeffs2localdispersion:nrhs","Exactly 5 input arguments are required");
    //=======================================================================================
    LDispesionIOData     io;
    LDispesionParameters params;
    //=======================================================================================
    SizeType ndim = mxGetNumberOfDimensions(prhs[0]);
    SizeType dims[4] = {1,1,1,1};
    for( unsigned int d=0; d<ndim; ++d )
        dims[d] = mxGetDimensions(prhs[0])[d];
    SizeType M = dims[0];
    SizeType fov[3];
    fov[0] = dims[1];
    fov[1] = dims[2];
    fov[2] = dims[3];
    SizeType nhood[3];
    nhood[0] = (SizeType)( mxGetDoubles(prhs[3])[0] );
    nhood[1] = (SizeType)( mxGetDoubles(prhs[3])[1] );
    nhood[2] = (SizeType)( mxGetDoubles(prhs[3])[2] );
    //=======================================================================================
    // Process inputs:
    // ------
    io.coeffs = mxGetDoubles(prhs[0]);
    io.mask   = mxGetDoubles(prhs[1]);
    // ------
    params.M     = M;
    params.wm    = mxGetDoubles(prhs[2]);
    params.fov   = (SizeBuffer)(fov);
    params.nhood = (SizeBuffer)(nhood);
    //=======================================================================================
    // Create the ouput
    plhs[0] = mxCreateNumericArray( 3, (SizeBuffer)(fov), mxDOUBLE_CLASS, mxREAL );
    io.disp = mxGetDoubles(plhs[0]);
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[4]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    SizeType N = (fov[0]) * (fov[1]) * (fov[2]);
    threader.setProcessSize( N, 20 );
    // Own values:
    threader.params = &params;
    threader.io     = &io;
    //=======================================================================================
    threader.threadedProcess( maxthreads, coeffs2localdispersion_process_fcn );
    //coeffs2localdispersion_process_fcn( (void*)(&threader) );
    //=======================================================================================
    return;
}

THFCNRET coeffs2localdispersion_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args                 = (ThArgs*)inargs;
    LDispesionIOData* io         = args->io;
    LDispesionParameters* params = args->params;
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads_thread(1);
    // ---------------------------------------------------------------
    // Create region iterators for the input coefficients volume
    // and for the mask image:
    dmriiters::NeighborhoodIterator coeffIt( 3, params->fov, params->M, io->coeffs, params->nhood );
    coeffIt.Begin();
    dmriiters::NeighborhoodIterator maskIt(  3, params->fov, 1,         io->mask,   params->nhood );
    maskIt.Begin();
    SizeType nneigh = coeffIt.GetNhoodSize();
    // ---------------------------------------------------------------
    // Create buffers where we will store the coefficients of each
    // neighboring voxel and at the center voxel:
    BufferType  neighbor = new ElementType[params->M];
    BufferType  voxel    = new ElementType[params->M];
    ElementType maskngh  = 0.0;
    // ---------------------------------------------------------------
    // Loop through the voxels
    IndexType start, end, i;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( maskIt.AbsolutePosition(start), coeffIt.AbsolutePosition(start), i=start; i<end; maskIt.Next(), coeffIt.Next(), ++i ){
            io->disp[i] = 0.0f;
            if( io->mask[i]<0.5 )
                continue;
            // Get the value of the center voxel and compute its norm:
            coeffIt.GetPixel(voxel);
            ElementType nrm0 = 0.0f;
            for( unsigned int j=0; j<params->M; ++j )
                nrm0 += (voxel[j])*(voxel[j])*(params->wm[j]);
            nrm0 = sqrt(nrm0);
            // Initialize the value of the dispersion
            ElementType dispersion = 0.0f;
            // For each neighbor...
            coeffIt.Rewind();
            maskIt.Rewind();
            SizeType  valid = 0;
            IndexType npos  = 0;
            bool proceed = true;
            while( proceed ){
                // Get the values of both the neighbor value
                // and the neighbor mask:
                maskIt.Play(&maskngh);
                proceed = coeffIt.Play(neighbor);
                // Only in case the mask value is non-zero, add
                // a new neighbor to the computations. Ignore
                // the center of the neighborhood
                if( (maskngh>0.5) && (npos!=(nneigh/2+1)) ){
                    // Add one more neighbor:
                    ++valid;
                    // Compute the dot product and the norm of
                    // the new neighbor:
                    ElementType dp  = 0.0f;
                    ElementType nrm = 0.0f;
                    for( unsigned int j=0; j<params->M; ++j ){
                        dp  += (voxel[j])*(neighbor[j])*(params->wm[j]);
                        nrm += (neighbor[j])*(neighbor[j])*(params->wm[j]);
                    }
                    nrm = sqrt(nrm);
                    dispersion += dp/(nrm*nrm0);
                }
                ++npos;
            }
            // Normalize:
            if( valid>0 )
                dispersion /= valid;
            else
                dispersion  = 1.0;
            // Set the final value:
            io->disp[i] = 1.0f - dispersion*dispersion;
        }
    }
    while( start < args->getN() );

    delete[] neighbor;
    delete[] voxel;
    blas_num_threads_thread(blas_threads);

    return (THFCNRET)NULL;
}

bool coeffs2localdispersion_callstack( const char* callerFunc )
{
    if( callerFunc==(char*)NULL )
        return false;
    else{
        if( !strcmp(callerFunc,"spectrum2fdi") )
            return true;
        if( !strcmp(callerFunc,"mapl2fdi") )
            return true;
        if( !strcmp(callerFunc,"hydidsi2fdi") )
            return true;
        if( !strcmp(callerFunc,"shodf2fdi") )
            return true;
    }
    return false;
}
