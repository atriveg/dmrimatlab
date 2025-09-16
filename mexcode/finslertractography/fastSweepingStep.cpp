/*==========================================================
 * fastSweepingStep.cpp
 *
 * This is a core function to fastSweeping, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2025 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
//#include "matrix.h"
#include "math.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../pixelIterators/iterators.h"

typedef struct FovDescription{
    unsigned int ndim;
    SizeBuffer fov;
    SizeType N;
} FovDescription;

void fastSweepingStep( const FovDescription&, const BufferType, const BufferType, const BufferType, const BufferType, BufferType, BufferType );

ElementType interpolateIncomingCost( const unsigned int, const SizeType, const IndexType, const BufferType, const BufferType, const BufferType );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed; in all cases ndim stands for the number of
     *          dimensions of the image):
     * 
     * prhs[0]:  costs, the map of arrival times, X_1 x X_2 x ... x X_ndim
     * prhs[1]:  dcosts, the map of directional costs, N x X_1 x X_2 x ... x X_ndim
     * prhs[2]:  mask, the mask of pixels to actually process, X_1 x X_2 x ... x X_ndim
     * prhs[3]:  neighbors, N x ndim, the neighbors used to interpolate each of the N directions
     * prhs[4]:  weights, N x ndim, the weigths to apply to each of the ndim neighbors
     *
     *  OUTPUTS:
     *
     * plhs[0]: costs, the updated map of costs X_1 x X_2 x ... x X_ndim
     * plhs[1]: dirs, the map of optimal (indexed) arrival directions X_1 x X_2 x ... x X_ndim
     * 
     * BEWARE:
     *
     *    - The calling function must ensure all pixels outside the mask are fixed to -1
     *    - The calling function must ensure all directional costs are non-negative
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:fastSweepingStep:callstack","This function should only be called from fastSweeping");
    else if( strcmp(callerFunc,"fastSweeping") )
        mexErrMsgIdAndTxt("MyToolbox:fastSweepingStep:callstack","This function should only be called from fastSweeping");
    //===========================  fastSweepingStep  ========================================
    if(nrhs!=5)
        mexErrMsgIdAndTxt("MyToolbox:fastSweepingStep:nrhs","Exactly 5 input arguments are required");
    if(nlhs!=2)
        mexErrMsgIdAndTxt("MyToolbox:fastSweepingStep:nlhs","Exactly 2 output argument is required");
    //=======================================================================================
    const unsigned int ndim = mxGetNumberOfDimensions(prhs[0]);
    SizeBuffer fov = new SizeType[ndim];
    for( unsigned int d=0; d<ndim; ++d )
        fov[d] = mxGetDimensions(prhs[0])[d];
    SizeType N = mxGetDimensions(prhs[1])[0];
    //=======================================================================================
    mwSize* arrSize = new mwSize[ndim];
    for( unsigned int d=0; d<ndim; ++d )
        arrSize[d] = (mwSize)fov[d];
    plhs[0] = mxCreateNumericArray( (mwSize)ndim, arrSize, mxDOUBLE_CLASS, mxREAL );
    plhs[1] = mxCreateNumericArray( (mwSize)ndim, arrSize, mxDOUBLE_CLASS, mxREAL );
    delete[] arrSize;
    //=======================================================================================
    FovDescription desc;
    desc.ndim = ndim;
    desc.fov  = fov;
    desc.N    = N;
    //=======================================================================================
    BufferType rcosts  = mxGetDoubles(prhs[0]);
    BufferType dcosts  = mxGetDoubles(prhs[1]);
    BufferType mask    = mxGetDoubles(prhs[2]);
    BufferType neighs  = mxGetDoubles(prhs[3]);
    BufferType weights = mxGetDoubles(prhs[4]);
    //=======================================================================================
    BufferType wcosts = mxGetDoubles(plhs[0]);
    BufferType odirs  = mxGetDoubles(plhs[1]);
    //=======================================================================================
    SizeType T = 1;
    for( unsigned int d=0; d<ndim; ++d )
        T *= fov[d];
    memcpy( wcosts, rcosts, T*sizeof(ElementType) );
    //=======================================================================================
    fastSweepingStep( desc, dcosts, mask, neighs, weights, wcosts, odirs );
    //=======================================================================================
    delete[] fov;
    //=======================================================================================
    return;
}

void fastSweepingStep(
    const FovDescription& desc,
    const BufferType dcosts,
    const BufferType mask,
    const BufferType neighs,
    const BufferType weights,
    BufferType wcosts,
    BufferType odirs )
{
    // ------------------------------------------------------------------------------------
    // Directional neighborhood iterator over the final costs (which are recursively updated):
    SizeBuffer radii = new SizeType[desc.ndim];
    for( unsigned int d=0; d<desc.ndim; ++d )
        radii[d] = 1;
    dmriiters::DirectionalNeighborhoodIterator wcosts_it( desc.ndim, desc.fov, 1, wcosts, radii );
    delete[] radii;
    // Infinite costs will be represented with -1, and the cost outside the FOV is always infinite:
    ElementType fillval = -1;
    wcosts_it.SetFillValue( &fillval );
    // ------------------------------------------------------------------------------------
    // Allocate buffers:
    // To store the entire neighborhood of the cost pixel being iterated, so that costs
    // coming from different directions can be interpolated:
    BufferType nhood = new ElementType[ wcosts_it.GetNhoodSize() ]; // The costs are scalar
    // To store the directional costs for the current pixel:
    BufferType dcost = new ElementType[ desc.N ];
    // ------------------------------------------------------------------------------------
    // Loop through the output costs in all possible causal and non-causal directions:
    for( wcosts_it.BeginDir(); !wcosts_it.EndDir(); wcosts_it.NextDir() ){
        for( wcosts_it.Begin(); !wcosts_it.End(); wcosts_it.Next() ){
            // Get the absolute position of the pixel being iterated within the
            // image buffer. We can use this value to index all images:
            IndexType apos = wcosts_it.AbsolutePosition();
            // If this position lays outside the mask, its cost must have been
            // fixed to -1 from the calling function. Nothing to do, so just
            // skip this pixel
            if( mask[apos] < 0.5 )
                continue;
            // If this is a seeding point, we should also keep on going without any
            // further processing, otherwise its value will be overriden and the
            // algorithm will never converge:
            if( wcosts[apos]==0 )
                continue;
            // Otherwise, we will need the neighborhood of this pixel, and also the set
            // of directional costs:
            wcosts_it.GetNeighborhood(nhood);
            memcpy( dcost, &(dcosts[apos*(desc.N)]), (desc.N)*sizeof(ElementType) );
            // Now, check all incoming directions and choose the one leading to
            // the smallest possible cost at the present pixel:
            ElementType cost = fillval;
            IndexType   odir = -1;
            ElementType incoming = -1.0;
            for( IndexType n=0; n<(desc.N); ++n ){
                // Get the cost at the n-th incoming direction by interpolation:
                incoming = interpolateIncomingCost( desc.ndim, desc.N, n, nhood, neighs, weights );
                // Proceed only if this cost is finite (i.e. positive)
                if(incoming>=0){
                    // Compute the new cost for the present direction
                    ElementType newCost = incoming + dcost[n];
                    // Update the cost at this voxel if the new cost is
                    // smaller than the previous one:
                    if( (cost<0.0) || (newCost<cost) ){
                        cost = newCost;
                        odir = n;
                    }
                }
            }
            // At this point, either the cost has been properly updated to the minimum
            // one reachable or it is still infinite because all its neighbors have
            // infinite cost.
            if( cost>=0.0 ){
                // We must update the cost in anycase:
                wcosts_it.SetPixel(&cost);
            }
            odirs[apos] = (ElementType)odir;
        }
    }
    // ------------------------------------------------------------------------------------
    // Delete allocated buffers:
    delete[] nhood;
    delete[] dcost;
    // ------------------------------------------------------------------------------------
    return;
}

ElementType interpolateIncomingCost( const unsigned int ndim, const SizeType N, const IndexType n, const BufferType nhood, const BufferType neighs, const BufferType weights )
{
    ElementType interpolated = 0.0;
    ElementType weightnorm   = 0.0;
    // Use as many neighbors as image dimensions:
    for( unsigned int d=0; d<ndim; ++d ){
        // What is the d-th neighbor for the n-th direction?
        IndexType neighbor = (IndexType)(neighs[n+N*d]);
        // If the value of this neighbor is infinite (i.e
        // negative), it cannot be used for interpolation
        ElementType value = nhood[neighbor];
        if(value<0.0)
            continue;
        // Otherwise, it will be taken into account with its
        // corresponding weight:
        ElementType weight = weights[n+N*d];
        // Actullay do interpolation:
        interpolated += value*weight;
        weightnorm   += weight;
    }
    if( weightnorm<1.0e-12 ){
        // If all the neighbors are infinite (i.e negative)
        // we cannot interpolate, and we can simply return
        // infinite (i.e negative)
        return -1.0;
    }
    else{
        // Otherwise, return the normalized weighted average:
        return (interpolated/weightnorm);
    }
}
