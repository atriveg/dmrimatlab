/*==========================================================
 * backTracing_.cpp
 *
 * This is a core function to backTracing.m, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2025 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <cstring>
#include <algorithm>
#include <vector>
#include "../mathsmex/mexToMathsTypes.h"
#include "../mathsmex/matrixCalculus.h"
#include "../threads/threadHelper.h"
#include "../pixelIterators/iterators.h"

#include <iostream>

#define MAXIMUM_NUMBER_OF_POINTS_PER_PATH 1e5

typedef struct FovDescription{
    unsigned int ndim;
    BufferType ijk2xyz;
    BufferType xyz2ijk;
    SizeBuffer fov;
} FovDescription;

typedef struct BackTracingArgs{
    SizeType N;
    BufferType  targets;
    ElementType step;
    ElementType Mlength;
    ElementType mcurv;
} BackTracingArgs;

typedef std::vector<ElementType> PointType;
typedef std::vector<PointType>   PathType;
typedef std::vector<ElementType> CostsPathType;

typedef struct BackTracingIO{
    BufferType     costs;
    BufferType     dirs;
    BufferType     mask;
    PathType*      paths;
    CostsPathType* pcosts;
    BufferType     stopconds;
} BackTracingIO;

typedef struct NewPathBuffers{
    BufferType  prev;
    BufferType  current;
    BufferType  next;
    BufferType  RK_k1;
    BufferType  RK_k2;
    BufferType  RK_k3;
    BufferType  RK_k4;
    BufferType  ijk;
    IndexBuffer neighs;
    BufferType  weights;
    dmriiters::NeighborhoodIterator* nhood;
    IndexBuffer position;
    IndexBuffer offset;
} NewPathBuffers;

class BackTracingThreader : public DMRIThreader
{
public:
    FovDescription* desc;
    BackTracingArgs* args;
    BackTracingIO* io;
};

THFCNRET backtracing_process_fcn( void* );

ElementType newPath( BackTracingThreader*, NewPathBuffers*, SizeType&, unsigned int&, PathType*, CostsPathType* );

int newPoint( BackTracingThreader*, NewPathBuffers*, ElementType&, ElementType& );

int dumbNewPoint( BackTracingThreader*, NewPathBuffers*, ElementType&, const ElementType& );

bool computeInterpolationFactors( const unsigned int, const SizeBuffer, BufferType, const BufferType, BufferType, IndexBuffer, BufferType );

void interpolateImage( const unsigned int, const SizeType, const BufferType, const IndexBuffer, const BufferType, BufferType );

void allocateNewPathBuffers( NewPathBuffers*, FovDescription*, BufferType );

void freeNewPathBuffers( NewPathBuffers* );

ElementType computeCurvatureRadius( const unsigned int, const BufferType, const BufferType, const BufferType );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed; in all cases ndim stands for the number of
     *          dimensions of the image):
     * 
     * prhs[0]:  costs, the map of arrival times, X_1 x X_2 x ... x X_ndim
     * prhs[1]:  dirs, the map of arrival directions, ndim x X_1 x X_2 x ... x X_ndim
     * prhs[2]:  mask, the mask of pixels to actually process, X_1 x X_2 x ... x X_ndim. MUST be 0 or 1, in double format
     * prhs[3]:  targets, ndim x N, the physical target points (xyz space) from which paths are traced
     * prhs[4]:  ijk2xyz, (ndim+1)x(ndim+1), the homogeneous matrix to go from pixel coords to physical coords
     * prhs[5]:  xyz2ijk, (ndim+1)x(ndim+1), the homogeneous matrix to go from physical coords to pixel coords
     * prhs[6]:  step, 1x1, the base step for Runge-Kutta integration
     * prhs[7]:  Mlength, 1x1, the maximum length of a path to be accepted
     * prhs[8]:  mcurv, 1x1, the minimum curvature radius before a path is aborted
     * prhs[9]:  maxthreads, 1x1, the maximum number of threads to use
     *
     *  OUTPUTS:
     *
     * plhs[0]: paths, a Nx1 cell array, each field is a ndim x M_n array describing the backtraced path
     * plhs[1]: pcosts, a Nx1 cell array, each field is a 1 x M_n array describing the costs at each point in the path
     * plhs[2]: stopconds, a Nx1 array of doublesm each value corresponding to the stop condition associated to each path
     *
     * BEWARE: the calling function must ensure that the mask is represented in ElementType format, and
     *         it has to evaluate to 1 (inside) or 0 (outside)
     * 
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:backTracing_:callstack","This function should only be called from backTracing");
    else if( strcmp(callerFunc,"backTracing") )
        mexErrMsgIdAndTxt("MyToolbox:backTracing_:callstack","This function should only be called from backTracing");
    //===========================  fastSweepingStep  ========================================
    if(nrhs!=10)
        mexErrMsgIdAndTxt("MyToolbox:backTracing_:nrhs","Exactly 10 input arguments are required");
    if(nlhs!=3)
        mexErrMsgIdAndTxt("MyToolbox:backTracing_:nlhs","Exactly 3 output arguments are required");
    //=======================================================================================
    const unsigned int ndim = mxGetNumberOfDimensions(prhs[0]);
    SizeBuffer fov = new SizeType[ndim];
    for( unsigned int d=0; d<ndim; ++d )
        fov[d] = mxGetDimensions(prhs[0])[d];
    //=======================================================================================
    FovDescription desc;
    desc.ndim = ndim;
    desc.ijk2xyz = mxGetDoubles(prhs[4]);
    desc.xyz2ijk = mxGetDoubles(prhs[5]);
    desc.fov = fov;
    //=======================================================================================
    SizeType N = mxGetN(prhs[3]);
    BackTracingArgs args;
    args.N = N;
    args.targets = mxGetDoubles(prhs[3]);
    args.step = mxGetScalar(prhs[6]);
    args.Mlength = mxGetScalar(prhs[7]);
    args.mcurv = mxGetScalar(prhs[8]);
    //=======================================================================================
    BackTracingIO io;
    io.costs = mxGetDoubles(prhs[0]);
    io.dirs = mxGetDoubles(prhs[1]);
    io.mask = mxGetDoubles(prhs[2]);
    PathType* paths = new PathType[N];
    CostsPathType* pcosts = new CostsPathType[N];
    io.paths = paths;
    io.pcosts = pcosts;
    //=======================
    plhs[2] = mxCreateDoubleMatrix( 1, N, mxREAL );
    io.stopconds = mxGetDoubles( plhs[2] );
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[9]) );
    unsigned long chunksz = 5;
    if(N<10)
        chunksz = 1;
    //=======================================================================================
    BackTracingThreader threader;
    threader.setProcessSize( N, chunksz );
    threader.desc = &desc;
    threader.args = &args;
    threader.io = &io;
    threader.threadedProcess( maxthreads, backtracing_process_fcn );
    //=======================================================================================
    delete[] fov;
    //=======================================================================================
    plhs[0] = mxCreateCellMatrix(N,1);
    for( IndexType n=0; n<N; ++n ){
        PathType path = io.paths[n];
        if(path.size()>0){
            mxArray* current = mxCreateDoubleMatrix( ndim, path.size(), mxREAL );
            IndexType p = 0;
            std::for_each(
                path.begin(),
                path.end(),
                [ndim,current,&p](const PointType& point){
                    for( unsigned int d=0; d<ndim; ++d )
                        mxGetDoubles(current)[d+ndim*p] = point[d];
                    ++p;
                }
            );
            mxSetCell( plhs[0], n, current );
        }
    }
    delete[] paths;
    //=======================================================================================
    plhs[1] = mxCreateCellMatrix(N,1);
    for( IndexType n=0; n<N; ++n ){
        CostsPathType costpath = io.pcosts[n];
        if(costpath.size()>0){
            mxArray* current = mxCreateDoubleMatrix( 1, costpath.size(), mxREAL );
            IndexType p = 0;
            std::for_each(
                costpath.begin(),
                costpath.end(),
                [current,&p](const ElementType& ccost){
                    mxGetDoubles(current)[p] = ccost;
                    ++p;
                }
            );
            mxSetCell( plhs[1], n, current );
        }
    }
    delete[] pcosts;
    //=======================================================================================
    return;
}

/**
 *  Threaded process function
 */
THFCNRET backtracing_process_fcn( void* inargs )
{
    // ----
    BackTracingThreader* btargs = (BackTracingThreader*)inargs;
    FovDescription* desc = btargs->desc;
    BackTracingArgs* args = btargs->args;
    BackTracingIO* io = btargs->io;
    // ----
    NewPathBuffers buffers;
    allocateNewPathBuffers( &buffers, desc, io->costs );
    // ----
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads(1);
    // ----
    IndexType start = 0; // First index of data blocks to process
    IndexType end   = 0; // Post-last index of data blocks to process
    do{
        // Claim a new block of data to process within this thread:
        btargs->claimNewBlock( &start, &end );
        // While there is still data not processed by any thread...
        for( IndexType i=start; i<end; ++i ){
            // Initialize the current path to an empty one
            (io->paths[i]).clear();
            (io->pcosts[i]).clear();
            // Find the i-th seeding point:
            memcpy( buffers.current, &(args->targets[i*(desc->ndim)]), (desc->ndim)*sizeof(ElementType) );
            // Try to trace a new path:
            SizeType numberOfPoints = 0;
            unsigned int causeOfStop = 0;
            ElementType length = newPath( (BackTracingThreader*)inargs, &buffers, numberOfPoints, causeOfStop, &(io->paths[i]), &(io->pcosts[i]) );
            // Store the stopping condition in the output buffer:
            io->stopconds[i] = causeOfStop;
        }
    }
    while( start < btargs->getN() );
    // ----
    // Revert BLAS threads usage to its default:
    blas_num_threads(blas_threads);
    // ----
    freeNewPathBuffers( &buffers );
    // ----
    return (THFCNRET)NULL;
}

/**
 *
 *  Backtraces a new complete streamline until a stop condition is reached,
 *  and this stop condition is stored in causeOfStop, with the following
 *  meaning:
 *     0: normal stop condition, reached seeding region (low cost reached)
 *     1: normal stop condition, reached seeding region (low jump reached)
 *     2: streamline reached a point outside the FoV
 *     3: streamline reached a point outside the mask
 *     4: maximum number of points per stream was reached
 *     5: maximum streamline length reached
 *     6: local curvature radius of the streamline below allowed minimum
 *
 */
ElementType newPath( BackTracingThreader* inargs, NewPathBuffers* buffers, SizeType& numberOfPoints, unsigned int& causeOfStop, PathType* path, CostsPathType* pcost )
{
    // ---------------------------------------
    BackTracingThreader* btargs = (BackTracingThreader*)inargs;
    FovDescription* desc = btargs->desc;
    BackTracingArgs* args = btargs->args;
    BackTracingIO* io = btargs->io;
    // ---------------------------------------
    memcpy( buffers->next, buffers->current, (desc->ndim)*sizeof(ElementType) );
    memcpy( buffers->prev, buffers->current, (desc->ndim)*sizeof(ElementType) );
    // ---------------------------------------
    PointType point(desc->ndim); // We will store here the succesive points in the path
    // ---------------------------------------
    // Tolerances used for the stopping criteria:
    ElementType fCostTol   = 0.0001;            // Minimum cost until we consider we reached the seeding region
    ElementType mjump      = 0.01*(args->step); // Minimum jump until we consider we reached the seeding region
    // ---------------------------------------
    // Initialize the algortihm:
    causeOfStop = 0;
    numberOfPoints = 0;
    ElementType dpath = 0.0;
    // ---------------------------------------
    // Find the new path point by point:
    do{
        // Actual implementation of the Runge-Kutta integration step:
        ElementType curvature; // The curcature in the current target
        ElementType cost;      // The interpolated cost in the current target
        int status = newPoint( inargs, buffers, curvature, cost );
        // If the status is 2 or 3, it means that the current target is
        // either outside the image buffer or outside the image mask.
        // In both cases, the path should be aborted without adding the
        // current target to the path
        if( status == 2 ){
            causeOfStop = status; // target out of bounds
            return(dpath);
        }
        if( status == 3 ){
            causeOfStop = status; // target outside the mask
            return(dpath);
        }
        if( status == -1 ){
            // Runge-Kutta step produced either:
            //   - An out-of-bounds point
            //   - An out-of-mask pixel
            //   - A point whose cost has increased w.r.t. the previous one
            // Let's try to fix it by explicitly moving towards the
            // neighbor of the curren pixel that has minimal cost (out-of
            // -bounds and out-of-mask pixels have infinite cost)
            causeOfStop = dumbNewPoint( inargs, buffers, curvature, cost );
            if( causeOfStop != 0 )
                return(dpath);
        }
        // ----------------------------
        // Otherwise, the current target must be included in the path...
        for( unsigned int d=0; d<desc->ndim; ++d )
            point[d] = buffers->current[d];
        path->push_back(point);
        pcost->push_back(cost);
        // ... so that we have a new point in the path:
        ++numberOfPoints;
        ElementType dist = 0.0;
        for( unsigned int d=0; d<desc->ndim; ++d )
            dist += (buffers->next[d]-buffers->current[d])*(buffers->next[d]-buffers->current[d]);
        dist   = sqrt(dist);
        dpath += dist;
        // ----------------------------
        // Iterate:
        memcpy( buffers->prev, buffers->current, (desc->ndim)*sizeof(ElementType) ); // No need to copy ndim+1 since it is 1
        memcpy( buffers->current, buffers->next, (desc->ndim)*sizeof(ElementType) ); // No need to copy ndim+1 since it is 1
        // ----------------------------
        // If this is the very first point added to the path, we can use it
        // to establish the tolerances for stopping criteria:
        if(numberOfPoints==1){
            fCostTol = std::max( cost/1000, 1.0e-3 );
            mjump    = std::max( 0.01*dist, 1.0e-3 );
        }
        // ----------------------------
        // We consider the path is all done if we have reached the seeding
        // region, i.e if the cost is small enough:
        if( cost < fCostTol ) // Exit value 0
            return(dpath);
        // ----------------------------
        // Check if additional stop conditions occur:
        if( dist < mjump ) // Note this is also a "good exit", meaning we are close to the seeding region
            causeOfStop = 1;
        if( numberOfPoints >= MAXIMUM_NUMBER_OF_POINTS_PER_PATH )
            causeOfStop = 4;
        if( dpath > args->Mlength )
            causeOfStop = 5;
        if( (curvature < args->mcurv) && (numberOfPoints>2) )
            causeOfStop = 6;
        // ----------------------------
    }
     while(causeOfStop==0);

     return(dpath);
}


/**
 * Inserts a new point in the streamline using a RK4 Runge-Kutta
 * integration step of the arrival directions map, which is linearly
 * interpolated. Returns an error condition if:
 *  - the current point is outside the image FoV (2)
 *  - the current point is outside the image mask (3)
 *  - any of the intermediate RK4 estimates is outside the FoV (-1)
 */
int newPoint( BackTracingThreader* inargs, NewPathBuffers* buffers, ElementType& curvature, ElementType& cost )
{
    // ---------------------------------------
    BackTracingThreader* btargs = (BackTracingThreader*)inargs;
    FovDescription* desc = btargs->desc;
    BackTracingArgs* args = btargs->args;
    BackTracingIO* io = btargs->io;
    // ---------------------------------------
    // Compute the interpolation factors for the current target:
    bool inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->current, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    // If we are outside the buffer, just return:
    if(!inside)
        return 2;
    // ---------------------------------------
    // Get the interpolated mask value to check if we are at a masked value:
    ElementType maskval;
    interpolateImage( desc->ndim, 1, io->mask, buffers->neighs, buffers->weights, &maskval );
    if( maskval<0.5 )
        return 3;
    // ---------------------------------------
    // If we reach here, the current target is a good one. We
    // can compute its interpolated cost:
    interpolateImage( desc->ndim, 1, io->costs, buffers->neighs, buffers->weights, &cost );
    curvature = 2*(args->mcurv); // temporary
    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    // From this point, we implement the RK4 method to find the next point
    // in the streamline:
    // k1 is just the estimate from the current position:
    interpolateImage( desc->ndim, desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->RK_k1 );
    // Move h/2 towards k1 to find a first estimate of the middle point:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + 0.5*(args->step)*(buffers->RK_k1[d]);
    // k2 will be computed as the interpolated value at the first estimate
    // for the middle point:
    inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    // The management of this condition is delegated to the next step:
    if(!inside)
        return -1;
    // Actually interpolate:
    interpolateImage( desc->ndim, desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->RK_k2 );
    // Move h/2 towards k2 to find a better estimate of the middle point:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + 0.5*(args->step)*(buffers->RK_k2[d]);
    // k3 will be computed as the interpolated value at the refined estimate
    // for the middle point:
    inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    // The management of this condition is delegated to the next step:
    if(!inside)
        return -1;
    // Actually interpolate:
    interpolateImage( desc->ndim, desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->RK_k3 );
    // Move h towards k3 to find an estimate of the final point:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + (args->step)*(buffers->RK_k3[d]);
    // k4 will be computed as the interpolated value at the estimate for
    // the final point:
    inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    // The management of this condition is delegated to the next step:
    if(!inside)
        return -1;
    // Actually interpolate:
    interpolateImage( desc->ndim, desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->RK_k4 );
    // If we reach here, we have the four estimates k1, k2, k3 and k4, and
    // the actual step to jump is a weighted average of them:
    ElementType w1 = (args->step)*0.16666666666666666667;
    ElementType w2 = (args->step)*0.33333333333333333333;
    ElementType w3 = (args->step)*0.33333333333333333333;
    ElementType w4 = (args->step)*0.16666666666666666667;
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + w1*(buffers->RK_k1[d]) + w2*(buffers->RK_k2[d]) + w3*(buffers->RK_k3[d]) + w4*(buffers->RK_k4[d]);
    // Last check if the new point is outside the buffer:
    inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    if(!inside)
        return -1;
    interpolateImage( desc->ndim, 1, io->mask, buffers->neighs, buffers->weights, &maskval );
    if( maskval<0.5 )
        return -1;
    // If we reach here, we have computed a new point which is both
    // inside the FoV and inside the mask. We can re-compute the cost
    // accordingly:
    ElementType newcost;
    interpolateImage( desc->ndim, 1, io->costs, buffers->neighs, buffers->weights, &newcost );
    // If the cost has increased, something went wrong:
    if( newcost>cost )
        return -1;
    // ---------------------------------------
    // It only remains to compute the curvature radius:
    curvature = computeCurvatureRadius( desc->ndim, buffers->prev, buffers->current, buffers->next );
    if( curvature<0.0 ) // Meaning: we have repeated points
        curvature = 2*(args->mcurv);
    // ---------------------------------------
    return 0;
}

/**
 * Inserts a new point in the streamline based only on the costs of
 * neighboring voxels, ignoring the arrival directions map (hence, no
 * numerical integration is performed). This is only used when the RK4
 * method produces a result either outside the FoV or outside the
 * image mask
 */
int dumbNewPoint( BackTracingThreader* inargs, NewPathBuffers* buffers, ElementType& curvature, const ElementType& cost )
{
    // ---------------------------------------
    BackTracingThreader* btargs = (BackTracingThreader*)inargs;
    FovDescription* desc = btargs->desc;
    BackTracingArgs* args = btargs->args;
    BackTracingIO* io = btargs->io;
    // ---------------------------------------
    // From the current physical point, get the corresponding pixel point:
    buffers->current[desc->ndim] = 1.0;
    mataux::multiplyMxArrays( desc->xyz2ijk, buffers->current, buffers->ijk, desc->ndim+1, desc->ndim+1, 1 );
    // Round this value and make sure it is in-bounds:
    for( unsigned int d=0; d<desc->ndim; ++d ){
        buffers->position[d] = (IndexType)(   round( buffers->ijk[d] )   );
        buffers->position[d] = ( buffers->position[d] >= 0           ? buffers->position[d] : 0              );
        buffers->position[d] = ( buffers->position[d] < desc->fov[d] ? buffers->position[d] : desc->fov[d]-1 );
    }
    // ---------------------------------------
    // Place the iterator at the current position so that we can retrieve
    // the arrival costs in its 1x1x...x1 neighborhood:
    buffers->nhood->SetIndex(buffers->position);
    // ---------------------------------------
    // Iterate through the neighborhood to find the voxel with
    // the minimum arrival cost:
    ElementType mval = std::numeric_limits<ElementType>::infinity();
    ElementType cval;
    IndexType   pos = 0;
    IndexType   mpos = 0;
    buffers->nhood->Rewind();
    while( buffers->nhood->Play(&cval) ){
        if( cval<mval ){
            mval = cval;
            mpos = pos;
        }
        ++pos;
    }
    // Now mpos contains the neighborhood position
    // corresponding to the minimum cost
    // ---------------------------------------
    // Get the offset associated to the neighbor
    // with minimum cost...
    buffers->nhood->GetOffset( mpos, buffers->offset );
    // ... so that the overall position of such
    // neighbor is:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->ijk[d] = buffers->position[d] + buffers->offset[d];
    // ---------------------------------------
    // Map this neighbor back to the physical space:
    buffers->ijk[desc->ndim] = 1.0;
    mataux::multiplyMxArrays( desc->ijk2xyz, buffers->ijk, buffers->next, desc->ndim+1, desc->ndim+1, 1 );
    // Save this position in case we need it later:
    memcpy( buffers->RK_k1, buffers->next, (desc->ndim)*sizeof(ElementType) );
    // And move "h" towards this point:
    ElementType jumpnorm = 0.0;
    for( unsigned int d=0; d<desc->ndim; ++d ){
        buffers->next[d] -= buffers->current[d];
        jumpnorm += (buffers->next[d])*(buffers->next[d]);
    }
    jumpnorm = (args->step)/sqrt(jumpnorm);
    for( unsigned int d=0; d<desc->ndim; ++d ){
        buffers->next[d] *= jumpnorm;
        buffers->next[d] += buffers->current[d];
    }
    // ---------------------------------------
    int inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    if( inside ){
        // There's still a chance that this jump has led us outside
        // the image mask:
        ElementType maskval;
        interpolateImage( desc->ndim, 1, io->mask, buffers->neighs, buffers->weights, &maskval );
        // Or it might have increased the cost:
        ElementType newcost;
        interpolateImage( desc->ndim, 1, io->costs, buffers->neighs, buffers->weights, &newcost );
        if(   (maskval<0.5) || (newcost>cost)   ) // Revert the value of buffers->next to a pixel-based position
            memcpy( buffers->next, buffers->RK_k1, (desc->ndim)*sizeof(ElementType) );
    }
    else // Something went terribly wrong
        return -1;
    // ---------------------------------------
    // It only remains to compute the curvature radius:
    curvature = computeCurvatureRadius( desc->ndim, buffers->prev, buffers->current, buffers->next );
    if( curvature<0.0 ) // Meaning: we have repeated points
        curvature = 2*(args->mcurv);
    // ---------------------------------------
    return 0;
}

/**
 * Computes the curvature radius of the streamline at a given point
 * from its previous and next neighbors in that streamline
 */
ElementType computeCurvatureRadius( const unsigned int ndim, const BufferType prev, const BufferType current , const BufferType next )
{
    ElementType kl1 = 0.0;
    ElementType kl2 = 0.0;
    for( unsigned int d=0; d<ndim; ++d ){
        kl1 += (current[d]-next[d])*(current[d]-next[d]);
        kl2 += (prev[d]-current[d])*(prev[d]-current[d]);
    }
    kl1 = sqrt(kl1);
    kl2 = sqrt(kl2);
    if( (kl1<1e-12) || (kl2<1e-12) )
        return -1.0;
    ElementType norm = 2.0/(kl1+kl2);
    ElementType curvature = 0.0;
    for( unsigned int d=0; d<ndim; ++d ){
        ElementType tmp = (prev[d]-current[d])/kl2 - (current[d]-next[d])/kl1;
        tmp *= norm;
        curvature += (tmp*tmp);
    }
    curvature = sqrt(curvature);
    if( curvature<1.0e-6 )
        return std::numeric_limits<ElementType>::infinity();
    else
        return 1.0/curvature;
}

/**
 * Computes the indices and associated weights for linear linterpolation
 * based on the (physical) coordinates of a pixel. Buffer sizes:
 *
 * fov:     ndim
 * pixel:   ndim
 * neighs:  2^ndim = 1 << ndim
 * weights: 2^ndim = 1 << ndim
 */
bool computeInterpolationFactors( const unsigned int ndim, const SizeBuffer fov, BufferType pixelxyz, const BufferType xyz2ijk, BufferType pixelijk, IndexBuffer neighs, BufferType weights )
{
    // Transform the pixel in physical coordinates xyz to pixel
    // coordinates ijk using homogeneous coordinates:
    pixelxyz[ndim] = 1.0;
    mataux::multiplyMxArrays( xyz2ijk, pixelxyz, pixelijk, ndim+1, ndim+1, 1 );
    // Check if the pixel is in-bounds:
    ElementType EPS_ = 1.0e-6; // Used to ensure that points in the limit of the FoV are not discarded (very rarely used)
    for( unsigned int d=0; d<ndim; ++d ){
        if( (pixelijk[d]<0.0-EPS_) || (pixelijk[d]>(ElementType)(fov[d]-1)+EPS_) )
            return false;
        pixelijk[d] = ( pixelijk[d]>=0.0                     ? pixelijk[d] : 0.0                     );
        pixelijk[d] = ( pixelijk[d]<=(ElementType)(fov[d]-1) ? pixelijk[d] : (ElementType)(fov[d]-1) );
    }
    // Compute the number of neighbors as a function of the number
    // of image dimensions, i.e. 2^ndim:
    SizeType nneighs = ( 1 << ndim );
    // Iterate trough these neighbors:
    for( IndexType i=0; i<nneighs; ++i ){
        weights[i] = 1.0;
        neighs[i]  = 0;
        SizeType stride  = 1;
        for( unsigned int d=0; d<ndim; ++d ){
            IndexType lower = (IndexType)(pixelijk[d]);
            if( (i>>d) & 1 ) { // This neighbor is in the "ceil" for dimension d
                IndexType upper = ( lower<(fov[d]-1) ? lower+1 : lower );
                weights[i] *= (pixelijk[d]-lower);
                neighs[i] += upper*stride;
            }
            else{ // This neighbor is in the "floor" for dimension d
                weights[i] *= (1.0-pixelijk[d]+lower);
                neighs[i] += lower*stride;
            }
            stride *= (fov[d]);
        }
    }
    return true;
}

/**
 * Uses the interpolation weigths and neighbors determined with computeInterpolationFactors
 * to interpolate N-component vector images. Any pixel containing infinite values
 * is ignored from the computation. If all neighbors are infinite, the output
 * will also be infinite.
 */
void interpolateImage( const unsigned int ndim, const SizeType N, const BufferType buffer, const IndexBuffer neighs, const BufferType weights, BufferType interpolated )
{
    SizeType nneighs = ( 1 << ndim );
    for( IndexType n=0; n<N; ++n )
        interpolated[n] = 0.0;
    ElementType trueweight = 0.0;
    for( IndexType i=0; i<nneighs; ++i ){
        bool anyinf = false;
        for( IndexType n=0; n<N; ++n )
            anyinf = ( anyinf || std::isinf(buffer[n+N*neighs[i]]) );
        if(!anyinf){
            trueweight += weights[i];
            for( IndexType n=0; n<N; ++n )
                interpolated[n] += weights[i]*(buffer[n+N*neighs[i]]);
        }        
    }
    if(trueweight==0.0){
        for( IndexType n=0; n<N; ++n )
            interpolated[n] = std::numeric_limits<ElementType>::infinity();
    }
    else{
        for( IndexType n=0; n<N; ++n )
            interpolated[n] /= trueweight;
    }
    return;
}

void allocateNewPathBuffers( NewPathBuffers* buffers, FovDescription* desc, BufferType costs )
{
    const unsigned int ndim = desc->ndim;
    buffers->prev = new ElementType[ndim+1];    // Plus 1 because we use homogeneous coords
    buffers->current = new ElementType[ndim+1]; // Plus 1 because we use homogeneous coords
    buffers->next = new ElementType[ndim+1];    // Plus 1 because we use homogeneous coords
    buffers->RK_k1 = new ElementType[ndim];
    buffers->RK_k2 = new ElementType[ndim];
    buffers->RK_k3 = new ElementType[ndim];
    buffers->RK_k4 = new ElementType[ndim];
    buffers->ijk = new ElementType[ndim+1];     // Plus 1 because we use homogeneous coords
    SizeType nneighs = ( 1 << ndim );           // 2^ndim
    buffers->neighs = new IndexType[nneighs];
    buffers->weights = new ElementType[nneighs];
    // Prepare the region iterator over the costs image, which
    // is used only if the RK4 method fails:
    SizeBuffer radii = new SizeType[ndim];
    for( unsigned int d=0; d<ndim; ++d )
        radii[d] = 1;
    buffers->nhood = new dmriiters::NeighborhoodIterator( ndim, desc->fov, 1, costs, radii );
    delete[] radii;
    ElementType fill = std::numeric_limits<ElementType>::infinity();
    buffers->nhood->SetFillValue( &fill );
    buffers->nhood->Begin();
    buffers->position = new IndexType[ndim];
    buffers->offset   = new IndexType[ndim];
    return;
}

void freeNewPathBuffers( NewPathBuffers* buffers )
{
    delete[] buffers->prev;
    delete[] buffers->current;
    delete[] buffers->next;
    delete[] buffers->RK_k1;
    delete[] buffers->RK_k2;
    delete[] buffers->RK_k3;
    delete[] buffers->RK_k4;
    delete[] buffers->ijk;
    delete[] buffers->neighs;
    delete[] buffers->weights;
    delete buffers->nhood;
    delete[] buffers->position;
    delete[] buffers->offset;
    return;
}

