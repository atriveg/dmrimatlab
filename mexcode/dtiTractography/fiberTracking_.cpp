/*==========================================================
 * fiberTracking_.cpp
 *
 * This is a core function to backTracing.m, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2025 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
//#include "matrix.h"
#include "math.h"
#include <cstring>
#include <algorithm>
#include <vector>
#include <cmath>
#include "../mathsmex/mexToMathsTypes.h"
#include "../mathsmex/matrixCalculus.h"
#include "../threads/threadHelper.h"

#define MAXIMUM_NUMBER_OF_POINTS_PER_PATH 1e5

typedef struct FovDescription{
    unsigned int ndim;
    BufferType ijk2xyz;
    BufferType xyz2ijk;
    SizeBuffer fov;
} FovDescription;

typedef struct BackTracingArgs{
    SizeType N;
    BufferType  seeds;
    ElementType threshold;
    ElementType step;
    ElementType Mlength;
    ElementType mcurv;
} BackTracingArgs;

typedef std::vector<ElementType> PointType;
typedef std::vector<PointType>   PathType;
typedef std::vector<ElementType> ScalarPathType;

typedef struct BackTracingIO{
    BufferType      scalar;
    BufferType      dirs;
    BufferType      mask;
    PathType*       paths;
    ScalarPathType* pvals;
    BufferType      stopconds;
} BackTracingIO;

typedef struct NewPathBuffers{
    BufferType  prev;
    BufferType  current;
    BufferType  next;
    BufferType  direction; // Used for deambiguation of the sign of the eigenvector
    BufferType  RK_k1;
    BufferType  RK_k2;
    BufferType  RK_k3;
    BufferType  RK_k4;
    BufferType  ijk;
    IndexBuffer neighs;
    BufferType  weights;
} NewPathBuffers;

class BackTracingThreader : public DMRIThreader
{
public:
    FovDescription* desc;
    BackTracingArgs* args;
    BackTracingIO* io;
};

THFCNRET backtracing_process_fcn( void* );

ElementType newPath( BackTracingThreader*, NewPathBuffers*, SizeType&, unsigned int&, PathType*, ScalarPathType*, const int );

int newPoint( BackTracingThreader*, NewPathBuffers*, ElementType&, ElementType& );

bool computeInterpolationFactors( const unsigned int, const SizeBuffer, BufferType, const BufferType, BufferType, IndexBuffer, BufferType );

void interpolateScalarImage( const unsigned int, const BufferType, const IndexBuffer, const BufferType, BufferType );

void interpolateVectorImage( const unsigned int, const BufferType, const IndexBuffer, const BufferType, const BufferType, BufferType );


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
     * prhs[0]:  scalar, the map of scalar indices used to stop fiber tracking, X_1 x X_2 x ... x X_ndim
     * prhs[1]:  dirs, a directions volume with size ndim x X_1 x X_2 x ... x X_ndim
     * prhs[2]:  mask, the mask of pixels to actually process, X_1 x X_2 x ... x X_ndim. MUST be 0 or 1, in double format
     * prhs[3]:  seeds, ndim x N, the physical seeding points (xyz space) from which paths are traced
     * prhs[4]:  ijk2xyz, (ndim+1)x(ndim+1), the homogeneous matrix to go from pixel coords to physical coords
     * prhs[5]:  xyz2ijk, (ndim+1)x(ndim+1), the homogeneous matrix to go from physical coords to pixel coords
     * prhs[6]:  threshold, 1x1, the lower threshold applied to scalar to determine the stop condition
     * prhs[7]:  step, 1x1, the base step for Runge-Kutta integration
     * prhs[8]:  Mlength, 1x1, the maximum length of a path to be accepted
     * prhs[9]:  mcurv, 1x1, the minimum curvature radius before a path is aborted
     * prhs[10]: maxthreads, 1x1, the maximum number of threads to use
     *
     *  OUTPUTS:
     *
     * plhs[0]: paths, a Nx1 cell array, each field is a ndim x M_n array describing the backtraced path
     * plhs[1]: pvals, a Nx1 cell array, each field is a 1 x M_n array describing the scalar values at each point in the path
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
        mexErrMsgIdAndTxt("MyToolbox:fiberTracking_:callstack","This function should only be called from dti2tractography");
    else if( strcmp(callerFunc,"dti2tractography") )
        mexErrMsgIdAndTxt("MyToolbox:fiberTracking_:callstack","This function should only be called from dti2tractography");
    //===========================  fastSweepingStep  ========================================
    if(nrhs!=11)
        mexErrMsgIdAndTxt("MyToolbox:fiberTracking_:nrhs","Exactly 11 input arguments are required");
    if(nlhs!=3)
        mexErrMsgIdAndTxt("MyToolbox:fiberTracking_:nlhs","Exactly 3 output arguments are required");
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
    args.seeds = mxGetDoubles(prhs[3]);
    args.threshold = mxGetScalar(prhs[6]);
    args.step = mxGetScalar(prhs[7]);
    args.Mlength = mxGetScalar(prhs[8]);
    args.mcurv = mxGetScalar(prhs[9]);
    //=======================================================================================
    BackTracingIO io;
    io.scalar = mxGetDoubles(prhs[0]);
    io.dirs  = mxGetDoubles(prhs[1]);
    io.mask = mxGetDoubles(prhs[2]);
    PathType* paths = new PathType[2*N];             // Two streamlines per seed
    ScalarPathType* pvals = new ScalarPathType[2*N]; // Two streamlines per seed
    io.paths = paths;
    io.pvals = pvals;
    //=======================
    plhs[2] = mxCreateDoubleMatrix( 1, 2*N, mxREAL );
    io.stopconds = mxGetDoubles( plhs[2] );
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[10]) );
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
    plhs[0] = mxCreateCellMatrix(2*N,1); // Two streamlines per seed
    for( IndexType n=0; n<2*N; ++n ){
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
    plhs[1] = mxCreateCellMatrix(2*N,1); // Two streamlines per seed
    for( IndexType n=0; n<2*N; ++n ){
        ScalarPathType scpath = io.pvals[n];
        if(scpath.size()>0){
            mxArray* current = mxCreateDoubleMatrix( 1, scpath.size(), mxREAL );
            IndexType p = 0;
            std::for_each(
                scpath.begin(),
                scpath.end(),
                [current,&p](const ElementType& scval){
                    mxGetDoubles(current)[p] = scval;
                    ++p;
                }
            );
            mxSetCell( plhs[1], n, current );
        }
    }
    return;
    delete[] pvals;
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
    allocateNewPathBuffers( &buffers, desc, io->scalar );
    // ----
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads_thread(1);
    // ----
    IndexType start = 0; // First index of data blocks to process
    IndexType end   = 0; // Post-last index of data blocks to process
    do{
        // Claim a new block of data to process within this thread:
        btargs->claimNewBlock( &start, &end );
        // While there is still data not processed by any thread...
        for( IndexType i=start; i<end; ++i ){
            // --------------------------------------------------
            // Initialize the current paths to an empty one
            (io->paths[2*i]).clear();
            (io->pvals[2*i]).clear();
            // Find the i-th seeding point:
            memcpy( buffers.current, &(args->seeds[i*(desc->ndim)]), (desc->ndim)*sizeof(ElementType) );
            // Try to trace a new path for the first direction:
            for(unsigned int d=0; d<desc->ndim; ++d )
                buffers.direction[d] = 0.0;
            SizeType numberOfPoints = 0;
            unsigned int causeOfStop = 0;
            ElementType length = newPath( (BackTracingThreader*)inargs, &buffers, numberOfPoints, causeOfStop, &(io->paths[2*i]), &(io->pvals[2*i]), 1 );
            // Store the stopping condition in the output buffer:
            io->stopconds[2*i] = causeOfStop;
            // --------------------------------------------------
            // Initialize the current paths to an empty one
            (io->paths[2*i+1]).clear();
            (io->pvals[2*i+1]).clear();
            // Find the i-th seeding point:
            memcpy( buffers.current, &(args->seeds[i*(desc->ndim)]), (desc->ndim)*sizeof(ElementType) );
            // Try to trace a new path for the second direction:
            for(unsigned int d=0; d<desc->ndim; ++d )
                buffers.direction[d] = 0.0;
            numberOfPoints = 0;
            causeOfStop = 0;
            length = newPath( (BackTracingThreader*)inargs, &buffers, numberOfPoints, causeOfStop, &(io->paths[2*i+1]), &(io->pvals[2*i+1]), -1 );
            // Store the stopping condition in the output buffer:
            io->stopconds[2*i+1] = causeOfStop;
            // --------------------------------------------------
        }
    }
    while( start < btargs->getN() );
    // ----
    // Revert BLAS threads usage to its default:
    blas_num_threads_thread(blas_threads);
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
ElementType newPath( BackTracingThreader* inargs, NewPathBuffers* buffers, SizeType& numberOfPoints, unsigned int& causeOfStop, PathType* path, ScalarPathType* pvals, const int direction )
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
    ElementType mjump = 0.01*(args->step); // Minimum jump until we consider we reached the seeding region
    // ---------------------------------------
    // Initialize the algortihm:
    causeOfStop = 0;
    numberOfPoints = 0;
    ElementType dpath = 0.0;
    // ---------------------------------------
    // Find the leading direction of the seeding point
    // by interpolation:
    bool inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->current, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    if(inside){
        // Just use the closest direction as the reference direction:
        SizeType nneighs = ( 1 << desc->ndim );
        ElementType mweight = -1.0;
        IndexType   mpos    = 0;
        for( IndexType i=0; i<nneighs; ++i ){
            if( buffers->weights[i]>mweight ){
                mweight = buffers->weights[i];
                mpos    = i;
            }
        }
        for(unsigned int d=0; d<desc->ndim; ++d )
            buffers->direction[d] = (
                io->dirs[ d + (desc->ndim)*(buffers->neighs[mpos]) ]
            )*direction;
    }
    else{
        for(unsigned int d=0; d<desc->ndim; ++d )
            buffers->direction[d] = 0.0;
        buffers->direction[desc->ndim-1] = (ElementType)direction;
    }
    // ---------------------------------------
    // Find the new path point by point:
    do{
        // Actual implementation of the Runge-Kutta integration step:
        ElementType curvature; // The curcature in the current target
        ElementType scval;     // The interpolated scalar value in the current target
        int status = newPoint( inargs, buffers, curvature, scval );
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
            causeOfStop = 2;
            return(dpath);
        }
        // ----------------------------
        // Otherwise, the current target must be included in the path...
        for( unsigned int d=0; d<desc->ndim; ++d )
            point[d] = buffers->current[d];
        path->push_back(point);
        pvals->push_back(scval);
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
        if(numberOfPoints==1)
            mjump    = std::max( 0.01*dist, 1.0e-3 );
        // ----------------------------
        // We consider the path is all done if we have reached a point
        // for which the scalar value is small enough
        if( scval < args->threshold ) // Exit value 0
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
int newPoint( BackTracingThreader* inargs, NewPathBuffers* buffers, ElementType& curvature, ElementType& scval )
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
    interpolateScalarImage( desc->ndim, io->mask, buffers->neighs, buffers->weights, &maskval );
    if( maskval<0.5 )
        return 3;
    // ---------------------------------------
    // If we reach here, the current target is a good one. We
    // can compute its interpolated scval:
    interpolateScalarImage( desc->ndim, io->scalar, buffers->neighs, buffers->weights, &scval );
    curvature = 2*(args->mcurv); // temporary
    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    // From this point, we implement the RK4 method to find the next point
    // in the streamline:
    // k1 is just the estimate from the current position:
    interpolateVectorImage( desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->direction, buffers->RK_k1 );
    // Move h/2 towards k1 to find a first estimate of the middle point:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + 0.5*(args->step)*(buffers->RK_k1[d]);
    // k2 will be computed as the interpolated value at the first estimate
    // for the middle point:
    inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    // The management of this condition is delegated to the next step:
    if(!inside)
        return -1;
    interpolateVectorImage( desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->direction, buffers->RK_k2 );
    // Move h/2 towards k2 to find a better estimate of the middle point:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + 0.5*(args->step)*(buffers->RK_k2[d]);
    // k3 will be computed as the interpolated value at the refined estimate
    // for the middle point:
    inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    // The management of this condition is delegated to the next step:
    if(!inside)
        return -1;
    interpolateVectorImage( desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->direction, buffers->RK_k3 );
    // Move h towards k3 to find an estimate of the final point:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + (args->step)*(buffers->RK_k3[d]);
    // k4 will be computed as the interpolated value at the estimate for
    // the final point:
    inside = computeInterpolationFactors( desc->ndim, desc->fov, buffers->next, desc->xyz2ijk, buffers->ijk, buffers->neighs, buffers->weights );
    // The management of this condition is delegated to the next step:
    if(!inside)
        return -1;
    interpolateVectorImage( desc->ndim, io->dirs, buffers->neighs, buffers->weights, buffers->direction, buffers->RK_k4 );
    // If we reach here, we have the four estimates k1, k2, k3 and k4, and
    // the actual step to jump is a weighted average of them:
    ElementType w1   = 0.16666666666666666667;
    ElementType w2   = 0.33333333333333333333;
    ElementType norm = 0.0;
    // Update the reference direction:
    for( unsigned int d=0; d<desc->ndim; ++d ){
        buffers->direction[d] = w1*(buffers->RK_k1[d]) + w2*(buffers->RK_k2[d]) + w2*(buffers->RK_k3[d]) + w1*(buffers->RK_k4[d]);
        norm += buffers->direction[d] * buffers->direction[d];
    }
    norm = sqrt(norm);
    // Normalize the reference direction:
    if(norm>1.0e-12){
        for( unsigned int d=0; d<desc->ndim; ++d )
            buffers->direction[d] /= norm;
    }
    // And, finally, step towards this direction:
    for( unsigned int d=0; d<desc->ndim; ++d )
        buffers->next[d] = buffers->current[d] + (buffers->direction[d])*(args->step);
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
 * to interpolate a scalar buffer. It is assumed that out-of-mask values are filled with
 * zeros.
 */
void interpolateScalarImage( const unsigned int ndim, const BufferType buffer, const IndexBuffer neighs, const BufferType weights, BufferType interpolated )
{
    SizeType nneighs = ( 1 << ndim );
    *interpolated = 0.0;
    for( IndexType i=0; i<nneighs; ++i )
        *interpolated += weights[i]*(buffer[neighs[i]]);
    return;
}

/**
 * Linearly interpolate a vector field with the same dimensions as the image
 * (tipically, the field of eigenvectors of the diffusion tensor) preserving
 * the coherence with a given reference direction
 */
void interpolateVectorImage( const unsigned int ndim, const BufferType dirs, const IndexBuffer neighs, const BufferType weights, const BufferType refdir, BufferType iDir )
{
    // Number of neighbors to use:
    SizeType nneighs = ( 1 << ndim );
    // Initialize and compute the norm of the reference direction:
    ElementType rnorm = 0.0;
    for( unsigned int d=0; d<ndim; ++d ){
        iDir[d] = 0.0;
        rnorm  += refdir[d] * refdir[d];
    }
    rnorm = sqrt(rnorm);
    // For each neighbor...
    for( IndexType i=0; i<nneighs; ++i ){
        // Compute the dot product between the neighboring
        // direction and the reference direction, as well
        // as the norm of the former:
        ElementType cosine = 0.0;
        ElementType nnorm  = 0.0;
        for( unsigned int d=0; d<ndim; ++d ){
            ElementType dcomp = dirs[d+ndim*neighs[i]];
            cosine += refdir[d] * dcomp;
            nnorm  += dcomp * dcomp;
        }
        nnorm = sqrt(nnorm);
        if(nnorm>1.0e-12)
            cosine /= (nnorm*rnorm);
        else
            cosine = 0.0;
        // If the cosine is positive, the neighboring direction
        // is pointing in the same semiplane as the reference
        // direction. Otherwise, we have to invert it to get a
        // more coherent orientation:
        ElementType sign = 1.0;
        if(cosine<0){
            sign   = -1.0;
            cosine = -cosine;
        }
        // Linearly interpolate:
        ElementType cWeight = weights[i] * cosine;
        for( unsigned int d=0; d<ndim; ++d )
            iDir[d] += sign * cWeight * dirs[d+ndim*neighs[i]];
    }
    // Compute the norm of the interpolated vector:
    rnorm = 0.0;
    for( unsigned int d=0; d<ndim; ++d )
        rnorm   += iDir[d] * iDir[d];
    rnorm = sqrt(rnorm);
    // Finally, normalize to get a unit-norm direction:
    if(rnorm>1.0e-12){
        for( unsigned int d=0; d<ndim; ++d )
            iDir[d] /= rnorm;
    }
    return;
}

/*
void interpolateVectorImage( const unsigned int ndim, const BufferType dirs, const IndexBuffer neighs, const BufferType weights, const BufferType refdir, BufferType iDir )
{
    // Number of neighbors to use:
    SizeType nneighs = ( 1 << ndim );
    // Initialize and compute the norm of the reference direction:
    ElementType rnorm = 0.0;
    for( unsigned int d=0; d<ndim; ++d ){
        iDir[d] = 0.0;
        rnorm  += refdir[d] * refdir[d];
    }
    rnorm = sqrt(rnorm);
    // For each neighbor...
    for( IndexType i=0; i<nneighs; ++i ){
        // Compute the dot product between the neighboring
        // direction and the reference direction, as well
        // as the norm of the former:
        ElementType cosine = 0.0;
        ElementType nnorm  = 0.0;
        for( unsigned int d=0; d<ndim; ++d ){
            ElementType dcomp = dirs[d+ndim*neighs[i]];
            cosine += refdir[d] * dcomp;
            nnorm  += dcomp * dcomp;
        }
        nnorm = sqrt(nnorm);
        // If the cosine is positive, the neighboring direction
        // is pointing in the same semiplane as the reference
        // direction. Otherwise, we have to invert it to get a
        // more coherent orientation:
        ElementType sign = ( cosine>=0 ? 1.0 : -1.0 );
        // Linearly interpolate:
        for( unsigned int d=0; d<ndim; ++d )
            iDir[d] += sign * weights[i] * dirs[d+ndim*neighs[i]];
    }
    // Compute the norm of the interpolated vector:
    rnorm = 0.0;
    for( unsigned int d=0; d<ndim; ++d )
        rnorm   += iDir[d] * iDir[d];
    rnorm = sqrt(rnorm);
    // Finally, normalize to get a unit-norm direction:
    if(rnorm>1.0e-12){
        for( unsigned int d=0; d<ndim; ++d )
            iDir[d] /= rnorm;
    }
    return;
}
 */

void allocateNewPathBuffers( NewPathBuffers* buffers, FovDescription* desc, BufferType costs )
{
    const unsigned int ndim = desc->ndim;
    buffers->prev = new ElementType[ndim+1];    // Plus 1 because we use homogeneous coords
    buffers->current = new ElementType[ndim+1]; // Plus 1 because we use homogeneous coords
    buffers->next = new ElementType[ndim+1];    // Plus 1 because we use homogeneous coords
    buffers->direction = new ElementType[ndim];
    buffers->RK_k1 = new ElementType[ndim];
    buffers->RK_k2 = new ElementType[ndim];
    buffers->RK_k3 = new ElementType[ndim];
    buffers->RK_k4 = new ElementType[ndim];
    buffers->ijk = new ElementType[ndim+1];     // Plus 1 because we use homogeneous coords
    SizeType nneighs = ( 1 << ndim );           // 2^ndim
    buffers->neighs = new IndexType[nneighs];
    buffers->weights = new ElementType[nneighs];
    return;
}

void freeNewPathBuffers( NewPathBuffers* buffers )
{
    delete[] buffers->prev;
    delete[] buffers->current;
    delete[] buffers->next;
    delete[] buffers->direction;
    delete[] buffers->RK_k1;
    delete[] buffers->RK_k2;
    delete[] buffers->RK_k3;
    delete[] buffers->RK_k4;
    delete[] buffers->ijk;
    delete[] buffers->neighs;
    delete[] buffers->weights;
    return;
}

