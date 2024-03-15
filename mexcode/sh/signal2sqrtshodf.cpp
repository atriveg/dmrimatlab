/*==========================================================
 * signal2sqrtshodf.c
 *
 * This function fits the squared root of a strictly 
 * positive ODF with unit mass to a set of measurements
 * according to a convolution model
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
#include "../mathsmex/posODFsMaths.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

typedef struct IOData{
    BufferType  E;
    BufferType  lambda;
    BufferType  pshell;
    BufferType  psi0;
    BufferType  psi;
    BufferType  nit;
    BufferType  mu;
    BufferType  Q;
    BufferType  gradnorm;
    BufferType  conderr;
} IOData;

class ThArgs : public DMRIThreader
{
public:
    SizeType                         M;      // The number of shells in the gradients scheme
    char                             algorithm;
    posODFs::ProblemFeatures*        features;
    posODFs::WignerSymbols*          wigner;
    posODFs::ODFAlgorithmParameters* parameters;
    IOData*                          io;
};


THFCNRET posodfsh_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]: E, the signal to fit, N x FOV
     * prhs[1]: lambda, the convolution factors (L+1) x M x FOV
     * prhs[2]: pshell, a poiter telling which of the M shells each sample belongs to, N x 1
     * prhs[3]: psi0, the initial guess for the (squared root of the) ODF, K x FOV, with K=(L+1)(L+2)/2
     * prhs[4]: Ylm, the SH matrix for gi and L, N x Kp, with Kp=(2L+1)(2L+2)/2
     * prhs[5]: opts, a strucutre with optional parameters:
     * 
     *       opts.nu:         the Laplace-Beltrami penalty (over the ODF itself)
     *       opts.algorithm:  'L' for Levenberg-Marquardt's, 'N' for Newton-Raphson's, or 'C' for combined
     *       opts.T:          maximum number of iterations used with either method
     *       opts.maxfails:   maximum number of consecutive failed iterations
     *       opts.thCost:     convergence threshold for the cost function (only with Lebenberg-Marquardt's)
     *       opts.thGrad:     convergence threshold for the modulus of the gradient (only with Newton-Raphson's)
     *       opts.thCond:     convergence threshold for the unit-norm constraint (only with Newton-Raphson's)
     *       opts.rho0:       the initial value of the damping factor for Levenberg-Marquardt's algorithm
     *       opts.minrcn:     minimum reciprocal condition number of a matrix before it is considered nearly-singular
     *       opts.psi0eps:    minimum allowed value for psi[0]^2 (only with Levenberg-Marquardt's)
     *       opts.maxthreads: the maximum number of threads allowed
     *
     *  OUTPUTS:
     *
     * plhs[0]: psi, the SH coefficients of the squared root of the ODF, K x FOV, with K=(L+1)(L+2)/2
     * plhs[1]: nit, the number of iterations used to fit the model, 1 x FOV
     * plhs[2]: mu, the Lagrange multiplier at each voxel if Newton-Raphson's method is used, 1 x FOV
     * plhs[3]: Q, the final cost/Lagrangian value, 1 x FOV
     * plhs[4]: gradnorm, the final norm of the gradient, 1 x FOV
     * plhs[5]: conderr, the final value of the error in the constraint, 1 x FOV
     * 
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:callstack","This function should only be called from micro2shsqrtodf");
    else if( strcmp(callerFunc,"micro2shsqrtodf") && strcmp(callerFunc,"test_signal2sqrtshodf") )
        mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:callstack","This function should only be called from micro2shsqrtodf");
    //=======================================================================================
    if(nrhs<4)
        mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:nrhs","At least 4 inptus are required");
    //=======================================================================================
    // Process the inputs, retrieve useful constants:
    SizeType     N   = mxGetM(prhs[0]); // Number of dMRI meaasurements per voxel
    SizeType     M;
    SizeType     FOV = mxGetN(prhs[0]); // Number of voxels to be processed
    unsigned int L   = (unsigned int)(mxGetDimensions(prhs[1])[0]-1); // Maximum SH order of psi (squared root of the ODF)
    unsigned int K   = (L+1)*(L+2)/2; // Number of SH coefficients of psi (squared root of the ODF)
    unsigned int Kp  = (2*L+1)*(2*L+2)/2; // Number of SH coefficients of phi (the ODF)
    //=======================================================================================
    // Sanity checks on input parameters
    // -- lambda
    if( FOV>1 ){
        if( mxGetNumberOfDimensions(prhs[1])<3 )
            mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:lambda","lambda should have size (L+1) x N x FOV");
        M = mxGetDimensions(prhs[1])[1]; // Number of dmRI shells
        if( mxGetDimensions(prhs[1])[2]!=FOV )
            mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:lambda","lambda should have size (L+1) x M x FOV");
    }
    else
        M = mxGetN(prhs[1]);
    // -- pshell
    if( mxGetM(prhs[2])!=N )
        mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:pshell","pshell should have size N x 1");
    for( IndexType n=0; n<(IndexType)N; ++n ){
        IndexType ps = (IndexType)(mxGetDoubles(prhs[2])[n]);
        if( (ps<0) || (ps>M-1) )
            mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:pshell","pshell value out of bounds");
    }
    // -- psi
    if( mxGetM(prhs[3])!=K || mxGetN(prhs[3])!=FOV )
        mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:psi","psi should have size K x FOV");
    // -- Ylm
    if( mxGetM(prhs[4])!=N || mxGetN(prhs[4])!=Kp )
        mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:Ylm","Ylm should have size N x Kp");
    //=======================================================================================
    // Retrieve the appropraite options from the options structure
    // ---
    ElementType nu = 0.001;
    mxArray* mxnu = mxGetField( prhs[5], 0, "nu");
    if( mxnu!=NULL ){ nu = mxGetScalar(mxnu); }
    // ---
    char algorithm[2] = "C";
    mxArray* mxalgorithm = mxGetField( prhs[5], 0, "algorithm");
    if( mxalgorithm!=NULL ){ mxGetString(mxalgorithm,algorithm,2); }
    switch(algorithm[0]){
        case 'L':
        case 'N':
        case 'C':
            break;
        default:
            mexErrMsgIdAndTxt("MyToolbox:signal2sqrtshodf:algorithm","algorithm must be one of 'L', 'N' or 'C'");
    }
    // ---
    unsigned int T = 100;
    mxArray* mxT = mxGetField( prhs[5], 0, "T");
    if( mxT!=NULL ){ T = (unsigned int)mxGetScalar(mxT); }
    // ---
    ElementType maxfails = 10;
    mxArray* mxmaxfails = mxGetField( prhs[5], 0, "maxfails");
    if( mxmaxfails!=NULL ){ maxfails = mxGetScalar(mxmaxfails); }
    // ---
    ElementType thCost = 0.0001;
    mxArray* mxthCost = mxGetField( prhs[5], 0, "thCost");
    if( mxthCost!=NULL ){ thCost = mxGetScalar(mxthCost); }
    // ---
    ElementType thGrad = 0.0001;
    mxArray* mxthGrad = mxGetField( prhs[5], 0, "thGrad");
    if( mxthGrad!=NULL ){ thGrad = mxGetScalar(mxthGrad); }
    // ---
    ElementType thCond = 0.0001;
    mxArray* mxthCond = mxGetField( prhs[5], 0, "thCond");
    if( mxthCond!=NULL ){ thCond = mxGetScalar(mxthCond); }
    // ---
    ElementType rho0 = 1.0;
    mxArray* mxrho0 = mxGetField( prhs[5], 0, "rho0");
    if( mxrho0!=NULL ){ rho0 = mxGetScalar(mxrho0); }
    // ---
    ElementType minrcn = 1.0e-6;
    mxArray* mxminrcn = mxGetField( prhs[5], 0, "minrcn");
    if( mxminrcn!=NULL ){ minrcn = mxGetScalar(mxminrcn); }
    // ---
    ElementType psi0eps = 0.01*0.01;
    mxArray* mxpsi0eps = mxGetField( prhs[5], 0, "psi0eps");
    if( mxpsi0eps!=NULL ){ psi0eps = mxGetScalar(mxpsi0eps); }
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( 1e6 );
    mxArray* mxmaxthreads = mxGetField( prhs[5], 0, "maxthreads");
    if( mxmaxthreads!=NULL )
        maxthreads = ( (unsigned int)mxGetScalar(mxmaxthreads)<maxthreads ? (unsigned int)mxGetScalar(mxmaxthreads) : maxthreads );
    //=======================================================================================
    // Set the problem features:
    posODFs::ProblemFeatures features;
    features.N = N;
    features.L = L;
    features.nu = nu;
    features.Y = mxGetDoubles(prhs[4]);
    //=======================================================================================
    // Compute Wigner's symbols
    posODFs::WignerSymbols wigner;
    posODFs::createWignerSymbols( &wigner, L );
    //=======================================================================================
    // Set the parameters of the algorithm
    posODFs::ODFAlgorithmParameters parameters;
    parameters.T = T;
    parameters.maxfails = (unsigned int)maxfails;
    parameters.thCost = thCost;
    parameters.thGrad = thGrad;
    parameters.thCond = thCond;
    parameters.rho0 = rho0;
    parameters.minrcn = minrcn;
    parameters.psi0eps = psi0eps;
    //=======================================================================================
    // Create the I/O  parameters
    IOData io;
    // Three mandatory inputs:
    io.E      = mxGetDoubles(prhs[0]);
    io.lambda = mxGetDoubles(prhs[1]);
    io.pshell = mxGetDoubles(prhs[2]);
    io.psi0   = mxGetDoubles(prhs[3]);
    if(nlhs<1)
        return;
    // All outputs are optional:
    //----- The SH coefficients of the squared root of the ODF
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( K, FOV, mxREAL ); // psi
        io.psi = mxGetDoubles(plhs[0]);
    }
    else
        return;
    //----- The number of iterations it took to converge
    if(nlhs>1){
        plhs[1] = mxCreateDoubleMatrix( 1, FOV, mxREAL ); // nit
        io.nit = mxGetDoubles(plhs[1]);
    }
    else
        io.nit = (BufferType)NULL;
    //----- The Lagrange multiplier
    if(nlhs>2){
        plhs[2] = mxCreateDoubleMatrix( 1, FOV, mxREAL ); // mu
        io.mu = mxGetDoubles(plhs[2]);
    }
    else
        io.mu = (BufferType)NULL;
    //----- The final cost/Lagrangian
    if(nlhs>3){
        plhs[3] = mxCreateDoubleMatrix( 1, FOV, mxREAL ); // Q
        io.Q = mxGetDoubles(plhs[3]);
    }
    else
        io.Q = (BufferType)NULL;
    //----- The final value of the gradient norm
    if(nlhs>4){
        plhs[4] = mxCreateDoubleMatrix( 1, FOV, mxREAL ); // gradnorm
        io.gradnorm = mxGetDoubles(plhs[4]);
    }
    else
        io.gradnorm = (BufferType)NULL;
    //----- The final value of the error in the fulfillment of the constraint
    if(nlhs>5){
        plhs[5] = mxCreateDoubleMatrix( 1, FOV, mxREAL ); // conderr
        io.conderr = mxGetDoubles(plhs[5]);
    }
    else
        io.conderr = (BufferType)NULL;

    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( FOV, 20 );
    // Own values:
    threader.M          = M;
    threader.algorithm  = algorithm[0];
    threader.features   = &features;
    threader.wigner     = &wigner;
    threader.parameters = &parameters;
    threader.io         = &io;
    //=======================================================================================
    threader.threadedProcess( maxthreads, posodfsh_process_fcn );
    //=======================================================================================
    
    posODFs::destroyWignerSymbols( &wigner );
    
    return;
}

THFCNRET posodfsh_process_fcn( void* inargs )
{
    // Retrive the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads(1);
    
    // Convenience constants:
    SizeType N = args->features->N;
    SizeType M = args->M;
    unsigned int L = args->features->L;
    unsigned int K   = (L+1)*(L+2)/2;
    // Declare cost/gradient values to be (optionally) returned:
    ElementType Q = 0.0; // The final cost/Lagrangian
    ElementType gradnorm = 0.0; // The final norm of the gradient
    // If we use Newton-Raphson's, the error in the constraint
    // will be overwritten; if we use Levenberg-Marquardt's,
    // it is null by construction:
    ElementType conderr = 0.0; // The final mismatch in the unit norm constraint
    // Create and allocate workspaces for the algorithms:
    posODFs::VoxelFeatures voxel;
    posODFs::NRWorkBuffers nrwork;
    posODFs::LMWorkBuffers lmwork;
    posODFs::createVoxelFeatures( &voxel, N, M, L );
    // Populate the shells pointer:
    for( IndexType n=0; n<(IndexType)N; ++n )
        voxel.pshell[n] = (IndexType)(args->io->pshell[n]);
    // If we use Newton-Raphson's, the Lagrange multiplier
    // will be overwritten; if we use Levenberg-Marquardt's,
    // it is not defined:
    voxel.mu = NAN;
    if( (args->algorithm=='L') || (args->algorithm=='C') )
        posODFs::createLMWorkBuffers( &lmwork, N, L );
    if( (args->algorithm=='N') || (args->algorithm=='C') )
        posODFs::createNRWorkBuffers( &nrwork, N, L );
    // ---------------------------------------------------------------
    int niters; // The number of iterations it took to converge
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( IndexType i=start; i<end; ++i ){
            // First, copy voxel-dependent values to the appropriate buffers:
            memcpy(
                voxel.E,
                &(args->io->E[i*N]),
                N*sizeof(ElementType)
            );
            memcpy(
                voxel.lambda,
                &(args->io->lambda[i*(L+1)*M]),
                (L+1)*M*sizeof(ElementType)
            );
            memcpy(
                voxel.psi,
                &(args->io->psi0[i*K]),
                K*sizeof(ElementType)
            );
            // Now we can call the appropriate optimization routine/routines
            if( (args->algorithm=='L') || (args->algorithm=='C') ){
                niters = posODFs::LMFitPositiveODF(
                    args->features,
                    args->wigner,
                    &voxel,
                    args->parameters,
                    &lmwork,
                    &Q,
                    &gradnorm
                );
                // No matter if the algorithm converges or not, we can securely copy
                // its result to the output buffer, since Levenberg-Marquardt's will
                // always check if each iteration improves the cost function, or will
                // be reverted otherwise:
                memcpy(
                    &(args->io->psi[i*K]),
                    voxel.psi,
                    K*sizeof(ElementType)
                );
            }
            if( (args->algorithm=='N') || (args->algorithm=='C') ){
                // Whether the algorithm is 'N' or 'C', voxel.psi contains
                // a proper initial iteration (the output of L-M's algorithm
                // or the one passed by the calling function
                niters = posODFs::NRFitPositiveODF(
                    args->features,
                    args->wigner,
                    &voxel,
                    args->parameters,
                    &nrwork,
                    &Q,
                    &gradnorm,
                    &conderr
                );
                // In this case we can copy the output back only if the iterations
                // converged, since Newton-Raphson's algorithm will not check if the
                // iterations succeeded or not:
                if(niters>=0){
                    memcpy(
                        &(args->io->psi[i*K]),
                        voxel.psi,
                        K*sizeof(ElementType)
                    );
                }
            }
            // If necessary, return the number of iterations it took to converge:
            if( args->io->nit != (BufferType)NULL )
                args->io->nit[i] = niters;
            // If necessary, return the Lagrange multiplier
            if( args->io->mu != (BufferType)NULL ){
                if(niters>=0)
                    args->io->mu[i] = voxel.mu;
                else
                    args->io->mu[i] = NAN;
            }
            // If necessary, return the final cost
            if( args->io->Q != (BufferType)NULL )
                args->io->Q[i] = Q;
            // If necessary, return the final value of the gradient
            if( args->io->gradnorm != (BufferType)NULL )
                args->io->gradnorm[i] = gradnorm;
            // If necessary, return the final mismatch in the unit norm constraint
            if( args->io->conderr != (BufferType)NULL )
                args->io->conderr[i] = conderr;
        }

    }
    while( start < args->getN() );

    posODFs::destroyVoxelFeatures( &voxel );
    if( (args->algorithm=='L') || (args->algorithm=='C') )
        posODFs::destroyLMWorkBuffers( &lmwork );
    if( (args->algorithm=='N') || (args->algorithm=='C') )
        posODFs::destroyNRWorkBuffers( &nrwork );

    blas_num_threads(blas_threads);
    
    return (THFCNRET)NULL;
}
