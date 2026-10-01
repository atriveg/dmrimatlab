/*==========================================================
 * atti2dkt_.c
 *
 * This is a core function to atti2dkt, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/sh2hot.h"
#include "../quadprog/dmriquadprog.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

const unsigned int NK = 15;
const unsigned int ND = 6;
const unsigned int NS = 1;
const unsigned int NR = NK+ND+NS;

typedef struct DTIIOData{
    // Input:
    BufferType signal;
    BufferType gi;
    BufferType bi;
    // Output:
    BufferType dkt;
    BufferType dti;
    BufferType qps;
    BufferType S0;
} DTIIOData;

typedef struct DTIParameters{
    char estimS0; // 'y'/'n'
    char mode;    // 'o'/'w'/'q'
    unsigned int wlsit;
    ElementType wsc;
    ElementType rcondth;
    SizeType nconst;
    char negkrt;  // 'y'/'n'
    BufferType At;
    BufferType Drec;
} DTIParameters;

class ThArgs : public DMRIThreader
{
public:
    SizeType G;          // The number of gradients
    DTIIOData* io;
    DTIParameters* params;
};


THFCNRET atti2dkt_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  signal, the signal to fit, G x N
     * prhs[1]:  gi, the gradients table, G x 3
     * prhs[2]:  bi, the q-values vector, G x 1
     * prhs[3]:  opts, a structure with the parameters to the algorithm
     *       opts.wlsit, the number of iterations with WLS, 1 x 1
     *       opts.wsc: the minimum weight to be applied in the WLS problem with
     *                 respect to the maximum one, 1 x 1
     *       opts.rcondth, for WLS, the minimum condition number so that the WLS matrix can be inverted, 1 x 1
     *       opts.mode, the algorithm to be used:
     *            'o': ordinary LS (log domain)
     *            'w', weigthed LS (log domain)
     *            'q', quadratic programming (log domain)
     *       opts.At, 21 x nconst the set of inequality constraints if a QP problem is solved (transposed)
     *       opts.Drec, nconst/3 x 6 the reconstruction matrix for the tensor model at the constraint directions
     *       opts.negkrt, for QP, whether ('y') or not ('n') allow kurtoses values in [-2,0], 1 x 1
     *       opts.estimS0, whether or not estimate the correction to S0.
     * prhs[4], maxthreads, the maximum number of threads in POSIX systems, 1 x 1
     *
     *  OUTPUTS:
     *
     * plhs[0]: dkt, 15 x N, the unique components of the Kurtosis tensor
     * plhs[1]: dti, 6 x N, the unique components of the diffusion tensor
     * plhs[2]: S0, 1 x N, the (normalized) baseline estimated, will be near 1 in all cases. Note
     *    it is necessarily normalize since we pass S_i/S_0 to the function instead of S_i AND S_0
     * plhs[3]: qps, 1 x N, info about the final status of the quadratic program
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:atti2dkt_:callstack","This function should only be called from atti2dkt");
    else if( strcmp(callerFunc,"atti2dkt") )
        mexErrMsgIdAndTxt("MyToolbox:atti2dkt_:callstack","This function should only be called from atti2dkt");
    //=======================================================================================
    if(nrhs!=5)
        mexErrMsgIdAndTxt("MyToolbox:atti2dkt_:nrhs","Exactly 5 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType G = mxGetM(prhs[0]); // The number of gradient directions per voxel
    //=======================================================================================
    DTIParameters params;
    // --
    params.estimS0 = 'n';
    mxArray* mxestimS0 = mxGetField( prhs[3], 0, "estimS0");
    if( mxestimS0!=NULL )
        params.estimS0 = (char)(mxGetChars(mxestimS0)[0]);
    // --
    params.mode = 'q';
    mxArray* mxmode = mxGetField( prhs[3], 0, "mode");
    if( mxmode!=NULL ){ params.mode = (char)(mxGetChars(mxmode)[0]); }
    // --
    params.wlsit = 5;
    mxArray* mxwlsit = mxGetField( prhs[3], 0, "wlsit");
    if( mxwlsit!=NULL ){ params.wlsit = (unsigned int)mxGetScalar(mxwlsit); }
    // --
    params.wsc = 0.01;
    mxArray* mxwsc = mxGetField( prhs[3], 0, "wsc");
    if( mxwsc!=NULL ){ params.wsc = mxGetScalar(mxwsc); }
    // --
    params.rcondth = 1.0e-6;
    mxArray* mxrcondth = mxGetField( prhs[3], 0, "rcondth");
    if( mxrcondth!=NULL ){ params.rcondth = mxGetScalar(mxrcondth); }
    // --
    params.At = NULL;
    params.nconst = 0;
    mxArray* mxAt = mxGetField( prhs[3], 0, "At");
    if( mxAt!=NULL ){
        params.At = (BufferType)mxGetDoubles(mxAt);
        params.nconst = mxGetN(mxAt);
    }
    // --
    params.Drec = NULL;
    mxArray* mxDrec = mxGetField( prhs[3], 0, "Drec");
    if( mxDrec!=NULL )
        params.Drec = (BufferType)mxGetDoubles(mxDrec);
    // --
    params.negkrt = 'n';
    mxArray* mxnegkrt = mxGetField( prhs[3], 0, "negkrt");
    if( mxnegkrt!=NULL )
        params.negkrt = (char)(mxGetChars(mxnegkrt)[0]);

    //=======================================================================================
    DTIIOData io;
    // ------ Inputs
    io.signal = mxGetDoubles(prhs[0]);
    io.gi = mxGetDoubles(prhs[1]);
    io.bi = mxGetDoubles(prhs[2]);
    // ------ Outputs
    if(nlhs!=4)
        mexErrMsgIdAndTxt("MyToolbox:atti2dti_:nlhs","This function accepts just 4 output arguments");
    plhs[0] = mxCreateDoubleMatrix( 15, N, mxREAL );
    io.dkt = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 6, N, mxREAL );
    io.dti = mxGetDoubles(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix( 1, N, mxREAL );
    io.S0 = mxGetDoubles(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix( 1, N, mxREAL );
    io.qps = mxGetDoubles(plhs[3]);
    mataux::setValueMxArray( io.qps, 1, N, 0.0f );
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[4]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, 20 );
    // Own values:
    threader.G      = G;
    threader.io     = &io;
    threader.params = &params;
    //=======================================================================================
    threader.threadedProcess( maxthreads, atti2dkt_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET atti2dkt_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    DTIIOData* io = args->io;
    DTIParameters* params = args->params;

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads_thread(1);
    
    // Convenience constants:
    SizeType G = args->G;
    unsigned int NR = ( params->estimS0=='y' ? NK+ND+1 : NK+ND );

    // Allocate auxiliar buffers for computations
    // ---------------------------------------------------------------
    ElementType dkt[NK];
    ElementType dti[ND];
    BufferType x  = new ElementType[NR];
    BufferType x2 = new ElementType[NR];
    ElementType eigval[3];
    ElementType eigvec[9];
    const BLAS_INT dim = 3;
    BLAS_INT info = 0;
    ElementType work[9]; // According to Lapack's docs for dspev 
    // ---------------------------------------------------------------
    // The gradients matrix:
    BufferType lSi  = new ElementType[G];
    BufferType lSiw = new ElementType[G];
    BufferType lSii = new ElementType[NR];
    BufferType A    = new ElementType[G*NR];
    BufferType A2   = new ElementType[G*NR];
    mataux::setValueMxArray( A, G, NR, 0.0f );
    // ---------------------------------------------------------------
    ElementType scale2 = 0.0;
    // ------------
    // Populate the columns of A related to the Kurtosis tensor (0-14)
    unsigned int* nx = new unsigned int[NK];
    unsigned int* ny = new unsigned int[NK];
    unsigned int* nz = new unsigned int[NK];
    sh2hot::computeHOTPowers(4,nx,ny,nz);
    unsigned long* mu = new unsigned long[NK];
    sh2hot::computeHOTMultiplicity(4,mu,nx,ny,nz);
    // ---
    for( IndexType g=0; g<(IndexType)G; ++g ){
        for( unsigned int c=0; c<NK; ++c ){
            ElementType gix = io->gi[g];
            ElementType giy = io->gi[g+G];
            ElementType giz = io->gi[g+2*G];
            A[G*c+g]  = mu[c] * pow(gix,nx[c]) * pow(giy,ny[c]) * pow(giz,nz[c]);
            A[G*c+g] *= (io->bi[g])*(io->bi[g])/6;
            // ---
            scale2 += (A[G*c+g])*(A[G*c+g]);
        }
    }
    delete[] nx;
    delete[] ny;
    delete[] nz;
    delete[] mu;
    // ------------
    // Populate the columns of A related to the DTI tensor (15-20):
    ElementType scale1 = 1000.0; // Since the Kurtosis encoding is multiplied by b^2
    for( IndexType g=0; g<(IndexType)G; ++g ){
        A[G*(NK+0)+g] = -(io->bi[g]) * (io->gi[g]) * (io->gi[g]) * scale1;
        A[G*(NK+1)+g] = -2.0f * (io->bi[g]) * (io->gi[g]) * (io->gi[G+g]) * scale1;
        A[G*(NK+2)+g] = -2.0f * (io->bi[g]) * (io->gi[g]) * (io->gi[G*2+g]) * scale1;
        A[G*(NK+3)+g] = -(io->bi[g]) * (io->gi[G+g]) * (io->gi[G+g]) * scale1;
        A[G*(NK+4)+g] = -2.0f * (io->bi[g]) * (io->gi[G+g]) * (io->gi[G*2+g]) * scale1;
        A[G*(NK+5)+g] = -(io->bi[g]) * (io->gi[G*2+g]) * (io->gi[G*2+g]) * scale1;
        // ---
        for( unsigned int c=0; c<ND; ++c )
            scale2 += (A[G*(NK+c)+g])*(A[G*(NK+c)+g]);
    }
    // ------------
    // Populate the last column of A, related to S0, if necessary:
    if( params->estimS0=='y' ){
        scale2 = sqrt( scale2 / (G*(ND+NK)) );
        for( IndexType g=0; g<(IndexType)G; ++g )
            A[G*(NK+ND)+g] = scale2;
    }
    // ---------------------------------------------------------------
    // Compute the fixed fitting matrix for ordinary LS:
    BufferType ATA = new ElementType[NR*NR];
    mataux::transposeMultiplyMxArray( A, (BufferType)ATA, G, NR );
    // ---
    IndexBuffer pivot1 = new IndexType[NR];
    IndexBuffer pivot2 = new IndexType[NR];
    BLAS_INT flag = 1;
    BLAS_INT K = NR;
    BLAS_INT lname = 6; // Length of "dgetri"
    BLAS_INT largs = 0; // Length of ""
#ifdef OCTAVE_BUILD
    // ilaenv seems not to be present in liblapack. Use BS=4 as a
    // "one size fits all" thing:
    BLAS_INT BS = 4;
#else
    BLAS_INT BS  = LAPACKCALLFCN(ilaenv)(
            &flag, "dgetri", "", 
            &K, &K, &K, &K, lname, largs );
    BS = ( BS>4 ? BS : 4 );
#endif
    SizeType lwork = BS*NR;
    BufferType work2 = new ElementType[lwork];
    mataux::checkAndInvertMxArray( (BufferType)ATA, NR, mxGetEps()*10, (IndexBuffer)pivot1, (IndexBuffer)pivot2, lwork, work2 );
    // ---
    BufferType AT  = new ElementType[NR*G];
    mataux::transposeMxArray( A, AT, G, NR );
    // ---
    BufferType piA = new ElementType[NR*G];
    mataux::multiplyMxArrays( ATA, AT, piA, NR, NR, G );
    // ---------------------------------------------------------------
    BufferType wi = new ElementType[G];
    // ---------------------------------------------------------------
    // If QP has to be solved:
    // Data structures for quadratic programming
    dmriqpp::QPProblem  qpproblem; // The problem itself
    dmriqpp::QPAuxiliar qpaux;     // Auxiliar buffers
    dmriqpp::QPParams   qpparams;  // Parameters for the algorithm
    // Initallize these structures. There are:
    //    - NR = NK+ND+NS variables to optimize,
    //    - CE = 0 equality constraints
    //    - CI equality constraints
    //    - CL = 0 lower bounds
    //    - CU = 0 upper bounds
    const unsigned int CE = 0;
    const unsigned int CI = params->nconst;
    const unsigned int CL = 0;
    const unsigned int CU = 0;
    dmriqpp::allocateQPProblem( qpproblem, NR, CE, CI, CL, CU );
    // Finally, set the parameters to the algorithm:
    qpproblem.step = 1.0; // Will be overwritten since qpparams.computes0 = true;
    qpparams.maxiters = 100000;
    qpparams.streak = 3;
    qpparams.steptol = 1.0e-18;
    qpparams.costtol = 1.0e-18;
    qpparams.normalize = true;
    qpparams.computes0 = true;
    // Now we can allocate all auxiliar buffers:
    dmriqpp::allocateQPAuxiliar( qpproblem,  qpparams, qpaux );
    // The inequality constraints are constant for all voxels:
    if(CI>0){
        mataux::setValueMxArray( qpproblem.b, CI, 1, 0.0 );
        memcpy( qpproblem.At, params->At, CI*NR*sizeof(ElementType) );
    }
    // -----------------------------------------------------------------------------
    // ---
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for(  IndexType i=start; i<end; ++i ){
            //---------------------------------------------------------------------
            // Regardless of the working mode, we will always begin with an OLS
            // approximation to the solution, which reduces to the product of the
            // pre-computed pseudo-inverse of A, piA, with the acquired signal:
            memcpy( lSi, &(io->signal[i*G]), G*sizeof(ElementType) );
            for( IndexType g=0; g<(IndexType)G; ++g )
                lSi[g] = log( lSi[g]>mxGetEps() ? lSi[g] : mxGetEps() );
            mataux::multiplyMxArrays( piA, lSi, (BufferType)x, NR, G, 1 );
            //---------------------------------------------------------------------
            IndexType result = 0;
            if( params->mode=='w' || params->mode=='q' ){
                mataux::multiplyMxArrays( AT, lSi, lSii, NR, G, 1 );
                // In any of these cases we need to iterate to compute the
                // WLS solution
                for( unsigned int n=0; n<params->wlsit; ++n ){
                    // 1- Compute the weights wi at this iteration:
                    mataux::multiplyMxArrays( A, (BufferType)x, wi, G, NR, 1 );
                    ElementType maxw = -1.0f;
                    bool anynan = false;
                    for( IndexType g=0; g<(IndexType)G; ++g ){
                        wi[g] = exp(2.0f*wi[g]);
                        wi[g] = ( wi[g]>1.0 ? 1.0 : wi[g] );
                        wi[g] = ( wi[g]<0.0 ? 0.0 : wi[g] );
                        maxw = ( wi[g]>maxw ? wi[g] : maxw );
                        anynan |= isnan(wi[g]);
                    }
                    if(anynan)
                        break;
                    for( IndexType g=0; g<(IndexType)G; ++g ){
                        if( wi[g] < maxw*(params->wsc) )
                            wi[g] = maxw*(params->wsc);
                    }
                    // 2- Create the WLS problem:
                    memcpy( A2, A, G*NR*sizeof(ElementType) );
                    for( IndexType g=0; g<(IndexType)G; ++g ){
                        lSiw[g] = lSi[g] * wi[g];
                        for( unsigned int c=0; c<NR; ++c )
                            A2[c*G+g] *= wi[g];
                    }
                    mataux::multiplyMxArrays( AT, A2, ATA, NR, G, NR );
                    mataux::multiplyMxArrays( AT, lSiw, lSii, NR, G, 1 );
                    // 3- Solve the WLS problem:
                    result = mataux::checkAndInvertMxArray( (BufferType)ATA, NR, params->rcondth,
                                                                     (IndexBuffer)pivot1, (IndexBuffer)pivot2, lwork, work2 );
                    if( result == 0 ){
                        mataux::multiplyMxArrays( (BufferType)ATA, lSii, (BufferType)x2, NR, NR, 1 );
                        for( unsigned int c=0; c<NR; ++c )
                            anynan |= isnan(x2[c]);
                        if(anynan){
                            result = 1;
                            break;
                        }
                        else
                            memcpy( x, x2, NR*sizeof(ElementType) );
                    }
                    else
                        break;
                }
            }
            //---------------------------------------------------------------------
            if( (params->mode=='q') && (info==0) ){
                // NOTE: if info!=0, the eigenvalue computation failed and
                // we cannot assure there is a feasible first iteration
                // Otherwise, we will find the solution by solving a quadratic program
                // Note this makes sense only in case the variable "result" indicates
                // the algorithm succeeded to invert A^T*W*A
                io->qps[i] = -4;
                if(result==0){
                    io->qps[i] = -2;
                    memcpy( qpproblem.Qi, ATA, NR*NR*sizeof(ElementType) );
                    memcpy( qpproblem.f, lSii, NR*sizeof(ElementType) );
                    if( (CI>0) && (params->negkrt=='y') ){
                        // Prepare the independent term of the inequality constraints:
                        const unsigned int CIPD = CI/3;
                        mataux::multiplyMxArrays( params->Drec, &x[NK], &qpproblem.b[CIPD], CIPD, ND, 1 );
                        for( unsigned int d=CIPD; d<2*CIPD; ++d )
                            qpproblem.b[d] = 2.0*(qpproblem.b[d])*(qpproblem.b[d]);
                    }
                    mataux::scalaropMxArray( qpproblem.f, NR, 1, -1.0f, mataux::MULTIPLY );
                    int qp_result = dmriqpp::solveQuadraticProgram( qpproblem,
                                                                    qpaux, qpparams );
                    if(qp_result>=0){
                        memcpy( x, qpproblem.x, NR*sizeof(ElementType) );
                        io->qps[i] = qp_result;
                    }
                }
            }
            //---------------------------------------------------------------------
            // Re-scale the solutions vector:
            for( unsigned int c=NK; c<NK+ND; ++c )
                x[c] *= scale1;
            ElementType md = (x[NK]+x[NK+3]+x[NK+5])/3;
            for( unsigned int c=0; c<NK; ++c )
                x[c] /= (md*md);
            // At this point, we can copy the outputs to the proper place:
            memcpy( &(io->dkt[NK*i]), x, NK*sizeof(ElementType) );
            memcpy( &(io->dti[ND*i]), &(x[NK]), ND*sizeof(ElementType) );
            // If S0 has to be estimated:
            if( params->estimS0=='y' )
                io->S0[i] = exp(x[NK+ND]*scale2);
            else
                io->S0[i] = 1.0;
        }
    }
    while( start < args->getN() );
    
    blas_num_threads_thread(blas_threads);
    
    // Free memory previously allocated
    delete[] x;
    delete[] x2;
    delete[] lSi;
    delete[] lSiw;
    delete[] lSii;
    delete[] A;
    delete[] A2;
    delete[] ATA;
    delete[] pivot1;
    delete[] pivot2;
    delete[] work2;
    delete[] AT;
    delete[] piA;
    delete[] wi;

    return (THFCNRET)NULL;
}
