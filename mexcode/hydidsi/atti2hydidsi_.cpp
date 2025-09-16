/*==========================================================
 * atti2hydidsi_.c
 *
 * This is a core function to atti2hydidsi, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include <cmath>
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/sphericalHarmonics.h" // Just for the definition of PI
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"
#include "../mathsmex/sanityCheckDTI.h"
#include "../gcv/compute_gcv.h"
#include "../quadprog/dmriquadprog.h"

typedef struct HYDIIOData{
    // Input:
    BufferType atti;
    BufferType dti;
    BufferType gi;
    BufferType qi;
    BufferType DFT;
    BufferType u;
    BufferType v;
    BufferType w;
    // Output:
    BufferType eap;
    BufferType Qx;
    BufferType Qy;
    BufferType Qz;
    BufferType resn;
    BufferType lapln;
    BufferType lopt;
} HYDIIOData;

typedef struct HYDIParameters{
    ElementType lambda;
    ElementType ADC0;
    ElementType lRth;
    ElementType tau;
    unsigned int miters;
    ElementType otol;
    ElementType stol;
    ElementType ctol;
} HYDIParameters;

class ThArgs : public DMRIThreader
{
public:
    SizeType G;          // The number of gradients
    SizeType K;          // The number of EAP points
    SizeType dftdim[2];  // The size of the DFT matrix
    SizeType lattice[3]; // The lattice size
    HYDIIOData* io;
    HYDIParameters* params;
};


THFCNRET atti2hydidsi_process_fcn( void* );

void setDTISolution( const BufferType rx, const BufferType ry, const BufferType rz,
                    const ElementType l1, const ElementType l2, const ElementType l3,
                    const ElementType tau, const SizeType NR, BufferType eap );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     *
     * prhs[0]:  atti, the attenuation signal to fit, G x N
     * prhs[1]:  dti, the tensor fit of that signal, 6 x N
     * prhs[2]:  gi, the gradients table, G x 3
     * prhs[3]:  qi, the q-values vector, G x 1
     * prhs[4]:  lattice, the lattice size, 3 x 1
     * prhs[5]:  lambda, the regularization parameter, 1 x 1
     * prhs[6]:  DFT, the DFT matrix to compute the Laplacian penalty, dftdim[0] x dftdim[1]
     * prhs[7]:  uvw, the locations to interpret DFT, dftdim[0] x 3
     * prhs[8]:  ADC0, the free-water diffusivity at body temperature, 1 x 1
     * prhs[9],  lRth, the cropping factor for the tails of the EAP, 1 x 1
     * prhs[10], tau, the effective diffusion time, 1 x 1
     * prhs[11], optim, a structure with the parameters to the optimization algorithm
     *           (see the help on dmriQuadprog)
     *       optim.miters, maximum number of iterations
     *       optim.otol, tolerance for the optimality condition, 1 x 1
     *       optim.stol, tolerance for the optimization step, 1 x 1
     *       optim.ctol, tolerance for the constraints, 1 x 1
     * prhs[12], maxthreads, the maximum number of threads in POSIX systems, 1 x 1
     *
     *  OUTPUTS:
     *
     * plhs[0]: eap, (prod(2*lattice+1)+1)/2 x N, the EAP computed at lattice points
     * plhs[1]: Qx, the bandwidth in the (rotated) x-axis, 1 x N
     * plhs[2]: Qy, the bandwidth in the (rotated) y-axis, 1 x N
     * plhs[3]: Qz, the bandwidth in the (rotated) z-axis, 1 x N
     * plhs[4]: resn, the residual in the LS fit, 1 x N
     * plhs[5]: lapln, the (normalized) energy of the Laplacian, 1 x N
     * plhs[6]: lopt, the optimal value of lambda using GCV, 1 x N
     *
     */

    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:atti2hydidsi_:callstack","This function should only be called from atti2hydidsi");
    else if( strcmp(callerFunc,"atti2hydidsi") )
        mexErrMsgIdAndTxt("MyToolbox:atti2hydidsi_:callstack","This function should only be called from atti2hydidsi");
    //=======================================================================================
    if(nrhs!=13)
        mexErrMsgIdAndTxt("MyToolbox:atti2hydidsi_:nrhs","Exactly 13 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType G = mxGetM(prhs[0]); // The number of gradient directions per voxel
    //=======================================================================================
    HYDIParameters params;
    params.lambda = mxGetScalar(prhs[5]);
    params.ADC0 = mxGetScalar(prhs[8]);
    params.lRth = mxGetScalar(prhs[9]);
    params.tau = mxGetScalar(prhs[10]);
    params.miters = 100;
    mxArray* mxmiters = mxGetField( prhs[11], 0, "miters");
    if( mxmiters!=NULL ){ params.miters = mxGetScalar(mxmiters); }
    params.otol = 1.0e-6;
    mxArray* mxotol = mxGetField( prhs[11], 0, "otol");
    if( mxotol!=NULL ){ params.otol = mxGetScalar(mxotol); }
    params.stol = 1.0e-6;
    mxArray* mxstol = mxGetField( prhs[11], 0, "stol");
    if( mxstol!=NULL ){ params.stol = mxGetScalar(mxstol); }
    params.ctol = 1.0e-8;
    mxArray* mxctol = mxGetField( prhs[11], 0, "ctol");
    if( mxctol!=NULL ){ params.ctol = mxGetScalar(mxctol); }
    //=======================================================================================
    HYDIIOData io;
    // ------ Inputs
    io.atti = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    io.gi = mxGetDoubles(prhs[2]);
    io.qi = mxGetDoubles(prhs[3]);
    io.DFT = mxGetDoubles(prhs[6]);
    io.u = &(mxGetDoubles(prhs[7])[0]);
    io.v = &(mxGetDoubles(prhs[7])[mxGetM(prhs[6])]);
    io.w = &(mxGetDoubles(prhs[7])[2*mxGetM(prhs[6])]);
    // ------ Outputs
    SizeType ltt[3];
    ltt[0] = (SizeType)(mxGetDoubles(prhs[4])[0]);
    ltt[1] = (SizeType)(mxGetDoubles(prhs[4])[1]);
    ltt[2] = (SizeType)(mxGetDoubles(prhs[4])[2]);
    SizeType K = (2*ltt[0]+1)*(2*ltt[1]+1)*(2*ltt[2]+1);
    K = (K+1)/2;
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( K, N, mxREAL );
        io.eap = mxGetDoubles(plhs[0]);
    }
    else
        return;
    if(nlhs>1){
        plhs[1] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.Qx = mxGetDoubles(plhs[1]);
    }
    else
        io.Qx = NULL;
    if(nlhs>2){
        plhs[2] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.Qy = mxGetDoubles(plhs[2]);
    }
    else
        io.Qy = NULL;
    if(nlhs>3){
        plhs[3] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.Qz = mxGetDoubles(plhs[3]);
    }
    else
        io.Qz = NULL;
    if(nlhs>4){
        plhs[4] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.resn = mxGetDoubles(plhs[4]);
    }
    else
        io.resn = NULL;
    if(nlhs>5){
        plhs[5] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.lapln = mxGetDoubles(plhs[5]);
    }
    else
        io.lapln = NULL;
    if(nlhs>6){
        plhs[6] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.lopt = mxGetDoubles(plhs[6]);
    }
    else
        io.lopt = NULL;
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[12]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, 1 );
    // Own values:
    threader.G          = G;
    threader.K          = K;
    threader.dftdim[0]  = mxGetM(prhs[6]);
    threader.dftdim[1]  = mxGetN(prhs[6]);
    threader.lattice[0] = ltt[0];
    threader.lattice[1] = ltt[1];
    threader.lattice[2] = ltt[2];
    threader.io         = &io;
    threader.params     = &params;
    //=======================================================================================
    // Do the threaded job:
    threader.threadedProcess( maxthreads, atti2hydidsi_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET atti2hydidsi_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    HYDIIOData* io = args->io;
    HYDIParameters* params = args->params;

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads_thread(1);

    // -----------------------------------------------------------------------------
    // Convenience constants:
    SizeType G = args->G;
    SizeType K = args->K;
    // Allocate auxiliar buffers for computations
    // -----------------------------------------------------------------------------
    // Auxiliary buffers to compute DTI spectrum:
    ElementType dti[6];
    ElementType eigval[3];
    ElementType eigvec[9];
    // -----------------------------------------------------------------------------
    // ---
    BufferType gi2 = new ElementType[(args->G)*3];
    // ---
    IndexBuffer used = new IndexType[args->G];
    BufferType qx = new ElementType[args->G];
    BufferType qy = new ElementType[args->G];
    BufferType qz = new ElementType[args->G];
    // ---
    const SizeType NR = (args->dftdim[1]);
    BufferType rx = new ElementType[NR];
    BufferType ry = new ElementType[NR];
    BufferType rz = new ElementType[NR];
    // ---
    BufferType encode = new ElementType[(args->G)*NR];
    BufferType atti = new ElementType[args->G];
    // ---
    BufferType dft = new ElementType[(args->dftdim[0])*(args->dftdim[1])];
    BufferType u = new ElementType[args->dftdim[0]];
    BufferType v = new ElementType[args->dftdim[0]];
    BufferType w = new ElementType[args->dftdim[0]];
    // ---
    BufferType H  = new ElementType[NR*NR];
    BufferType Hi = new ElementType[NR*NR];
    // -----------------------------------------------------------------------------
    // If GCV is required, we will need some extra stuff:
    BLAS_INT info = 0;
    int result_gcv;
    double lambdagcv;
    double costgcv;
    BufferType pinv = NULL;
    bool use_gcv = ( (params->lambda<=0.0f) || mxIsNaN(params->lambda) || mxIsInf(params->lambda) );
    gcv::GCVParams gcvparams;
    if(use_gcv){
        gcv::allocateGCVMemory( G, NR, &gcvparams );
        pinv = new ElementType[NR*NR];
        gcvparams.lambda0 = 0.1f;
        gcvparams.lambdastep = exp(-log(10)/5.0f);
        gcvparams.maxvals = 20;
    }
    // -----------------------------------------------------------------------------
    // Data structures for quadratic programming
    dmriqpp::QPProblem  qpproblem; // The problem itself
    dmriqpp::QPAuxiliar qpaux;     // Auxiliar buffers
    dmriqpp::QPParams   qpparams;  // Parameters for the algorithm
    // Initallize these structures. There are NR variables to optimize,
    // with CE=1 equality constraints, CI=0 inequality constraints, CL=NR
    // lower bounds and CU=0 upper bounds:
    dmriqpp::allocateQPProblem( qpproblem, NR, 1, 0, NR, 0 );
    // The lower bounds are all zeros, so that we can fix them here:
    mataux::setValueMxArray( qpproblem.lb, NR, 1, 0.0 );
    // Finally, set the parameters to the algorithm:
    qpproblem.step = 0.1;
    qpparams.maxiters = 1000000;
    qpparams.streak = 3;
    qpparams.steptol = 1.0e-6;
    qpparams.costtol = 1.0e-6;
    qpparams.normalize = true;
    qpparams.computes0 = true;
    // Now we can allocate all auxiliar buffers:
    dmriqpp::allocateQPAuxiliar( qpproblem,  qpparams, qpaux );
    // NOTE: The voxel-specific buffers to fill are:
    //   The pseudo-inverse of Q, qpproblem.Qi (NR x NR)
    //   The linear term of the cost function, qpproblem.f (NR x 1)
    //   The factors of the equality constraints, qpproblem.Aeqt (NR x 1)
    // -----------------------------------------------------------------------------
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for(  IndexType i=start; i<end; ++i ){
            // ---------------------------------------------------------------------
            // 1- Check the tensor model and compute eigenvalues and eigenvectors
            //    to shape the transformed (anatomical) space (will use Lapack's
            //    dsyev)
            memcpy( dti, &(io->dti[6*i]), 6*sizeof(ElementType) );
            dtisc::sanityCheckDTI( (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, args->params->ADC0, 'A' );
            //---------------------------------------------------------------------
            // 2- Determine the supports of the q- and the R-spaces
            mataux::multiplyMxArrays( args->io->gi, (BufferType)eigvec, gi2,
                                     args->G, 3, 3 );
            ElementType Qx = (ElementType)(args->lattice[0]) / sqrt( args->params->lRth * eigval[0] );
            ElementType Qy = (ElementType)(args->lattice[1]) / sqrt( args->params->lRth * eigval[1] );
            ElementType Qz = (ElementType)(args->lattice[2]) / sqrt( args->params->lRth * eigval[2] );
            ElementType Q  = Qx*Qy*Qz;
            if( args->io->Qx != NULL )
                args->io->Qx[i] = Qx;
            if( args->io->Qy != NULL )
                args->io->Qy[i] = Qy;
            if( args->io->Qz != NULL )
                args->io->Qz[i] = Qz;
            //---------------------------------------------------------------------
            // 3- Determine which q-samples will be actually used (i.e. eliminate
            //    out-of-bandwidth values):
            SizeType total_used = 0;
            for( IndexType g=0; g<args->G; ++g ){
                qx[g] = (gi2[g])             * (args->io->qi[g]);
                qy[g] = (gi2[g+(args->G)])   * (args->io->qi[g]);
                qz[g] = (gi2[g+2*(args->G)]) * (args->io->qi[g]);
                if( (qx[g]<Qx/2 && qx[g]>-Qx/2) && (qy[g]<Qy/2 && qy[g]>-Qy/2) && (qz[g]<Qz/2 && qz[g]>-Qz/2) ){
                    atti[(IndexType)total_used] = args->io->atti[(args->G)*i+g];
                    used[(IndexType)total_used] = g;
                    total_used++;
                }
            }
            // If the total number of used voxels is just 0, we have
            // nothing to do. We should skip this voxel and continue
            if(total_used==0)
                continue; // Should never be reached
            //---------------------------------------------------------------------
            // 4- Sample the R-domain within a regular lattice
            // The (somehow) weird ordering of the loops in nz, then nx, then ny
            // is for consistency with matlab's meshgrid
            IndexType pos = 0;
            for( IndexType nz=0; nz<=(IndexType)(args->lattice[2]); ++nz ){
                IndexType nx = ( nz==0 ? 0 : -(IndexType)(args->lattice[0]) );
                for( ; nx<=(IndexType)(args->lattice[0]); ++nx ){
                    IndexType ny = ( (nx==0 && nz==0) ? 0 : -(IndexType)(args->lattice[1]) );
                    for( ; ny<=(IndexType)(args->lattice[1]); ++ny ){
                        rx[pos] = (ElementType)nx/Qx;
                        ry[pos] = (ElementType)ny/Qy;
                        rz[pos] = (ElementType)nz/Qz;
                        ++pos;
                    }
                }
            }
            //---------------------------------------------------------------------
            // 5- Create the encoding matrix, total_used x dftdim[1]
            // Each row corresponds to a measured value in the q-space, and each
            // column corresponds to a lattice point. Hence, the actual number
            // of "useful" rows is variable (total_used rows are usable)
            for( IndexType r=0; r<total_used; ++r ){
                for( IndexType c=0; c<NR; ++c ){
                    ElementType carg  = ( qx[used[r]] ) * ( rx[c] );
                    carg             += ( qy[used[r]] ) * ( ry[c] );
                    carg             += ( qz[used[r]] ) * ( rz[c] );
                    carg             *= 2*PI;
                    if(c==0) // First column, DC component of the EAP
                        encode[total_used*c+r] = cos(carg)/Q;
                    else // all non-DC components
                        encode[total_used*c+r] = 2.0*cos(carg)/Q;
                }
            }
            //---------------------------------------------------------------------
            // 6- Create the Laplacian penalty, dftdim[0] x dftdim[1]
            memcpy( dft, args->io->DFT, (args->dftdim[0])*(args->dftdim[1])*sizeof(ElementType) );
            ElementType iQ = pow(Q,-2.0/3.0);
            for( IndexType r=0; r<(IndexType)(args->dftdim[0]); ++r ){
                u[r] = Qx*Qx*(args->io->u[r])*iQ;
                v[r] = Qy*Qy*(args->io->v[r])*iQ;
                w[r] = Qz*Qz*(args->io->w[r])*iQ;
                for( IndexType c=0; c<(IndexType)(args->dftdim[1]); ++c ){
                    ElementType norm = (u[r]+v[r]+w[r])/Q;
                    dft[c*(args->dftdim[0])+r] *= norm;
                }
            }
            //---------------------------------------------------------------------
            // 7- Compute the two addends of the quadratic term of the cost
            //    function, i.e. of matrix H
            mataux::transposeMultiplyMxArray( dft, H, args->dftdim[0], NR );
            mataux::transposeMultiplyMxArray( encode, Hi, total_used, NR );
            //---------------------------------------------------------------------
            // 8- If necessary, run Generalized Cross-Validation
            lambdagcv  = params->lambda;
            result_gcv = -1;
            if(use_gcv){
                result_gcv = gcv::computeGCV(encode,Hi,H,atti,pinv,
                        lambdagcv,costgcv,total_used,NR,&gcvparams);
                // In case this call succeeded (the usual case), pinv
                // already contains the inverse of the quadratic term of
                // the cost function, Qi = H^{-1}, which we can use
                if(result_gcv<0) // Matrix inversion failed. Revert...
                    lambdagcv = gcvparams.lambda0;
            }
            if( args->io->lopt != NULL )
                args->io->lopt[i] = (ElementType)lambdagcv;
            //---------------------------------------------------------------------
            // 9- Set up the final shape of the quadratic programming
            //    problem, i.e.:
            //        H = encode'*encode + lambda*frt'*ftr
            //        b = -encode'*atti
            mataux::scalaropMxArray( H, NR, NR, lambdagcv, mataux::MULTIPLY );
            mataux::addMxArrays( Hi, H, H, NR, NR );
            mataux::multiplyMxArrays( atti, encode, qpproblem.f, 1, total_used, NR );
            mataux::scalaropMxArray( qpproblem.f, NR, 1, -1.0f, mataux::MULTIPLY );
            // If we already didn't (during the computation of the GCV-
            // derived lambda), we need to invert H to obtain Qi:
            info = result_gcv;
            if( info<0 ){
                // Set the identity matrix:
                mataux::setValueMxArray( qpproblem.Qi, NR, NR, 0.0 );
                for( IndexType n=0; n<(IndexType)NR; ++n )
                    qpproblem.Qi[n+NR*n] = 1.0;
                memcpy( Hi, H, NR*NR*sizeof(ElementType) );
                BLAS_INT NR_   = NR;
                BLAS_INT nrhs_ = NR;
                LAPACKCALLFCN(dposv)( "L", &NR_, &nrhs_, Hi, &NR_, qpproblem.Qi, &NR_, &info
#ifdef LAPACK_FORTRAN_STRLEN_END
                                       , 1
#endif
                );
            }
            else
                memcpy( qpproblem.Qi, pinv, NR*NR*sizeof(double) );
            if(info<0){
                // If info is not 0, it means that H cannot be inverted,
                // likely because it is not positive definite. The best we
                // can do in this case is just setting the DTI solution:
                setDTISolution( rx, ry, rz,
                        eigval[0], eigval[1], eigval[2],
                        params->tau, NR, &(args->io->eap[i*NR]) );
            }
            else{
                // H is PD and we can solve the problem. It only remains to
                // set the equality constraint:
                mataux::setValueMxArray( qpproblem.Aeqt, NR, 1, 2.0/Q );
                qpproblem.Aeqt[0] = 1.0/Q;
                qpproblem.beq[0] = 1.0;
            }
            //---------------------------------------------------------------------
            // 10- If H is actually invertible, try to solve the problem:
            if(info==0){
                int qp_result = dmriqpp::solveQuadraticProgram( qpproblem,
                        qpaux, qpparams );
                if(qp_result<0){
                    // If qp_result is < 0, the algorithm failed to solve
                    // the problem, and the best we can do, once again, is
                    // setting the DTI solution of the problem:
                    setDTISolution( rx, ry, rz,
                        eigval[0], eigval[1], eigval[2],
                        params->tau, NR, &(args->io->eap[i*NR]) );
                }
                else{
                    // Otherwise, the solution to the problem is stored in
                    // qpproblem.x, which we can directly pass as the
                    // output:
                    memcpy( &(args->io->eap[i*NR]), qpproblem.x,
                            NR*sizeof(ElementType) );
                }
            }
            //---------------------------------------------------------------------
            // 11- If necessary, compute the residual of the fitting and
            //    the Laplacian penalty
            if( args->io->resn != NULL ){
                // In MATLAB, we do:
                //  resn  = ( atti - encode*eap );
                //  resn  = (resn')*resn;
                // which we can translate to calls to mataux::
                //   atti:   total_used x 1
                //   encode: total_used x NR
                //   eap:    NR x 1
                // Use qx as a temporary buffer to store encode*eap:
                mataux::multiplyMxArrays( encode, &(args->io->eap[i*NR]), qx, total_used, NR, 1 );
                // Subtract the arrays (atti can be re-used safely):
                mataux::subtractMxArrays( qx, atti, atti, total_used, 1 );
                // Compute the norm:
                ElementType resn = mataux::normMxNDArray( atti, total_used );
                args->io->resn[i] = resn;
            }
            if( args->io->lapln != NULL ){
                // In MATLAB, we do:
                //  lapln = dft*eap;
                //  lapln = (lapln')*lapln;
                // which we can translate to calls to mataux::
                //   dft: args->dftdim[0] x NR
                //   eap: NR x 1
                mataux::multiplyMxArrays( dft, &(args->io->eap[i*NR]), u, args->dftdim[0], NR, 1 );
                ElementType lapln = mataux::normMxNDArray( u, args->dftdim[0] );
                args->io->lapln[i] = lapln;
            }
            //---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Free memory previously allocated
    // --------------------------------------------------------------------
    delete[] gi2;
    delete[] used;
    // --------------------------------------------------------------------
    delete[] qx;
    delete[] qy;
    delete[] qz;
    delete[] rx;
    delete[] ry;
    delete[] rz;
    // --------------------------------------------------------------------
    delete[] encode;
    delete[] atti;
    delete[] dft;
    delete[] u;
    delete[] v;
    delete[] w;
    // --------------------------------------------------------------------
    delete[] H;
    delete[] Hi;
    // --------------------------------------------------------------------
    if(use_gcv){
        delete[] pinv;
        gcv::freeGCVMemory( &gcvparams );
    }
    // --------------------------------------------------------------------
    dmriqpp::freeQPProblem( qpproblem );
    dmriqpp::freeQPAuxiliar( qpaux );
    // --------------------------------------------------------------------

    // Revert BLAS threads usage to its default:
    blas_num_threads_thread(blas_threads);

    return (THFCNRET)NULL;
}


void setDTISolution( const BufferType rx, const BufferType ry, const BufferType rz,
                    const ElementType l1, const ElementType l2, const ElementType l3,
                    const ElementType tau, const SizeType NR, BufferType eap )
{
    ElementType fptau = 4.0*PI*tau;
    ElementType den   = sqrt( std::abs(l1*l2*l3)*fptau*fptau*fptau );
    for( IndexType n=0; n<(IndexType)NR; ++n ){
        ElementType x = rx[n] * rx[n] / std::abs(l1);
        ElementType y = ry[n] * ry[n] / std::abs(l2);
        ElementType z = rz[n] * rz[n] / std::abs(l3);
        //---
        ElementType earg = -0.25*(x+y+z)/tau;
        eap[n] = exp(earg)/den;
    }
    return;
}
