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
#include "matrix.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/sphericalHarmonics.h" // Just for the definition of PI
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"
#include "../gcv/compute_gcv.h"

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

typedef struct QuadProgMem{
    BufferType H;
    BufferType Hi;
    BufferType g;
    BufferType dir;
    BufferType x;
    BufferType x0;
    BufferType dx;
    ElementType Q;
    IndexBuffer bnd1;
    IndexBuffer bnd1n;
    IndexType bnd2;
    IndexType bnd2n;
    // ---
    IndexBuffer pivot;
    ptrdiff_t llwork;
    BufferType work;
} QuadProgMem;

void allocateQuadProgMem( QuadProgMem*, const SizeType );

void setupQuadProgMem( QuadProgMem*, const BufferType, const BufferType, const BufferType, const ElementType, const SizeType );

void freeQuadProgMem( QuadProgMem* );

IndexType hydidsiQuadProg( QuadProgMem*, const unsigned int, const ElementType, const ElementType, const ElementType, const SizeType );

THFCNRET atti2hydidsi_process_fcn( void* );

void hydidsiComputeBounds( QuadProgMem*, const SizeType, const ElementType );

void hydidsiPredictBounds( QuadProgMem*, const SizeType );

ElementType hydidsiComputeResidual( const QuadProgMem*, const SizeType );

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
     * prhs[10], optim, a structure with the parameters to the optimization algorithm
     *           (see the help on dmriQuadprog)
     *       optim.miters, maximum number of iterations
     *       optim.otol, tolerance for the optimality condition, 1 x 1
     *       optim.stol, tolerance for the optimization step, 1 x 1
     *       optim.ctol, tolerance for the constraints, 1 x 1
     * prhs[11], maxthreads, the maximum number of threads in POSIX systems, 1 x 1
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
    if(nrhs!=12)
        mexErrMsgIdAndTxt("MyToolbox:atti2hydidsi_:nrhs","Exactly 12 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType G = mxGetM(prhs[0]); // The number of gradient directions per voxel
    //=======================================================================================
    HYDIParameters params;
    params.lambda = mxGetScalar(prhs[5]);
    params.ADC0 = mxGetScalar(prhs[8]);
    params.lRth = mxGetScalar(prhs[9]);
    params.miters = 100;
    mxArray* mxmiters = mxGetField( prhs[10], 0, "miters");
    if( mxmiters!=NULL ){ params.miters = mxGetScalar(mxmiters); }
    params.otol = 1.0e-6;
    mxArray* mxotol = mxGetField( prhs[10], 0, "otol");
    if( mxotol!=NULL ){ params.otol = mxGetScalar(mxotol); }
    params.stol = 1.0e-6;
    mxArray* mxstol = mxGetField( prhs[10], 0, "stol");
    if( mxstol!=NULL ){ params.stol = mxGetScalar(mxstol); }
    params.ctol = 1.0e-8;
    mxArray* mxctol = mxGetField( prhs[10], 0, "ctol");
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
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[11]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, 20 );
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
    unsigned int blas_threads = blas_num_threads(1);

    // Convenience constants:
    SizeType G = args->G;
    SizeType K = args->K;
    // Allocate auxiliar buffers for computations
    // ---
    ElementType dti[6];
    ElementType eigval[3];
    ElementType eigvec[9];
    const ptrdiff_t dim = 3;
    ptrdiff_t info = 0;
    ElementType work[9]; // According to Lapack's docs for dspev 
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
    BufferType b = new ElementType[NR];
    BufferType eap = new ElementType[NR];
    // If GCV is required, we will need some extra stuff:
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
    // ---
    IndexBuffer pivot = new IndexType[NR];
    ptrdiff_t llwork2 = -1;
    ptrdiff_t nrhs = 1;
    ptrdiff_t NR_ = NR;
    ptrdiff_t total_used_ = args->G;
    ptrdiff_t unit_ = 1;
    ElementType alpha_ = 1.0f; // For use with dgemm
    ElementType beta_  = 1.0f; // For use with dgemm
    ElementType dumb;
    dsysv( "U", &NR_, &nrhs, H, &NR_, pivot, b, &NR_, &dumb, &llwork2, &info );
    llwork2 = (ptrdiff_t)dumb;
    BufferType work2 = new ElementType[ llwork2 ];
    // ---
    QuadProgMem mem;
    allocateQuadProgMem( &mem, NR );
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
            // 1- Check the tensor model and compute eigenvalues and eigenvectors
            //    to shape the transformed space (will use Lapack's dsyev)
            memcpy( dti, &(io->dti[6*i]), 6*sizeof(ElementType) );
            dspev( "V", "L", &dim, (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, &dim, (BufferType)work, &info );
            if(info!=0){
                // Computation failed for some reason. Fill eigenvalues with
                // free-water value:
                eigval[0] = eigval[1] = eigval[2] = args->params->ADC0;
                eigvec[0] = eigvec[4] = eigvec[8] = 1.0f;
                eigvec[1] = eigvec[2] = eigvec[3] = eigvec[5] = eigvec[6] =  eigvec[7] = 0.0f;
            }
            else{
                eigval[0] = abs(eigval[0]);
                eigval[1] = abs(eigval[1]);
                eigval[2] = abs(eigval[2]);
                // Sanity checks...
                for( unsigned int p=0; p<3; ++p ){
                    eigval[p] = abs(eigval[p]);
                    if( eigval[p]<(args->params->ADC0)/60 )
                        eigval[p] = (args->params->ADC0)/60;
                    if( eigval[p]>args->params->ADC0 )
                        eigval[p] = args->params->ADC0;
                }
                // Ordering:
                ElementType eval;
                ElementType evec[3];
                for( unsigned int p=0; p<2; ++p ){
                    for( unsigned int q=p+1; q<3; ++q ){
                        if(eigval[p]>eigval[q]){ // Must reorder
                            eval = eigval[p];
                            evec[0] = eigvec[3*p+0];
                            evec[1] = eigvec[3*p+1];
                            evec[2] = eigvec[3*p+2];
                            eigval[p] = eigval[q];
                            eigval[q] = eval;
                            eigvec[3*p+0] = eigvec[3*q+0];
                            eigvec[3*p+1] = eigvec[3*q+1];
                            eigvec[3*p+2] = eigvec[3*q+2];
                            eigvec[3*q+0] = evec[0];
                            eigvec[3*q+1] = evec[1];
                            eigvec[3*q+2] = evec[2];
                        }
                    }
                }
                // Make sure eigvec is a rotation matrix
                // by making e1 = e2 x e3
                eigvec[0] = eigvec[4]*eigvec[8]-eigvec[5]*eigvec[7];
                eigvec[1] = eigvec[5]*eigvec[6]-eigvec[3]*eigvec[8];
                eigvec[2] = eigvec[3]*eigvec[7]-eigvec[4]*eigvec[6];
            }
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
            // 7- Set up the constrained, regularized problem
            // This is something like:
            //    H = encode'*encode + lambda*frt'*ftr
            //    b = -encode'*atti
            // ----
            mataux::transposeMultiplyMxArray( dft, H, args->dftdim[0], NR );
            mataux::transposeMultiplyMxArray( encode, Hi, total_used, NR );
            /** GENERALIZED CROSS-VALIDATION IS CARRIED OUT HERE */
            // Note this piece of code, if run, is especially slow, since
            // the bottle neck of the algorithm (the inversion of the PD
            // matrix) has to be repeated many times
            lambdagcv = params->lambda;
            if(use_gcv){
                result_gcv = gcv::computeGCV(encode,Hi,H,atti,pinv,lambdagcv,costgcv,total_used,NR,&gcvparams);
                // In GCV we will return the final value of the inverted matrix, so that we can replace
                // the "dposv" call later on with a product GCV*eap
                if(result_gcv<0) // Matrix inversion failed. Revert...
                    lambdagcv = gcvparams.lambda0;
            }
            if( args->io->lopt != NULL )
                args->io->lopt[i] = (ElementType)lambdagcv;
            /** END OF GENERALIZED CROSS-VALIDATION */
            mataux::scalaropMxArray( H, NR, NR, lambdagcv, mataux::MULTIPLY );
            mataux::addMxArrays( Hi, H, H, NR, NR );
            mataux::multiplyMxArrays( atti, encode, b, 1, total_used, NR );
            mataux::scalaropMxArray( b, NR, 1, -1.0f, mataux::MULTIPLY );
            memcpy( eap, b, NR*sizeof(ElementType) );
            // -----
            // We need to invert the matrix H; the way it is computed, it is necessarily
            // symmetric and positive SEMI definite. But, if it is invertible, then it
            // is positive definite. Summarizing, we will call Lapack's dposv for a
            // Cholesky factorization based solution; if that fails, we cannot proceed
            //
            if( use_gcv && (result_gcv>=0) ){
                // NOTE: if result_gcv==0, the optimal value of lambda
                // was found and the matrix inversion succeeded. If
                // result_gcv>0, we didn't get a local minimum, but still
                // the GCV cost decreased for all iterations and the
                // matrix inversion was successful, so it makes sense to
                // use the final value obtained
                mataux::multiplyMxArrays( pinv, b, eap, NR, NR, 1 );
            }
            else{
                // If either no GCV was required or result_gcv<0 (matrix
                // inversion failed), we have to explicitly invert the matrix:
                memcpy( Hi, H, NR*NR*sizeof(ElementType) );
                dposv( "L", &NR_, &nrhs, Hi, &NR_, eap, &NR_, &info );
                if(info!=0)
                    continue;
            }
            // -----
            ElementType csum = 0.0f;
            for( IndexType r=0; r<NR; ++r ){
                if(eap[r]<0.0)
                    eap[r] = -eap[r];
                else
                    eap[r] = 0.0f;
                if(r>0)
                    csum += 2.0f * eap[r] / Q;
                else
                    csum += eap[r] / Q;
            }
            if(csum>mxGetEps()){
                for( IndexType r=0; r<NR; ++r )
                    eap[r] /= csum;
            }
            else
                eap[0] = Q;
            //---------------------------------------------------------------------
            // 8- Solve the constrained, regularized problem
            IndexType nit = 0;
            if(params->miters>0){
                setupQuadProgMem( &mem, H, b, eap, Q, NR );
                nit = hydidsiQuadProg( &mem, params->miters, params->otol, params->stol, params->ctol, NR );
                memcpy( &(eap[1]), mem.x0, (unsigned int)(NR-1)*sizeof(ElementType) ); // The (unisgned int) cast prevents recent gcc's to throw an overflow warning
                eap[0] = Q;
                for( IndexType r=1; r<NR; ++r )
                    eap[0] -= 2.0f * eap[r];
                memcpy( &(args->io->eap[i*NR]), eap, NR*sizeof(ElementType) );
            }
            //---------------------------------------------------------------------
            // 9- If necessary, compute the residual of the fitting and
            //    the Laplacian penalty
            if( args->io->resn != NULL ){
                // In MATLAB, we do:
                //  resn  = ( atti - encode*eap );
                //  resn  = (resn')*resn;
                // which we can translate to direct Lapack's calls
                //   atti:   total_used x 1
                //   encode: total_used x NR
                //   eap:    NR x 1
                total_used_ = total_used;
                alpha_ = 1.0f;
                beta_  = -1.0f;
                dgemv( "N", &total_used_, &NR_, &alpha_, encode,
                      &total_used_, eap, &unit_, &beta_, atti, &unit_ );
                ElementType resn = ddot( &total_used_, atti, &unit_, atti, &unit_ );
                /*
                ElementType resn = 0.0f;
                for( IndexType r=0; r<total_used; ++r ){
                ElementType temp = args->io->atti[i*(args->G)+used[r]];
                for( IndexType c=0; c<NR; ++c )
                temp -= encode[total_used*c+r] * (args->io->eap[i*NR+c]);
                resn += (temp*temp);
                }
                */
                args->io->resn[i] = resn;
            }
            if( args->io->lapln != NULL ){
                // In MATLAB, we do:
                //  lapln = dft*eap;
                //  lapln = (lapln')*lapln;
                // which we can translate to direct Lapack's calls
                //   dft: args->dftdim[0] x NR
                //   eap: NR x 1
                total_used_ = args->dftdim[0]; // Re-use the variable
                alpha_ = 1.0f;
                beta_  = 0.0f;
                dgemv( "N", &total_used_, &NR_, &alpha_, dft,
                      &total_used_, eap, &unit_, &beta_, u, &unit_ ); // Re-use "u"
                ElementType lapln = ddot( &total_used_, u, &unit_, u, &unit_ );
                /*
                ElementType lapln = 0.0f;
                for( IndexType r=0; r<(IndexType)(args->dftdim[0]); ++r ){
                ElementType temp = 0.0f;
                for( IndexType c=0; c<(IndexType)NR; ++c )
                temp += dft[c*(args->dftdim[0])+r] * (args->io->eap[i*NR+c]);
                lapln += (temp*temp);
                }
                */
                args->io->lapln[i] = lapln;
            }
            //---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );
   
    // Free memory previously allocated
    delete[] gi2;
    delete[] used;
    delete[] qx;
    delete[] qy;
    delete[] qz;
    delete[] rx;
    delete[] ry;
    delete[] rz;
    delete[] encode;
    delete[] atti;
    delete[] dft;
    delete[] u;
    delete[] v;
    delete[] w;
    delete[] H;
    delete[] Hi;
    delete[] b;
    delete[] eap;
    delete[] pivot;
    delete[] work2;
    if(use_gcv){
        delete[] pinv;
        gcv::freeGCVMemory( &gcvparams );
    }
    freeQuadProgMem( &mem );

    // Revert BLAS threads usage to its default:
    blas_num_threads(blas_threads);
    
    return (THFCNRET)NULL;
}

void allocateQuadProgMem( QuadProgMem* mem, const SizeType NR ){
    // --------------------------------------------
    mem->H = new ElementType[(NR-1)*(NR-1)];
    mem->Hi = new ElementType[(NR-1)*(NR-1)];
    mem->g = new ElementType[NR-1];
    mem->dir = new ElementType[NR-1];
    mem->x = new ElementType[NR-1];
    mem->x0 = new ElementType[NR-1];
    mem->dx = new ElementType[NR-1];
    mem->bnd1 = new IndexType[NR-1];
    mem->bnd1n = new IndexType[NR-1];
    // ----
    mem->pivot = new IndexType[NR-1];
    mem->llwork = -1;
    IndexType info = 0;
    ElementType dumb;
    ptrdiff_t nrhs = 1;
    ptrdiff_t NR_ = NR-1;
    // This is a query call, just to check the proper size of the work buffer:
    dsysv( "L", &NR_, &nrhs, mem->H, &NR_, mem->pivot, mem->g, &NR_, &dumb, &(mem->llwork), &info );
    mem->llwork = (ptrdiff_t)dumb;
    mem->work = new ElementType[ mem->llwork ];
    return;
}

void setupQuadProgMem( QuadProgMem* mem, const BufferType H, const BufferType b, const BufferType eap, const ElementType Q, const SizeType NR ){
    // --------------------------------------------
    mem->Q = Q;
    for( IndexType r=0; r<NR-1; ++r ){
        for( IndexType c=r; c<NR-1; ++c ){
            mem->H[c*(NR-1)+r] = H[(c+1)*NR+(r+1)] + 4.0f*H[0] - 2.0f*H[r+1] - 2.0f*H[c+1];
            mem->H[r*(NR-1)+c] = mem->H[c*(NR-1)+r];
        }
    }
    for( IndexType r=0; r<NR-1; ++r ){
        mem->g[r] = b[r+1] + H[r+1]*Q - 2.0f*(H[0]*Q+b[0]);
    }
    memcpy( mem->x0, &eap[1], (NR-1)*sizeof(ElementType) );
    ptrdiff_t NR_ = NR-1;
    ptrdiff_t nrhs = 1;
    ptrdiff_t info = 0;
    memcpy( mem->Hi, mem->H, (NR-1)*(NR-1)*sizeof(ElementType) );
    memcpy( mem->dir, mem->g, (NR-1)*sizeof(ElementType) );
    // We need to invert mem->H. We will first try Cholesky factorization-based
    // solving, and if that fails we will use a generic solver
    dposv( "L", &NR_, &nrhs, mem->Hi, &NR_, mem->dir, &NR_, &info );
    if(info!=0){
        // Restore the problem:
        memcpy( mem->Hi, mem->H, (NR-1)*(NR-1)*sizeof(ElementType) );
        memcpy( mem->dir, mem->g, (NR-1)*sizeof(ElementType) );
        dsysv( "L", &NR_, &nrhs, mem->Hi, &NR_, mem->pivot, mem->dir, &NR_, mem->work, &(mem->llwork), &info );
    }
}

void freeQuadProgMem( QuadProgMem* mem ){
    delete[] mem->H;
    delete[] mem->Hi;
    delete[] mem->g;
    delete[] mem->dir;
    delete[] mem->x;
    delete[] mem->x0;
    delete[] mem->dx;
    delete[] mem->bnd1;
    delete[] mem->bnd1n;
    delete[] mem->pivot;
    delete[] mem->work;
    return;
}

IndexType hydidsiQuadProg( QuadProgMem* mem, const unsigned int miters, const ElementType otol, const ElementType stol, const ElementType ctol, const SizeType NR ){
    // ------------
    memcpy( mem->x, mem->x0, (NR-1)*sizeof(ElementType) );
    hydidsiComputeBounds( mem, NR, ctol );
    memcpy( mem->x0, mem->x, (NR-1)*sizeof(ElementType) );
    memcpy( mem->bnd1, mem->bnd1n, (NR-1)*sizeof(IndexType) );
    mem->bnd2 = mem->bnd2n;
    ElementType res0 = hydidsiComputeResidual( mem, NR);
    ElementType mu = 0.1f;
    IndexType it;
    for( it=1; it<=(IndexType)miters; ++it ){ // We do this weird indexing for consistency with matlab code
        mataux::addMxArrays( mem->x0, mem->dir, mem->dx, NR-1, 1 );
        mataux::scalaropMxArray( mem->dx, NR-1, 1, -mu, mataux::MULTIPLY );
        // ------------
        hydidsiPredictBounds( mem, NR );
        // ------------
        mataux::addMxArrays( mem->x0, mem->dx, mem->x, NR-1, 1 );
        // ------------
        hydidsiComputeBounds( mem, NR, ctol );
        // ------------
        // res = 0.5*(x'*H*x)+g'*x;
        ElementType res = hydidsiComputeResidual( mem, NR);
        // ------------
        if( res<res0 ){
            memcpy( mem->x0, mem->x, (NR-1)*sizeof(BufferType) );
            memcpy( mem->bnd1, mem->bnd1n, (NR-1)*sizeof(IndexType) );
            mem->bnd2 = mem->bnd2n;
            mu = 2.0f * mu;
            if( abs((res0-res)/res0) < otol ){
                res0 = res;
                break;
            }
            else
                res0 = res;
        }
        else
            mu = 0.5f*mu;
        // ------------
        if( mu<stol )
            break;
        // ------------
    }
    // ------------
    return it;
}

void hydidsiPredictBounds( QuadProgMem* mem, const SizeType NR )
{
    ElementType csum = 0.0f;
    SizeType goods = 0;
    for( IndexType r=0; r<(IndexType)NR-1; ++r ){
        if( (mem->bnd1[r]>0) && (mem->dx[r]<0) ){
            mem->bnd1n[r] = 1;
            mem->dx[r] = 0.0f;
        }
        else{
            mem->bnd1n[r] = 0;
            csum += mem->dx[r];
            ++goods;
        }
        if( (mem->bnd2>0) && (csum>0.0f) && (goods>0) ){
            csum /= (ElementType)goods;
            for( IndexType r=0; r<NR-1; ++r ){
                if(mem->bnd1n[r]==0)
                    mem->dx[r] -= csum;
            }
        }
    }
    return;
}

void hydidsiComputeBounds( QuadProgMem* mem, const SizeType NR, const ElementType ctol )
{
    ElementType csum = 0.0f;
    for( IndexType r=0; r<(IndexType)NR-1; ++r ){
        if( mem->x[r]<ctol ){
            mem->bnd1n[r] = 1;
            mem->x[r] = 0.0f;
        }
        else{
            mem->bnd1n[r] = 0;
            csum += mem->x[r];
        }
    }
    csum *= (2.0f/(mem->Q));
    mem->bnd2n = ( csum>1.0f-ctol ? 1 : 0 );
    if( mem->bnd2n != 0 )
        mataux::scalaropMxArray( mem->x, (NR-1), 1, csum, mataux::DIVIDE );
    return;
}

ElementType hydidsiComputeResidual( const QuadProgMem* mem, const SizeType NR )
{
    // Compute:
    //    res = 0.5*(x'*H*x)+g'*x; % This is the original residual
    // I can use mem->dx as an intermediate buffer. Besides, the
    // operation can be done as:
    //   res = x'*(0.5*H*x+g),
    // where the second term can be directly computed with a
    // single call to Lapack's dsymv
    ptrdiff_t N_ = NR-1;
    ptrdiff_t inc_ = 1;
    ElementType alpha = 0.5f;
    ElementType beta = 1.0f;
    memcpy( mem->dx, mem->g, (NR-1)*sizeof(ElementType) );
    dsymv( "U", &N_, &alpha, mem->H, &N_, mem->x, &inc_, &beta, mem->dx, &inc_ );
    ElementType res = ddot( &N_, mem->dx, &inc_, mem->x, &inc_ );
    /*
    ElementType res = 0.0f;
    for( IndexType r=0; r<(IndexType)NR-1; ++r ){
        res += (mem->g[r])*(mem->x[r]);
        res += 0.5f * (mem->x[r]) * (mem->x[r]) * mem->H[r*(NR-1)+r];
        for( IndexType c=r+1; c<NR-1; ++c )
            res += (mem->x[r]) * (mem->x[c]) * mem->H[c*(NR-1)+r];
    }
    */
    return res;
}
