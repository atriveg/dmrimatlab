/*==========================================================
 * atti2mapl_.c
 *
 * This is a core function to atti2mapl, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024- Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
//#include "matrix.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../threads/threadHelper.h"
#include "hermitePols.h"
#include "maplMaths.h"
#include "../quadprog/dmriquadprog.h"
#include "../mathsmex/sanityCheckDTI.h"
#include "../gcv/compute_gcv.h"

typedef struct MAPLIOData{
    // Input:
    SizeType G; // The number of gradients
    BufferType atti;
    BufferType dti;
    BufferType gi;
    BufferType bi;
    // Outputs:
    BufferType mapl;
    BufferType lambdaopt;
} MAPLIOData;

typedef struct MAPLParameters{
    unsigned int Nmax;  // The maximum degree of Hermite polynomials
    ElementType lambda; // Regularization parameter
    ElementType ADC0;   // Free-water diffusivity
    ElementType tau;    // The diffusion time
    BufferType xyz;     // ci x 3
    SizeType ci;        // number of inequality constraints
} MAPLParameters;

class ThArgs : public DMRIThreader
{
public:
    MAPLIOData* io;
    MAPLParameters* params;
};

void setDTISolution( BufferType mapl, const SizeType M );

THFCNRET atti2mapl_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  atti, the attenuation signal to fit, G x N
     * prhs[1]:  dti, the tensor fit of that signal, 6 x N
     * prhs[2]:  gi, the gradients table, G x 3
     * prhs[3]:  bi, the b-values vector, G x 1
     * prhs[4]:  Nmax, the maximum (integer, even, >=0) degree of Hermite polynomials, 1 x 1
     * prhs[5]:  lambda, the regularization parameter, 1 x 1
     * prhs[6]:  ADC0, the free-water diffusivity at body temperature, 1 x 1
     * prhs[7]:  tau, diffusion time in seconds, 1 x 1
     * prhs[8]:  xyz, the normalized coordinates to impose constraints, ci x 3 or empty
     * prhs[9],  maxthreads, the maximum number of threads, 1x1
     *
     *  OUTPUTS:
     *
     * plhs[0]: mapl, M x N, with M = (Nmax+2)*(Nmax+4)*(2*Nmax+3)/24
     * 
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:atti2mapl_:callstack","This function should only be called from atti2mapl");
    else if( strcmp(callerFunc,"atti2mapl") )
        mexErrMsgIdAndTxt("MyToolbox:atti2mapl_:callstack","This function should only be called from atti2mapl");
    //=======================================================================================
    if(nrhs!=10)
        mexErrMsgIdAndTxt("MyToolbox:atti2mapl_:nrhs","Exactly 9 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType G = mxGetM(prhs[0]); // The number of gradient directions per voxel
    //=======================================================================================
    MAPLParameters params;
    params.Nmax = (unsigned int)( mxGetScalar(prhs[4]) );
    params.lambda = mxGetScalar(prhs[5]);
    params.ADC0 = mxGetScalar(prhs[6]);
    params.tau = mxGetScalar(prhs[7]);
    if( mxGetNumberOfElements(prhs[8])==0 ){
        params.xyz = NULL; 
        params.ci = 0;
    }
    else{
        params.xyz = mxGetDoubles(prhs[8]);
        params.ci  = mxGetM( prhs[8] );
    }
    //=======================================================================================
    MAPLIOData io;
    // ------ Inputs
    io.atti = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    io.gi = mxGetDoubles(prhs[2]);
    io.bi = mxGetDoubles(prhs[3]);
    io.G = G;
    // ------ Outputs
    SizeType M = mapl::numBasisFunctions(params.Nmax);
    if(nlhs==0)
        return;
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( M, N, mxREAL );
        io.mapl = mxGetDoubles(plhs[0]);
    }
    else
        io.mapl = NULL;
    if(nlhs>1){
        plhs[1] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.lambdaopt = mxGetDoubles(plhs[1]);
    }
    else
        io.lambdaopt = NULL;
    if(nlhs>2)
        mexErrMsgIdAndTxt("MyToolbox:atti2mapl_:nlhs","Up to 2 outputs can be returned");
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[9]) );
    unsigned long chunksz = 5;
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, chunksz );
    // Own values:
    threader.io         = &io;
    threader.params     = &params;
    //=======================================================================================
    // Do the threaded job:
    threader.threadedProcess( maxthreads, atti2mapl_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET atti2mapl_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    MAPLIOData* io = args->io;
    const MAPLParameters* params = args->params;
    const unsigned int Nmax = params->Nmax;
    const unsigned int ci = params->ci; // Number of ineq. constraints, if any 
    const unsigned int ce = ( ci>0 ? 1 : 0 ); // Either 1 or 0 equality constraints
    const SizeType G = io->G;
    const SizeType M = mapl::numBasisFunctions( Nmax );
    ElementType tau = params->tau; // seconds, since eigenvalues are mm^2/s
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads_thread(1);
    // -----------------------------------------------------------------------------
    // Buffer to store Q:
    BufferType Q = new double[M*M];
    // -----------------------------------------------------------------------------
    // Auxiliary buffers to compute DTI spectrum:
    ElementType dti[6];
    ElementType eigval[3];
    ElementType eigvec[9];
    // -----------------------------------------------------------------------------
    // Set-up the evaluation of Hermite polinomials
    BufferType coeffsHn;
    hermpols::allocatePolynomialCoeffs( coeffsHn, Nmax ); // (Nmax+1) x (Nmax+1)
    hermpols::computePols( coeffsHn, Nmax );
    // Prepare the computation of the Phi and Psi (if required)
    // dictionaries
    // qx, qy, qz:   G x 1
    // q{xyz}HnVals: G x (Nmax+1)
    // Phi:          G x M
    BufferType qx, qy, qz, qxHnVals, qyHnVals, qzHnVals, Phi;
    mapl::allocateDictionaryBuffers( qx, qy, qz, qxHnVals, qyHnVals,
                                     qzHnVals, Phi, G, Nmax );
    // -----------------------------------------------------------------------------
    // x, y, z:     ci x 1
    // {xyz}HnVals: ci x (Nmax+1)
    // Psi:         ci x M
    BufferType x,  y,  z,  xHnVals,  yHnVals,  zHnVals,  Psi;
    if( ci!=0 ){ // There are positivity constraints
        mapl::allocateDictionaryBuffers( x, y, z, xHnVals, yHnVals,
                                         zHnVals, Psi, ci, Nmax );
        // The constraints are imposed over normalized coordinates xyz, hence
        // we can compute this dictionary just once:
        memcpy( x, params->xyz,          ci*sizeof(ElementType) );
        memcpy( y, &(params->xyz[ci]),   ci*sizeof(ElementType) );
        memcpy( z, &(params->xyz[2*ci]), ci*sizeof(ElementType) );
        hermpols::evaluateAllPols( coeffsHn, x, xHnVals, Nmax, ci );
        hermpols::evaluateAllPols( coeffsHn, y, yHnVals, Nmax, ci );
        hermpols::evaluateAllPols( coeffsHn, z, zHnVals, Nmax, ci );
        mapl::computePsiDictionary( x, y, z, xHnVals, yHnVals, zHnVals,
                                    Psi, 1.0, 1.0, 1.0,  ci, Nmax );
    }
    // -----------------------------------------------------------------------------
    BufferType Bnnn;
    if( ci!=0 ){
        // If there are inequality constraints, we must solve the constrained
        // problem, and we need to compute Bnnn so that E(0) = 1:
        Bnnn = new ElementType[M];
        mapl::computePhi0Values(  Bnnn, Nmax );
    }
    // -----------------------------------------------------------------------------
    // Prepare buffers for the regularization matrix,
    // U:             M x M
    // nx, ny, nz:    M x 1
    // Snm, Tnm, Unm: (Nmax+1) x (Nmax+1)
    BufferType U, Snm, Tnm, Unm;
    unsigned int *nx, *ny, *nz;
    mapl::allocateRegularizationMatrix( U, Snm, Tnm, Unm,
                                        nx, ny, nz, Nmax );
    // nx, ny, and nz are data-independent, hence they can be pre-computed here:
    mapl::computenxyz( nx, ny, nz, Nmax );
    // The same goes for Snm, Tnmm, and Unm
    mapl::computeSTUnm( Snm, Tnm, Unm, Nmax );
    // -----------------------------------------------------------------------------
    // Data structures for quadratic programming (NOTE: if the problem is
    // unconstrained, this will return the unconstrained solution)
    dmriqpp::QPProblem  qpproblem;
    dmriqpp::QPAuxiliar qpaux;
    dmriqpp::QPParams   qpparams;
    // Initallize them:
    dmriqpp::allocateQPProblem( qpproblem, M, ce, ci, 0, 0 );
    qpproblem.step = 0.1;
    qpparams.maxiters = 10000000;
    qpparams.streak = 2;
    qpparams.steptol = 1.0e-6;
    qpparams.costtol = 1.0e-6;
    qpparams.normalize = true;
    qpparams.computes0 = true;
    dmriqpp::allocateQPAuxiliar( qpproblem,  qpparams, qpaux );
    // NOTE: if necessary, both the equality and the inequality constraints
    // can be pre-computed here, since they are always the same thanks to
    // the normalization.
    if(ci!=0){
        // Set the equality constraint that ensures E(0) = 1:
        dmriqpp::setAeq( qpproblem, Bnnn );
        qpproblem.beq[0] = 1.0;
        // Set the inequality constraints: the EAP must be non-negative
        // within the set of pre-defined normalized coordinates x,y,z,
        // where the Psi basis functions evaluate to Psi, with size ci x M
        // Note that dmriqpp admits constraints of the form A*x<=b, hence
        // we need to multiply Psi times -1:
        mataux::scalaropMxArray( Psi, ci, M, -1.0, mataux::MULTIPLY );
        dmriqpp::setA( qpproblem, Psi ); // This stores a copy of Psi
        // Fill the vector of inequalities with zeros:
        mataux::setValueMxArray( qpproblem.b, ci, 1, 0.0 );
        // No longer needed:
        delete[] Bnnn;
        mapl::destroyDictionaryBuffers( x, y, z, xHnVals, yHnVals, zHnVals, Psi );
    }
    // -----------------------------------------------------------------------------
    // In case GCV is needed, some extra buffers are required:
    bool use_gcv = ( (params->lambda<0.0) || mxIsInf(params->lambda) || mxIsNaN(params->lambda) );
    gcv::GCVParams gcv;
    BufferType pinv;
    if( use_gcv ){
        use_gcv = true;
        gcv::allocateGCVMemory( G, M, &gcv );
        pinv = new ElementType[M*M];
        // These values are desgined so that a value of lambda ranging [2.0,0.00002] is
        // swept in descending order logarithmically, using 41 steps, so that lambda(n+1) 
        // approx. = 0.8*lambda(n)
        // As soon as the cost increases from n to n+1, the algorithm stops and
        // lambda(n) is chosen as the optimum
        gcv.lambda0 = 2.0;
        gcv.lambdastep = exp(-log(10)/10.0);
        gcv.maxvals = 41;
    }
    // -----------------------------------------------------------------------------
    // Other auxiliary buffers:
    BufferType gir = new ElementType[G*3]; // Rotated gradients table
    // -----------------------------------------------------------------------------
    // That's all we can pre-compute. From now on, we will perform 
    // voxel-especific operations:
    IndexType start = 0; // First index of data blocks to process
    IndexType end   = 0; // Post-last index of data blocks to process
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // While there is still data not processed by any thread...
        for( IndexType i=start; i<end; ++i ){
            // ---------------------------------------------------------------------
            // 1- Check the tensor model and compute eigenvalues and eigenvectors
            //    to shape the transformed (anatomical) space (will use Lapack's
            //    dsyev)
            memcpy( dti, &(io->dti[6*i]), 6*sizeof(ElementType) );
            dtisc::sanityCheckDTI( (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, args->params->ADC0, 'D' );
            // ---------------------------------------------------------------------
            // 2- Compute the normalized qx, qy, and qz vectors. This implies:
            //    - Rotating the gradients table so that the orientations aligned
            //      with the main diffusion direction (first eigen-vector) become
            //      aligned with the 'x' axis in the anatomical frame
            //          Note we do (gi')*U = ((U')*gi)', as desired
            mataux::multiplyMxArrays( args->io->gi,
                                      (BufferType)eigvec, gir, G, 3, 3 );
            //    - Scaling the coordinates so that q{xyz} becomes 2*pi*u{xyz}*q,
            //      for u{xyz} = sqrt(2*tau*eigval{xyz})
            //         This should be, according to eq. (12) in Ozarslan's paper:
            //           q = sqrt(b)/(2*PI)/sqrt(tau), in [mm^{-1}]
            //          q' = 2*PI*u*q 
            //             = 2*PI*sqrt(2*lambda*tau)*sqrt(b)/(2*PI)/sqrt(tau)
            //             = sqrt(2*lambda*b), [adimensional]
            for( IndexType g=0; g<G; ++g ){
                qx[g] = gir[g]     * sqrt( 2.0 * eigval[0] * (args->io->bi[g]) );
                qy[g] = gir[g+G]   * sqrt( 2.0 * eigval[1] * (args->io->bi[g]) );
                qz[g] = gir[g+2*G] * sqrt( 2.0 * eigval[2] * (args->io->bi[g]) );
            }
            // ---------------------------------------------------------------------
            // 3- Create the dictionary that will be used to span the signal from
            //    the Hermite polynomials:
            hermpols::evaluateAllPols( coeffsHn, qx, qxHnVals, Nmax, G );
            hermpols::evaluateAllPols( coeffsHn, qy, qyHnVals, Nmax, G );
            hermpols::evaluateAllPols( coeffsHn, qz, qzHnVals, Nmax, G );
            mapl::computePhiDictionary( qx, qy, qz,
                                        qxHnVals, qyHnVals, qzHnVals,
                                        Phi, G, Nmax ); // Phi is G x M
            // ---------------------------------------------------------------------
            // 4- Create the regularization matrix:
            double ux = sqrt(2.0*eigval[0]*tau);
            double uy = sqrt(2.0*eigval[1]*tau);
            double uz = sqrt(2.0*eigval[2]*tau);
            mapl::computeRegularizationMatrix( U, Snm, Tnm, Unm,
                                       nx, ny, nz, ux, uy, uz, Nmax ); // M x M
            // ---------------------------------------------------------------------
            // 5- Create the Q = Phi'*Phi matrix
            mataux::transposeMultiplyMxArray( Phi, Q, G, M );    // Phi^T*Phi, M x M
            // ---------------------------------------------------------------------
            // 6- In case GCV is requested, run it here-
            //    NOTE: this piece of code, if run, is especially slow, since
            //          the bottle neck of the algorithm (the inversion of the
            //          Q matrix) has to be repeated many times.
            double lambdagcv = params->lambda;
            double cost = mxGetInf();
            int result_gcv = -1;
            if( use_gcv ){
                result_gcv = gcv::computeGCV( 
                                 Phi,              // G x M
                                 Q,                // M x M
                                 U,                // M x M
                                 &(io->atti[i*G]), // G x 1
                                 pinv,             // M x M
                                 lambdagcv,        // final value of lambda
                                 cost,             // final cost
                                 G,                // number of equations
                                 M,                // number of unknowns
                                 &gcv );           // auxiliar value
                if(result_gcv<0) // Matrix inversion failed. Use a "magical value"
                    lambdagcv = 0.2;
            }
            if( io->lambdaopt != NULL )
                io->lambdaopt[i] = lambdagcv;
            // ---------------------------------------------------------------------
            // 7- Create the Linear Least Squares problem:
            //      min 1/2 || Phi*c - E ||^2 + 1/2 lambda c^TUc
            // =>   min 1/2 c^T*(Phi^T*Phi+lambda*U)*c - (E^T*Phi)*c
            mataux::scalaropMxArray( U, M, M, lambdagcv, mataux::MULTIPLY ); // U*lambda, M x M
            mataux::addMxArrays( Q, U, Q, M, M ); // Phi^T*Phi + U*lambda, M x M
            mataux::multiplyMxArrays( &(io->atti[i*G]), Phi, qpproblem.f, 1, G, M );  // E^T*Phi, 1 x M
            mataux::scalaropMxArray( qpproblem.f, 1, M, -1.0, mataux::MULTIPLY );     // -E^T*Phi, 1 x M
            // ---------------------------------------------------------------------
            // 8- Compute the inverse of Q, Qi, using Cholesky decomposition
            //    This will be used by the quadratic programming algorithm.
            //    NOTE: in case result_gcv==0, Qi is already available in 
            //    pinv as a result of GCV
            BLAS_INT info = result_gcv;
            if( info<0 ){
                BLAS_INT nrhs = M;
                BLAS_INT NR   = M;
                // Set the identity matrix:
                mataux::setValueMxArray( qpproblem.Qi, M, M, 0.0 );
                for( IndexType n=0; n<M; ++n )
                    qpproblem.Qi[n+M*n] = 1.0;
                LAPACKCALLFCN(dposv)( "L", &NR, &nrhs, Q, &NR, qpproblem.Qi, &NR, &info
#ifdef LAPACK_FORTRAN_STRLEN_END
                                      , 1
#endif
                );
            }
            else
                memcpy( qpproblem.Qi, pinv, M*M*sizeof(double) );
            // ---------------------------------------------------------------------
            // 9- Actually solve the problem
            if( info!=0 ){
                // The matrix is not P.D. The best we can do is keep
                // the DTI-base solution
                setDTISolution( qpproblem.x, M );
            }
            else{
                // --------------------------------------------------------
                // Whether the problem is constrained or not, the same 
                // function can be used:
                int qp_result = dmriqpp::solveQuadraticProgram( qpproblem, qpaux, qpparams );
                // If the algorithm succeeded, the solution is already in
                // problem.x. Otherwise, the best we can do is keep the
                // DTI-based solution
                if(qp_result<0)
                    setDTISolution( qpproblem.x, M );
            }
            // ---------------------------------------------------------------------
            // 10- Put the output in place
            memcpy( &(io->mapl[i*M]), qpproblem.x, M*sizeof(ElementType) );
            // ---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Revert BLAS threads usage to its default:
    blas_num_threads_thread(blas_threads);

    // Free memory previously allocated
    delete[] Q;
    hermpols::destroyPolynomialCoeffs( coeffsHn );
    mapl::destroyDictionaryBuffers( qx, qy, qz,
                          qxHnVals, qyHnVals, qzHnVals, Phi );
    mapl::destroyRegularizationMatrix( U, Snm, Tnm, Unm, nx, ny, nz );
    dmriqpp::freeQPProblem( qpproblem );
    dmriqpp::freeQPAuxiliar( qpaux );
    if( use_gcv ){
        gcv::freeGCVMemory( &gcv );
        delete[] pinv;
    }
    delete[] gir;

    return (THFCNRET)NULL;
}

void setDTISolution( BufferType mapl, const SizeType M )
{
    mapl[0] = 1.0;
    for( IndexType m=1; m<M; ++m )
        mapl[m] = 0.0;
    return;
}
