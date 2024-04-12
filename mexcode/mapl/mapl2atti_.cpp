/*==========================================================
 * mapl2atti_.c
 *
 * This is a core function to mapl2atti, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024- Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../threads/threadHelper.h"
#include "hermitePols.h"
#include "maplMaths.h"
#include "../mathsmex/sanityCheckDTI.h"

typedef struct MAPLIOData{
    // Input:
    BufferType mapl;
    BufferType dti;
    SizeType G; // The number of gradients
    BufferType gi;
    BufferType bi;
    // Output:
    BufferType atti;
    
} MAPLIOData;

typedef struct MAPLParameters{
    unsigned int Nmax;  // The maximum degree of Hermite polynomials
    ElementType ADC0;   // Free-water diffusivity
} MAPLParameters;

class ThArgs : public DMRIThreader
{
public:
    MAPLIOData* io;
    MAPLParameters* params;
};

THFCNRET mapl2atti_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  mapl, the MAPL cofficients, M x N
     * prhs[1]:  dti, the tensor fit of that signal, 6 x N
     * prhs[2]:  gi, the gradients table, G x 3
     * prhs[3]:  bi, the b-values vector, G x 1
     * prhs[4]:  ADC0, the free-water diffusivity at body temperature, 1 x 1
     * prhs[5],  maxthreads, the maximum number of threads, 1x1
     *
     *  OUTPUTS:
     *
     * plhs[0]: atti, G x N, the reconstructed E(q)
     * 
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:mapl2atti_:callstack","This function should only be called from mapl2atti");
    else if( strcmp(callerFunc,"mapl2atti") )
        mexErrMsgIdAndTxt("MyToolbox:mapl2atti_:callstack","This function should only be called from mapl2atti");
    //=======================================================================================
    if(nrhs!=6)
        mexErrMsgIdAndTxt("MyToolbox:mapl2atti_:nrhs","Exactly 6 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType M = mxGetM(prhs[0]); // The number of basis functions
    SizeType G = mxGetM(prhs[2]); // The number of gradient directions per voxel
    unsigned int Nmax = 0;
    while( mapl::numBasisFunctions(Nmax) < M )
        Nmax += 2;
    if( mapl::numBasisFunctions(Nmax) != M )
        mexErrMsgIdAndTxt("MyToolbox:mapl2atti_:dims","Input 0 doesn't seem a MAPL volume");
    //=======================================================================================
    MAPLParameters params;
    params.Nmax = Nmax;
    params.ADC0 = mxGetScalar(prhs[4]);
    //=======================================================================================
    MAPLIOData io;
    // ------ Inputs
    io.mapl = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    io.gi = mxGetDoubles(prhs[2]);
    io.bi = mxGetDoubles(prhs[3]);
    io.G = G;
    // ------ Outputs
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( G, N, mxREAL );
        io.atti = mxGetDoubles(plhs[0]);
    }
    else
        return;
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[5]) );
    unsigned long chunksz = N/maxthreads/20;
    if(chunksz<1)
        chunksz = 1;
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, chunksz );
    // Own values:
    threader.io         = &io;
    threader.params     = &params;
    //=======================================================================================
    // Do the threaded job:
    threader.threadedProcess( maxthreads, mapl2atti_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET mapl2atti_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    MAPLIOData* io = args->io;
    const MAPLParameters* params = args->params;
    const unsigned int Nmax = params->Nmax;
    const SizeType G = io->G;
    const SizeType M = mapl::numBasisFunctions( Nmax );
        
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads(1);
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

    // Prepare the computation of the Phi
    // dictionaries
    // qx, qy, qz:   G x 1
    // q{xyz}HnVals: G x (Nmax+1)
    // Phi:          G x M
    BufferType qx, qy, qz, qxHnVals, qyHnVals, qzHnVals, Phi;
    mapl::allocateDictionaryBuffers( qx, qy, qz, qxHnVals, qyHnVals,
                                     qzHnVals, Phi, G, Nmax );
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
        for(  IndexType i=start; i<end; ++i ){ // Process the current chunk
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
            // 4- Span the signal from its coefficients
            mataux::multiplyMxArrays( Phi, &(io->mapl[i*M]), &(io->atti[i*G]), G, M, 1 );
            // ---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Revert BLAS threads usage to its default:
    blas_num_threads(blas_threads);

    // Free memory previously allocated
    hermpols::destroyPolynomialCoeffs( coeffsHn );
    mapl::destroyDictionaryBuffers( qx, qy, qz,
                          qxHnVals, qyHnVals, qzHnVals, Phi );
    delete[] gir;
    
    return (THFCNRET)NULL;
}
