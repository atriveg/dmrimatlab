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
//#include "matrix.h"
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
    BufferType ui;
    BufferType ri;
    // Output:
    BufferType eap;
    
} MAPLIOData;

typedef struct MAPLParameters{
    unsigned int Nmax;  // The maximum degree of Hermite polynomials
    ElementType tau;    // Effective diffusion time
    ElementType ADC0;   // Free-water diffusivity
} MAPLParameters;

class ThArgs : public DMRIThreader
{
public:
    MAPLIOData* io;
    MAPLParameters* params;
};

THFCNRET mapl2eap_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  mapl, the MAPL cofficients, M x N
     * prhs[1]:  dti, the tensor fit of that signal, 6 x N
     * prhs[2]:  ui, the directions table, G x 3
     * prhs[3]:  ri, the distances-to-origin vector, G x 1
     * prhs[4]:  tau, the effective diffusion time, 1 x 1
     * prhs[5]:  ADC0, the free-water diffusivity at body temperature, 1 x 1
     * prhs[6],  maxthreads, the maximum number of threads, 1x1
     *
     *  OUTPUTS:
     *
     * plhs[0]: eap, G x N, the reconstructed P(R)
     * 
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:mapl2eap_:callstack","This function should only be called from mapl2eap");
    else if( strcmp(callerFunc,"mapl2eap") )
        mexErrMsgIdAndTxt("MyToolbox:mapl2eap_:callstack","This function should only be called from mapl2eap");
    //=======================================================================================
    if(nrhs!=7)
        mexErrMsgIdAndTxt("MyToolbox:mapl2eap_:nrhs","Exactly 7 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType M = mxGetM(prhs[0]); // The number of basis functions
    SizeType G = mxGetM(prhs[2]); // The number of gradient directions per voxel
    unsigned int Nmax = 0;
    while( mapl::numBasisFunctions(Nmax) < M )
        Nmax += 2;
    if( mapl::numBasisFunctions(Nmax) != M )
        mexErrMsgIdAndTxt("MyToolbox:mapl2eap_:dims","Input 0 doesn't seem a MAPL volume");
    //=======================================================================================
    MAPLParameters params;
    params.Nmax = Nmax;
    params.tau = mxGetScalar(prhs[4]);
    params.ADC0 = mxGetScalar(prhs[5]);
    //=======================================================================================
    MAPLIOData io;
    // ------ Inputs
    io.mapl = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    io.ui = mxGetDoubles(prhs[2]);
    io.ri = mxGetDoubles(prhs[3]);
    io.G = G;
    // ------ Outputs
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( G, N, mxREAL );
        io.eap = mxGetDoubles(plhs[0]);
    }
    else
        return;
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[6]) );
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
    threader.threadedProcess( maxthreads, mapl2eap_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET mapl2eap_process_fcn( void* inargs )
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
    unsigned int blas_threads = blas_num_threads_thread(1);
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

    // Prepare the computation of the Psi
    // dictionaries
    // rx, ry, rz:   G x 1
    // r{xyz}HnVals: G x (Nmax+1)
    // Psi:          G x M
    BufferType rx, ry, rz, rxHnVals, ryHnVals, rzHnVals, Psi;
    mapl::allocateDictionaryBuffers( rx, ry, rz, rxHnVals, ryHnVals,
                                     rzHnVals, Psi, G, Nmax );
    // -----------------------------------------------------------------------------
    // Other auxiliary buffers:
    BufferType uir = new ElementType[G*3]; // Rotated directions table
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
            // 2- Compute the normalization factors ux, uy, and uz
            double ux = sqrt(2.0*eigval[0]*(params->tau));
            double uy = sqrt(2.0*eigval[1]*(params->tau));
            double uz = sqrt(2.0*eigval[2]*(params->tau));
            // ---------------------------------------------------------------------
            // 3- Compute the normalized rx, ry, and rz vectors. This implies:
            //    - Rotating the directions table so that the orientations aligned
            //      with the main diffusion direction (first eigen-vector) become
            //      aligned with the 'x' axis in the anatomical frame
            //          Note we do (ui')*U = ((U')*ui)', as desired
            mataux::multiplyMxArrays( args->io->ui,
                                      (BufferType)eigvec, uir, G, 3, 3 );
            //    - Scaling the coordinates so that r{xyz} becomes r{xyz}/u{xyz},
            //      for u{xyz} = sqrt(2*tau*eigval{xyz})
            for( IndexType g=0; g<G; ++g ){
                rx[g] = uir[g]     * (args->io->ri[g]) / ux;
                ry[g] = uir[g+G]   * (args->io->ri[g]) / uy;
                rz[g] = uir[g+2*G] * (args->io->ri[g]) / uz;
            }
            // ---------------------------------------------------------------------
            // 4- Create the dictionary that will be used to span the signal from
            //    the Hermite polynomials:
            hermpols::evaluateAllPols( coeffsHn, rx, rxHnVals, Nmax, G );
            hermpols::evaluateAllPols( coeffsHn, ry, ryHnVals, Nmax, G );
            hermpols::evaluateAllPols( coeffsHn, rz, rzHnVals, Nmax, G );
            mapl::computePsiDictionary( rx, ry, rz,
                                        rxHnVals, ryHnVals, rzHnVals,
                                        Psi,
                                        ux, uy, uz,
                                        G, Nmax ); // Psi is G x M
            // ---------------------------------------------------------------------
            // 5- Span the signal from its coefficients
            mataux::multiplyMxArrays( Psi, &(io->mapl[i*M]), &(io->eap[i*G]), G, M, 1 );
            // ---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Revert BLAS threads usage to its default:
    blas_num_threads_thread(blas_threads);

    // Free memory previously allocated
    hermpols::destroyPolynomialCoeffs( coeffsHn );
    mapl::destroyDictionaryBuffers( rx, ry, rz,
                          rxHnVals, ryHnVals, rzHnVals, Psi );
    delete[] uir;
    
    return (THFCNRET)NULL;
}
