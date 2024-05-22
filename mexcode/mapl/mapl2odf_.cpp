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
    BufferType ui;
    // Output:
    BufferType odf;
    
} MAPLIOData;

typedef struct MAPLParameters{
    unsigned int Nmax;    // The maximum degree of Hermite polynomials
    ElementType tau;      // Effective diffusion time
    ElementType ADC0;     // Free-water diffusivity
    ElementType contrast; // The exponent of the radial coordinate in the integral
} MAPLParameters;

class ThArgs : public DMRIThreader
{
public:
    MAPLIOData* io;
    MAPLParameters* params;
};

THFCNRET mapl2odf_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  mapl, the MAPL cofficients, M x N
     * prhs[1]:  dti, the tensor fit of that signal, 6 x N
     * prhs[2]:  ui, the directions table, G x 3
     * prhs[3]:  tau, the effective diffusion time, 1 x 1
     * prhs[4]:  ADC0, the free-water diffusivity at body temperature, 1 x 1
     * prhs[5]:  contrast, the exponent q of r in the int( r^q * P(r*ux,r*uy,r*uz), r=0..infty ) integral, 1 x 1
     * prhs[6],  maxthreads, the maximum number of threads, 1 x 1
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
        mexErrMsgIdAndTxt("MyToolbox:mapl2odf_:callstack","This function should only be called from mapl2odf");
    else if( strcmp(callerFunc,"mapl2odf") )
        mexErrMsgIdAndTxt("MyToolbox:mapl2odf_:callstack","This function should only be called from mapl2odf");
    //=======================================================================================
    if(nrhs!=7)
        mexErrMsgIdAndTxt("MyToolbox:mapl2odf_:nrhs","Exactly 7 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType M = mxGetM(prhs[0]); // The number of basis functions
    SizeType G = mxGetM(prhs[2]); // The number of gradient directions per voxel
    unsigned int Nmax = 0;
    while( mapl::numBasisFunctions(Nmax) < M )
        Nmax += 2;
    if( mapl::numBasisFunctions(Nmax) != M )
        mexErrMsgIdAndTxt("MyToolbox:mapl2odf_:dims","Input 0 doesn't seem a MAPL volume");
    //=======================================================================================
    MAPLParameters params;
    params.Nmax = Nmax;
    params.tau = mxGetScalar(prhs[3]);
    params.ADC0 = mxGetScalar(prhs[4]);
    params.contrast = mxGetScalar(prhs[5]);
    //=======================================================================================
    MAPLIOData io;
    // ------ Inputs
    io.mapl = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    io.ui = mxGetDoubles(prhs[2]);
    io.G = G;
    // ------ Outputs
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( G, N, mxREAL );
        io.odf = mxGetDoubles(plhs[0]);
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
    threader.threadedProcess( maxthreads, mapl2odf_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET mapl2odf_process_fcn( void* inargs )
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
    // Allocate memory for the coefficients of all polynomials:
    BufferType coeffsHn, coeffsHnx, coeffsHny, coeffsHnz;
    hermpols::allocatePolynomialCoeffs( coeffsHn,  Nmax ); // (Nmax+1) x (Nmax+1)
    hermpols::allocatePolynomialCoeffs( coeffsHnx, Nmax ); // (Nmax+1) x (Nmax+1)
    hermpols::allocatePolynomialCoeffs( coeffsHny, Nmax ); // (Nmax+1) x (Nmax+1)
    hermpols::allocatePolynomialCoeffs( coeffsHnz, Nmax ); // (Nmax+1) x (Nmax+1)
    // We only need to compute the actual values for coeffsHn, since
    // all other cofficients will be over-written later on:
    hermpols::computePols( coeffsHn, Nmax );
    // We need to additional buffers to compute polynomials
    // products as vector convolutions:
    BufferType pConv1 = new ElementType[Nmax+1];
    BufferType pConv2 = new ElementType[Nmax+1];
    // Allocate memory for the ODF dictionary:
    BufferType Psi = new ElementType[G*M];
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
            dtisc::sanityCheckDTI( (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, params->ADC0, 'D' );
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
            // ---------------------------------------------------------------------
            // 4- Create the dictionary that will be used to span the ODF from
            //    the Hermite polynomials:
            mapl::computeODFDictionary( &uir[0], &uir[G], &uir[2*G],
                                       coeffsHn, coeffsHnx, coeffsHny, coeffsHnz,
                                       pConv1, pConv2,
                                       Psi, ux, uy, uz, G, Nmax,
                                       params->contrast ); // Psi is G x M
            // ---------------------------------------------------------------------
            // 5- Span the signal from its coefficients
            mataux::multiplyMxArrays( Psi, &(io->mapl[i*M]), &(io->odf[i*G]), G, M, 1 );
            // ---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Revert BLAS threads usage to its default:
    blas_num_threads(blas_threads);

    // Free memory previously allocated
    hermpols::destroyPolynomialCoeffs( coeffsHn  );
    hermpols::destroyPolynomialCoeffs( coeffsHnx );
    hermpols::destroyPolynomialCoeffs( coeffsHny );
    hermpols::destroyPolynomialCoeffs( coeffsHnz );
    delete[] pConv1;
    delete[] pConv2;
    delete[] Psi;
    delete[] uir;
    
    return (THFCNRET)NULL;
}
