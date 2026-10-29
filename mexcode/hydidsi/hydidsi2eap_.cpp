/*==========================================================
 * hydidsi2eap.cpp
 *
 * This is a core function to mapl2atti, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024- Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/DSIInterp.h"
#include "../threads/threadHelper.h"
#include "../mathsmex/sanityCheckDTI.h"

typedef struct HYDIDSIIOData{
    // Input:
    BufferType dsi;
    BufferType dti;
    SizeBuffer lattice;
    BufferType Qx;
    BufferType Qy;
    BufferType Qz;
    SizeType G; // The number of gradients
    SizeType M; // The number of samples in the "collapsed" grid
    BufferType ui;
    BufferType ri;
    // Output:
    BufferType eap;
    
} HYDIDSIIOData;

typedef struct HYDIDSIParameters{
    BufferType  nodeCoeffs;
    char        interp;
    ElementType ADC0;   // Free-water diffusivity
} HYDIDSIParameters;

class ThArgs : public DMRIThreader
{
public:
    HYDIDSIIOData* io;
    HYDIDSIParameters* params;
};

THFCNRET hydidsi2eap_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  dsi, the HYDI-DSI cofficients, M x N
     * prhs[1]:  dti, the tensor fit of the signal, 6 x N
     * prhs[2]:  lattice, the lattice, 3 x 1
     * prhs[3]:  Qx, the x-bandwidth, 1 x N
     * prhs[4]:  Qy, the y-bandwidth, 1 x N
     * prhs[5]:  Qz, the z-bandwidth, 1 x N
     * prhs[6]:  ui, the directions table, G x 3
     * prhs[7]:  ri, the distances-to-origin vector, G x 1
     * prhs[8]:  interp, the interpolation type, 1 x 1
     * prhs[9]:  ADC0, the free-water diffusivity at body temperature, 1 x 1
     * prhs[10], maxthreads, the maximum number of threads, 1x1
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
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2eap_:callstack","This function should only be called from hydidsi2eap");
    else if( strcmp(callerFunc,"hydidsi2eap") )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2eap_:callstack","This function should only be called from hydidsi2eap");
    //=======================================================================================
    if(nrhs!=11)
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2eap_:nrhs","Exactly 11 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType M = mxGetM(prhs[0]); // The number of EAP samples in the Cartesian grid functions
    SizeType G = mxGetM(prhs[6]); // The number of gradient directions per voxel
    SizeType lattice[3];
    lattice[0] = (SizeType)(round(mxGetDoubles(prhs[2])[0]));
    lattice[1] = (SizeType)(round(mxGetDoubles(prhs[2])[1]));
    lattice[2] = (SizeType)(round(mxGetDoubles(prhs[2])[2]));
    SizeType M2 = ( (2*lattice[0]+1)*(2*lattice[1]+1)*(2*lattice[2]+1) + 1 )/2;
    if( M2 != M )
        mexErrMsgIdAndTxt("MyToolbox:hydidsi2eap_:dims","The DSI volume in input 0 does not match the lattice size you provided");
    //=======================================================================================
    HYDIDSIParameters params;
    params.interp = mxGetChars(prhs[8])[0];
    params.ADC0 = mxGetScalar(prhs[9]);
    //=======================================================================================
    HYDIDSIIOData io;
    // ------ Inputs
    io.dsi = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    io.lattice = (SizeBuffer)lattice;
    io.Qx = mxGetDoubles(prhs[3]);
    io.Qy = mxGetDoubles(prhs[4]);
    io.Qz = mxGetDoubles(prhs[5]);
    io.ui = mxGetDoubles(prhs[6]);
    io.ri = mxGetDoubles(prhs[7]);
    io.G = G;
    io.M = M;
    // ------ Outputs
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( G, N, mxREAL );
        io.eap = mxGetDoubles(plhs[0]);
    }
    else
        return;
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[10]) );
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
    threader.threadedProcess( maxthreads, hydidsi2eap_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET hydidsi2eap_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    HYDIDSIIOData* io = args->io;
    const HYDIDSIParameters* params = args->params;

    const SizeType G = io->G;
    const SizeType M = io->M;

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
    // If necessary, allocate memory for Thomas' algorithm
    SizeType sizes[3];
    sizes[0] = 2*(io->lattice[1]) + 5; // The 1 <-> 0 swap is not a bug, but consistent with meshgrid
    sizes[1] = 2*(io->lattice[0]) + 5; // The 1 <-> 0 swap is not a bug, but consistent with meshgrid
    sizes[2] = 2*(io->lattice[2]) + 5;
    BufferType* gamma;
    BufferType* theta;
    if(params->interp=='s')
        dsiinterp::allocate_thomas3D( gamma, theta, (SizeBuffer)sizes );
    // -----------------------------------------------------------------------------
    // Other auxiliary buffers:
    BufferType rx  = new ElementType[G];
    BufferType ry  = new ElementType[G];
    BufferType rz  = new ElementType[G];
    BufferType uir = new ElementType[G*3]; // Rotated directions table
    BufferType nodeCoeffs;   // Node coefficients to find
    BufferType nodeValues;   // Node values to interpolate
    const SizeType K = dsiinterp::allocate_extended_grid( nodeValues, io->lattice[0], io->lattice[1], io->lattice[2] );
    dsiinterp::allocate_extended_grid( nodeCoeffs, io->lattice[0], io->lattice[1], io->lattice[2] );
    mataux::setValueMxArray(nodeCoeffs,K,1,0.0);
    mataux::setValueMxArray(nodeValues,K,1,0.0);
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
            dtisc::sanityCheckDTI( (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, args->params->ADC0, 'A' );
            // ---------------------------------------------------------------------
            // 2- Compute the normalized rx, ry, and rz vectors. This implies:
            //    - Rotating the directions table so that the orientations aligned
            //      with the main diffusion direction (first eigen-vector) become
            //      aligned with the 'z' axis in the anatomical frame
            //          Note we do (ui')*U = ((U')*ui)', as desired
            mataux::multiplyMxArrays( args->io->ui,
                                      (BufferType)eigvec, uir, G, 3, 3 );
            //    - Scaling the coordinates along each axis
            for( IndexType g=0; g<G; ++g ){
                rx[g] = uir[g]     * (args->io->ri[g]) * io->Qx[i];
                ry[g] = uir[g+G]   * (args->io->ri[g]) * io->Qy[i];
                rz[g] = uir[g+2*G] * (args->io->ri[g]) * io->Qz[i];
            }
            // ---------------------------------------------------------------------
            // 3- Fill the central part of the grid with the EAP passed as
            //    an input:
            dsiinterp::fill_extended_grid( &(io->dsi[M*i]), nodeValues, io->lattice[0], io->lattice[1], io->lattice[2] );
            // ---------------------------------------------------------------------
            // 4- In case it is necessary, solve for the node coefficients
            //    of the B-spline grid:
            if( params->interp=='s' )
                dsiinterp::thomas3D( nodeValues, gamma, theta, nodeCoeffs, (SizeBuffer)sizes );
            else
                memcpy( nodeCoeffs, nodeValues, K*sizeof(ElementType) );
            // ---------------------------------------------------------------------
            // 5- Interpolate:
            switch(params->interp){
                case 'n':
                    for( IndexType g=0; g<G; ++g ){
                        dsiinterp::interpolate_NN( nodeCoeffs,
                                               io->lattice[0], io->lattice[1], io->lattice[2],
                                               rx[g], ry[g], rz[g], io->eap[G*i+g] );
                    }
                    break;
                case 'l':
                    for( IndexType g=0; g<G; ++g ){
                        dsiinterp::interpolate_LN( nodeCoeffs,
                                               io->lattice[0], io->lattice[1], io->lattice[2],
                                               rx[g], ry[g], rz[g], io->eap[G*i+g] );
                    }
                    break;
                case 's':
                    for( IndexType g=0; g<G; ++g ){
                        dsiinterp::interpolate_BS( nodeCoeffs,
                                               io->lattice[0], io->lattice[1], io->lattice[2],
                                               rx[g], ry[g], rz[g], io->eap[G*i+g] );
                    }
                    break;

            }
            // ---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Revert BLAS threads usage to its default:
    blas_num_threads_thread(blas_threads);

    // Free memory previously allocated
    if(params->interp=='s')
        dsiinterp::free_thomas3D( gamma, theta );
    delete[] rx;
    delete[] ry;
    delete[] rz;
    delete[] uir;
    delete[] nodeCoeffs;
    delete[] nodeValues;
    
    return (THFCNRET)NULL;
}
