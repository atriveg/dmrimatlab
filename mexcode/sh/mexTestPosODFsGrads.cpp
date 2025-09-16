/*==========================================================
 * mexTestPosODFsGrads.c
 *
 * Code for testing maths of positive ODFs
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "math.h"
#include "../mathsmex/sphericalHarmonics.h"
//#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/posODFsMaths.h"
#include "../mathsmex/mexToMathsTypes.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    //=======================================================================================
    // Inputs:
    //
    // prhs[0]: E, the dMRI measurements, Nx1
    // prhs[1]: psi, the SH coefficients, Kx1
    // prhs[2]: lambda, the convolution factors, (L+1)xN
    // prhs[3]: pshell, a poiter telling which of the M shells each sample belongs to, N x 1
    // prhs[4]: mu, the Lagrange multipleir, 1x1
    // prhs[5]: nu, the Laplacian penalty weight, 1x1
    // prhs[6]: gi, the gradients table, Nx3
    //
    // Outputs:
    //
    // plhs[0]: the Lagrangian, 1x1
    // plhs[1]: the gradient, (K+1)x1
    // plhs[2]: the Hessian, (K+1)x1
    //
    // NO ARGUMENTS-CHECKING IS PERFORMED!!!
    //=======================================================================================
    posODFs::ProblemFeatures features;
    features.N  = mxGetM(prhs[0]);
    features.L  = (unsigned int)( (::sqrt(8*((double)(mxGetM(prhs[1]))-1)+9)-3)/2 );
    mwSize K  = ((features.L)+1)*((features.L)+2)/2;
    mwSize Kp = (2*(features.L)+1)*(2*(features.L)+2)/2;
    features.nu = mxGetScalar(prhs[5]);
    features.Y = new ElementType[(features.N)*Kp];
    //=======================================================================================
    BufferType Gi = mxGetDoubles(prhs[6]);
    double* gx = new double[features.N];
    double* gy = new double[features.N];
    double* gz = new double[features.N];
    for( unsigned int k=0; k<(unsigned int)features.N; ++k ){
        gx[k] = Gi[k];
        gy[k] = Gi[k+features.N];
        gz[k] = Gi[k+2*(features.N)];
    }
    double* theta = new double[features.N];
    double* phi   = new double[features.N];
    shmaths::computeSphericalCoordsFromCartesian( gx, gy, gz, theta, phi, features.N );
    delete[] gx;
    delete[] gy;
    delete[] gz;
    shmaths::computeSHMatrixSymmetric( features.N, theta, phi, 2*(features.L), features.Y );
    delete[] theta;
    delete[] phi;
    //=======================================================================================
    posODFs::WignerSymbols wigner;
    posODFs::createWignerSymbols(&wigner,features.L);
    //=======================================================================================
    IndexBuffer pshell = new IndexType[features.N];
    for( IndexType n=0; n<(IndexType)(features.N); ++n )
        pshell[n] = (IndexType)(mxGetDoubles(prhs[3])[n]);
    posODFs::VoxelFeatures voxel;
    voxel.E = mxGetDoubles(prhs[0]);
    voxel.lambda = mxGetDoubles(prhs[2]);
    voxel.pshell = pshell;
    voxel.psi = mxGetDoubles(prhs[1]);
    voxel.mu = mxGetScalar(prhs[4]);
    //=======================================================================================
    posODFs::GradientInputs ginputs;
    posODFs::createGradientInputs( &ginputs, features.N, features.L );
    //=======================================================================================
    posODFs::HessianInputs  hinputs;
    posODFs::createHessianInputs( &hinputs, features.N, features.L );
    //=======================================================================================
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( 1, 1, mxREAL);
        posODFs::ComputeLagrangian( &features, &wigner, &voxel, &ginputs, mxGetDoubles(plhs[0]) );
    }
    //=======================================================================================
    if(nlhs>1){
        plhs[1] = mxCreateDoubleMatrix( K+1, 1, mxREAL);
        posODFs::ComputeGradient( &features, &wigner, &voxel, &ginputs, &hinputs, mxGetDoubles(plhs[1]) );
    }
    //=======================================================================================
    if(nlhs>2){
        plhs[2] = mxCreateDoubleMatrix( K+1, K+1, mxREAL);
        posODFs::ComputeHessian( &features, &wigner, &voxel, &ginputs, &hinputs, mxGetDoubles(plhs[2]) );
    }
    //=======================================================================================
    BufferType work;
    if(nlhs>3)
        work = new ElementType[features.N+Kp];
    //=======================================================================================
    if(nlhs>3){
        plhs[3] = mxCreateDoubleMatrix( 1, 1, mxREAL);
        posODFs::ComputeResidual( &features, &wigner, &voxel, mxGetDoubles(plhs[3]), work );
    }
    //=======================================================================================
    if(nlhs>4){
        plhs[4] = mxCreateDoubleMatrix( features.N+Kp, K-1, mxREAL);
        posODFs::ComputeJacobian( &features, &wigner, &voxel, mxGetDoubles(plhs[4]), work );
    }
    //=======================================================================================
    if(nlhs>5){
        plhs[5] = mxCreateDoubleMatrix( features.N, Kp, mxREAL);
        memcpy( mxGetDoubles(plhs[5]), features.Y, (features.N)*Kp*sizeof(ElementType) );
    }
    //=======================================================================================
    delete[] features.Y;
    delete[] pshell;
    posODFs::destroyWignerSymbols(&wigner);
    posODFs::destroyGradientInputs( &ginputs );
    posODFs::destroyHessianInputs( &hinputs );
    if(nlhs>3)
        delete[] work;
    //=======================================================================================
}
