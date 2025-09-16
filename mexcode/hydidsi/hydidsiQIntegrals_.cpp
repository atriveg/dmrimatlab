/*==========================================================
 * hydidsiQIntegrals_.cpp
 *
 * This is a core function to hydidsi2shodf, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
//#include "matrix.h"
#include "math.h"
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/sphericalHarmonics.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

typedef struct HYDIIntegralsIOData{
    // Input:
    BufferType eap;
    BufferType dti;
    BufferType Qx;
    BufferType Qy;
    BufferType Qz;
    // Output:
    BufferType integrals;
} HYDIIntegralsIOData;

class ThArgs : public DMRIThreader
{
public:
    SizeType NR;         // The number of EAP points, must match the lattice
    SizeType G;          // The number of gradient directions to evaluate
    SizeType lattice[3]; // The lattice size
    ElementType ADC0;    // Free-water diffusivity
    BufferType grads;    // The unit-norm gradient directions
    HYDIIntegralsIOData* io;
};

THFCNRET hydidsiQIntegrals_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     *
     * prhs[0]:  eap, the EAP as computed with atti2hydidsi, NR x N
     * prhs[1]:  dti, the tensor fit of the signal as returned by atti2hydidsi, 6 x N
     * prhs[2]:  Qx, the x-bandwidth as computed with atti2hydidsi, 1 x N
     * prhs[3]:  Qy, the y-bandwidth as computed with atti2hydidsi, 1 x N
     * prhs[4]:  Qz, the z-bandwidth as computed with atti2hydidsi, 1 x N
     * prhs[5]:  lattice, the lattice size, 3 x 1, must be the same used with atti2hydidsi
     * prhs[6],  ADC0, the free-water diffusivity, 1 x 1
     * prhs[7]:  Gi, the gradient directions where the integrals will be evaluated, G x 3, unit norm
     * prhs[8],  maxthreads, the maximum number of threads in POSIX systems, 1 x 1
     *
     *  OUTPUTS:
     *
     * plhs[0]: integrals, the integrals at the desired directions in Gi, G x N
     *
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:hydidsiQIntegrals_:callstack","This function should only be called from hydidsi2shodf");
    else if( strcmp(callerFunc,"hydidsi2shodf") )
        mexErrMsgIdAndTxt("MyToolbox:hydidsiQIntegrals_:callstack","This function should only be called from hydidsi2shodf");
    //=====================================hydidsiQIntegrals_==================================================
    if(nrhs!=9)
        mexErrMsgIdAndTxt("MyToolbox:hydidsiQIntegrals_:nrhs","Exactly 9 input arguments are required");
    if(nlhs!=1)
        mexErrMsgIdAndTxt("MyToolbox:hydidsiQIntegrals_:nlhs","Exactly 1 output argument is required");
    //=======================================================================================
    SizeType N  = mxGetN(prhs[0]); // The number of voxels to process
    SizeType NR = mxGetM(prhs[0]); // The number of EAP points
    SizeType G  = mxGetM(prhs[7]); // The number of gradient directions to evalaute
    //=======================================================================================
    HYDIIntegralsIOData io;
    // ------ Inputs
    io.eap = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    io.Qx  = mxGetDoubles(prhs[2]);
    io.Qy  = mxGetDoubles(prhs[3]);
    io.Qz  = mxGetDoubles(prhs[4]);
    // ------ Outputs
    plhs[0] = mxCreateDoubleMatrix( G, N, mxREAL );
    io.integrals = mxGetDoubles(plhs[0]);
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[8]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, 20 );
    // Own values:
    threader.NR         = NR;
    threader.G          = G;
    threader.lattice[0] = (SizeType)mxGetDoubles(prhs[5])[0];
    threader.lattice[1] = (SizeType)mxGetDoubles(prhs[5])[1];
    threader.lattice[2] = (SizeType)mxGetDoubles(prhs[5])[2];
    threader.ADC0       = mxGetDoubles(prhs[6])[0];
    threader.grads      = mxGetDoubles(prhs[7]);
    threader.io         = &io;
    //=======================================================================================
    threader.threadedProcess( maxthreads, hydidsiQIntegrals_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET hydidsiQIntegrals_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    HYDIIntegralsIOData* io = args->io;

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads_thread(1);

    // Convenience constants:
    SizeType NR = args->NR;
    SizeType G  = args->G;
    // Allocate auxiliar buffers for computations
    // ---
    ElementType dti[6];
    ElementType eigval[3];
    ElementType eigvec[9];
    const BLAS_INT dim = 3;
    BLAS_INT info = 0;
    ElementType work[9]; // According to Lapack's docs for dspev
    // ---
    BufferType rx = new ElementType[NR];
    BufferType ry = new ElementType[NR];
    BufferType rz = new ElementType[NR];
    BufferType rs = new ElementType[NR];
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
            LAPACKCALLFCN(dspev)( "V", "L", &dim, (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, &dim, (BufferType)work, &info
#ifdef LAPACK_FORTRAN_STRLEN_END
                                  , 1, 1
#endif
            );
            if(info!=0){
                // Computation failed for some reason. Fill eigenvalues with
                // free-water value:
                eigval[0] = eigval[1] = eigval[2] = args->ADC0;
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
                    if( eigval[p]<(args->ADC0)/60 )
                        eigval[p] = (args->ADC0)/60;
                    if( eigval[p]>args->ADC0 )
                        eigval[p] = args->ADC0;
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
            // 2- Compute the points where the EAP is evaluated at these voxels
            ElementType Qx = args->io->Qx[i];
            ElementType Qy = args->io->Qy[i];
            ElementType Qz = args->io->Qz[i];
            ElementType Q = Qx*Qy*Qz;
            IndexType pos = 0;
            for( IndexType nz=0; nz<=(IndexType)(args->lattice[2]); ++nz ){
                IndexType nx = ( nz==0 ? 0 : -(IndexType)(args->lattice[0]) );
                for( ; nx<=(IndexType)(args->lattice[0]); ++nx ){
                    IndexType ny = ( (nx==0 && nz==0) ? 0 : -(IndexType)(args->lattice[1]) );
                    for( ; ny<=(IndexType)(args->lattice[1]); ++ny ){
                        rx[pos] = (ElementType)nx/Qx;
                        ry[pos] = (ElementType)ny/Qy;
                        rz[pos] = (ElementType)nz/Qz;
                        rs[pos] = rx[pos]*rx[pos] + ry[pos]*ry[pos] + rz[pos]*rz[pos];
                        ++pos;
                    }
                }
            }
            //---------------------------------------------------------------------
            // 3- Loop through the gradient directions to compute the integrals:
            for( IndexType g=0; g<(IndexType)G; ++g ){
                // 3.1. Rotate the gradient according to
                // the tensor model at this voxel
                ElementType u[3];
                ElementType up[3];
                u[0] = args->grads[g];
                u[1] = args->grads[G+g];
                u[2] = args->grads[2*G+g];
                mataux::multiplyMxArrays(u,eigvec,up,1,3,3);
                // 3.2. Determine the support of the integral:
                ElementType T = (Q/2)/mxGetEps();
                ElementType lim;
                if( abs(up[0])>mxGetEps() ){
                    lim = (Qx/2)/abs(up[0]);
                    if(lim<T)
                        T = lim;
                }
                if( abs(up[1])>mxGetEps() ){
                    lim = (Qy/2)/abs(up[1]);
                    if(lim<T)
                        T = lim;
                }
                if( abs(up[2])>mxGetEps() ){
                    lim = (Qz/2)/abs(up[2]);
                    if(lim<T)
                        T = lim;
                }
                // 3.3. Loop trough the EAP samples to compute the signal
                //  (the first EAP sample can be skipped since it corresponds
                //  to the origin, where ra[0]=0).
                ElementType value = 0.0f;
                for( IndexType k=1; k<(IndexType)NR; ++k ){
                    // This stands for:
                    //     integral [from 0 to T ] R_k^2 * 2/Q * q * cos( 2*Pi*q*(u^T*R_k) ) dq
                    //   = integral [from 0 to T ] A * q * cos( B*q ) dq
                    //   = (A/B^2) ( B*T*sin(B*T) + cos(B*T) - 1 )
                    ElementType A = rs[k] * 2.0f / Q;
                    ElementType B = 2.0f * PI * ( up[0]*rx[k] + up[1]*ry[k] + up[2]*rz[k] );
                    value += (args->io->eap[i*NR+k]) * (A/(B*B)) * ( B*T*sin(B*T) + cos(B*T) - 1.0f );
                }
                // 3.4. Assign the output:
                args->io->integrals[G*i+g] = value;
            }
            //---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Free memory previously allocated
    delete[] rx;
    delete[] ry;
    delete[] rz;
    delete[] rs;

    // Revert BLAS threads usage to its default:
    blas_num_threads_thread(blas_threads);

    return (THFCNRET)NULL;
}
