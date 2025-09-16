/*==========================================================
 * dti2spectrum_.c
 *
 * This is a core function to dti2spectrum, and should only
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
#include "../mathsmex/mexToMathsTypes.h"
#include "../threads/threadHelper.h"

typedef struct EigIOData{
    // Input:
    BufferType dti;
    // Output:
    BufferType l0;
    BufferType l1;
    BufferType l2;
    BufferType e0;
    BufferType e1;
    BufferType e2;
} EigIOData;

class ThArgs : public DMRIThreader
{
public:
    char* mode;          // Compute (or not) the eigenvectors
    EigIOData* io;
};

THFCNRET dti2spectrum_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     *
     * prhs[0]: dti, the signal to fit, 6 x N
     * prhs[1]: maxthreads, 1x1
     *
     *  OUTPUTS:
     *
     * plhs[0]: l0, 1 x N, the first eigenvalue
     * plhs[1]: l1, 1 x N, the second eigenvalue
     * plhs[2]: l2, 1 x N, the third eigenvalue
     * plhs[3]: e0, 3 x N, the first eigenvector
     * plhs[4]: e1, 3 x N, the second eigenvector
     * plhs[5]: e2, 3 x N, the third eigenvector
     * *
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:dti2spectrum_:callstack","This function should only be called from dti2spectrum");
    else if( strcmp(callerFunc,"dti2spectrum") )
        mexErrMsgIdAndTxt("MyToolbox:dti2spectrum_:callstack","This function should only be called from dti2spectrum");
    //=======================================================================================
    if(nrhs!=2)
        mexErrMsgIdAndTxt("MyToolbox:dti2spectrum_:nrhs","Exactly 2 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    //=======================================================================================
    EigIOData io;
    // ------ Inputs
    io.dti = mxGetDoubles(prhs[0]);
    // ------ Outputs
    if( (nlhs!=3) && (nlhs!=6) )
        mexErrMsgIdAndTxt("MyToolbox:dti2spectrum_:nlhs","This function accepts either 3 or 6 output arguments");
    plhs[0] = mxCreateDoubleMatrix( 1, N, mxREAL );
    plhs[1] = mxCreateDoubleMatrix( 1, N, mxREAL );
    plhs[2] = mxCreateDoubleMatrix( 1, N, mxREAL );
    io.l0 = mxGetDoubles(plhs[0]);
    io.l1 = mxGetDoubles(plhs[1]);
    io.l2 = mxGetDoubles(plhs[2]);
    char mode[2] = "V";
    if( nlhs==6 ){
        plhs[3] = mxCreateDoubleMatrix( 3, N, mxREAL );
        plhs[4] = mxCreateDoubleMatrix( 3, N, mxREAL );
        plhs[5] = mxCreateDoubleMatrix( 3, N, mxREAL );
        io.e0 = mxGetDoubles(plhs[3]);
        io.e1 = mxGetDoubles(plhs[4]);
        io.e2 = mxGetDoubles(plhs[5]);
        mode[0] = 'V';
    }
    else{
        io.e0 = (BufferType)NULL;
        io.e1 = (BufferType)NULL;
        io.e2 = (BufferType)NULL;
        mode[0] = 'N';
    }
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[1]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, 20 );
    // Own values:
    threader.mode = mode;
    threader.io   = &io;
    //=======================================================================================
    threader.threadedProcess( maxthreads, dti2spectrum_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET dti2spectrum_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    EigIOData* io = args->io;

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads_thread(1);

    // Allocate auxiliar buffers for computations
    // ---------------------------------------------------------------
    // Allocate memory to compute eigenvalues and eigenvectors
    ElementType dti[6];
    ElementType eigval[3];
    ElementType eigvec[9];
    const BLAS_INT dim = 3;
    BLAS_INT info = 0;
    ElementType work[9]; // According to Lapack's docs for dspev
    ElementType nanval = NAN;
    // ---------------------------------------------------------------
    // Loop through the voxels
    IndexType start = 0;
    IndexType end   = 0;
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // Process all pixels in the block:
        for( IndexType i=start; i<end; ++i ){
            //---------------------------------------------------------------------
            // Copy the current block:
            memcpy( (BufferType)dti, &(io->dti[6*i]), 6*sizeof(ElementType) );
            //---------------------------------------------------------------------
            // Call Lapack's dspev:
            LAPACKCALLFCN(dspev)( args->mode, "L", &dim, (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, &dim, (BufferType)work, &info
#ifdef LAPACK_FORTRAN_STRLEN_END
                                  , 1, 1
#endif
            );
            //---------------------------------------------------------------------
            // Fill the outputs:
            // NOTE: dspev return the eigenvalues in ascending order, but the DTI
            //       convention is the opposite
            if(info==0){
                io->l0[i] = eigval[2];
                io->l1[i] = eigval[1];
                io->l2[i] = eigval[0];
                if(io->e0!=NULL)
                    memcpy( &(io->e0[3*i]), (BufferType)&(eigvec[6]), 3*sizeof(ElementType) );
                if(io->e1!=NULL)
                    memcpy( &(io->e1[3*i]), (BufferType)&(eigvec[3]), 3*sizeof(ElementType) );
                if(io->e2!=NULL)
                    memcpy( &(io->e2[3*i]), (BufferType)&(eigvec[0]), 3*sizeof(ElementType) );
            }
            else{
                io->l0[i] = nanval;
                io->l1[i] = nanval;
                io->l2[i] = nanval;
                if(io->e0!=NULL)
                    io->e0[3*i] = io->e0[3*i+1] = io->e0[3*i+2] = nanval;
                if(io->e1!=NULL)
                    io->e1[3*i] = io->e1[3*i+1] = io->e1[3*i+2] = nanval;
                if(io->e2!=NULL)
                    io->e2[3*i] = io->e2[3*i+1] = io->e2[3*i+2] = nanval;
            }
            //---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    blas_num_threads_thread(blas_threads);

    return (THFCNRET)NULL;
}
