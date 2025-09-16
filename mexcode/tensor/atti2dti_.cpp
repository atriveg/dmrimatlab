/*==========================================================
 * atti2dti_.c
 *
 * This is a core function to atti2dti, and should only
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

typedef struct DTIIOData{
    // Input:
    BufferType signal;
    BufferType gi;
    BufferType bi;
    // Output:
    BufferType dti;
    BufferType s0;
} DTIIOData;

typedef struct DTIParameters{
    unsigned int wlsit;
    unsigned int maxiters;
    ElementType tol;
    ElementType wsc;
    ElementType rcondth;
    char fixmode; // 'n'/'z'/'a'
    char mode;    // 'o'/'w'/'p'
} DTIParameters;

class ThArgs : public DMRIThreader
{
public:
    SizeType G;          // The number of gradients
    DTIIOData* io;
    DTIParameters* params;
};


THFCNRET signal2dti_process_fcn( void* );

void squareDTI( const BufferType, BufferType );

void DTIFromEigs( const ElementType, const ElementType, const ElementType, const BufferType, BufferType );

void squaredDTIDerivatives( const BufferType, BufferType );

IndexType dtiLMFit( const SizeType, const BufferType, BufferType, BufferType,
                    const BufferType, BufferType, BufferType, BufferType,
                    BufferType, BufferType, BufferType,
                    IndexBuffer, IndexBuffer, SizeType, BufferType,
                    const DTIParameters* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  signal, the signal to fit, G x N
     * prhs[1]:  gi, the gradients table, G x 3
     * prhs[2]:  bi, the q-values vector, G x 1
     * prhs[3]:  opts, a structure with the parameters to the algorithm
     *       opts.wlsit, the number of iterations with WLS, 1 x 1
     *       opts.maxiters, the maximum number of iterations for the nonlinear fitting, 1 x 1
     *       opts.wsc: the minimum weight to be applied in the WLS problem with
     *                 respect to the maximum one, 1 x 1
     *       opts.tol, the minimum allowed realtive change in the solution for the nonlinear fitting, 1 x 1
     *       opts.rcondth, for WLS, the minimum condition number so that the WLS matrix can be inverted, 1 x 1
     *       opts.fixmode, the way negative eigenvalues are fixed, char:
     *            'n': none, nothing is done
     *            'z': zero-clip, they are changed to 0
     *            'a': absolute value
     *       opts.mode, the algorithm to be used:
     *            'o': ordinary LS (log domain)
     *            'w', weigthed LS (log domain)
     *            'p', (semi)positive define tensor (natural domain)
     * prhs[4], maxthreads, the maximum number of threads in POSIX systems, 1 x 1
     *
     *  OUTPUTS:
     *
     * plhs[0]: dti, 6 x N, the unique components of the diffusion tensor
     * plhs[1]: s0, 1 x N, the (normalized) baseline estimated, will be near 1 in all cases. Note
     *    it is necessarily normalize since we pass S_i/S_0 to the function instead of S_i AND S_0
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:atti2dti_:callstack","This function should only be called from atti2dti");
    else if( strcmp(callerFunc,"atti2dti") )
        mexErrMsgIdAndTxt("MyToolbox:atti2dti_:callstack","This function should only be called from atti2dti");
    //=======================================================================================
    if(nrhs!=5)
        mexErrMsgIdAndTxt("MyToolbox:atti2dti_:nrhs","Exactly 5 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType G = mxGetM(prhs[0]); // The number of gradient directions per voxel
    //=======================================================================================
    DTIParameters params;
    // --
    params.wlsit = 5;
    mxArray* mxwlsit = mxGetField( prhs[3], 0, "wlsit");
    if( mxwlsit!=NULL ){ params.wlsit = (unsigned int)mxGetScalar(mxwlsit); }
    // --
    params.maxiters = 100;
    mxArray* mxmaxiters = mxGetField( prhs[3], 0, "maxiters");
    if( mxmaxiters!=NULL ){ params.maxiters = (unsigned int)mxGetScalar(mxmaxiters); }
    // --
    params.wsc = 0.01;
    mxArray* mxwsc = mxGetField( prhs[3], 0, "wsc");
    if( mxwsc!=NULL ){ params.wsc = mxGetScalar(mxwsc); }
    // --
    params.tol = 1.0e-6;
    mxArray* mxtol = mxGetField( prhs[3], 0, "tol");
    if( mxtol!=NULL ){ params.tol = mxGetScalar(mxtol); }
    // --
    params.rcondth = 1.0e-6;
    mxArray* mxrcondth = mxGetField( prhs[3], 0, "rcondth");
    if( mxrcondth!=NULL ){ params.rcondth = mxGetScalar(mxrcondth); }
    // --
    params.fixmode = 'n';
    mxArray* mxfixmode = mxGetField( prhs[3], 0, "fixmode");
    if( mxfixmode!=NULL ){ params.fixmode = (char)(mxGetChars(mxfixmode)[0]); }
    // --
    params.mode = 'w';
    mxArray* mxmode = mxGetField( prhs[3], 0, "mode");
    if( mxmode!=NULL ){ params.mode = (char)(mxGetChars(mxmode)[0]); }
    //=======================================================================================
    DTIIOData io;
    // ------ Inputs
    io.signal = mxGetDoubles(prhs[0]);
    io.gi = mxGetDoubles(prhs[1]);
    io.bi = mxGetDoubles(prhs[2]);
    // ------ Outputs
    if(nlhs!=2)
        mexErrMsgIdAndTxt("MyToolbox:atti2dti_:nlhs","This function accepts just 2 output arguments");
    plhs[0] = mxCreateDoubleMatrix( 6, N, mxREAL );
    io.dti = mxGetDoubles(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix( 1, N, mxREAL );
    io.s0 = mxGetDoubles(plhs[1]);
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[4]) );
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, 20 );
    // Own values:
    threader.G      = G;
    threader.io     = &io;
    threader.params = &params;
    //=======================================================================================
    threader.threadedProcess( maxthreads, signal2dti_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET signal2dti_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    DTIIOData* io = args->io;
    DTIParameters* params = args->params;

    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance. In non-POSIX systems, however, we don't
    // externally create threads and we can let Open MP do its
    // stuff.
    unsigned int blas_threads = blas_num_threads_thread(1);
    
    // Convenience constants:
    SizeType G = args->G;
    // Allocate auxiliar buffers for computations
    // ---------------------------------------------------------------
    ElementType dti[6];
    ElementType x[7];
    ElementType x2[7];
    ElementType xr[7];
    ElementType xr2[7];
    ElementType eigval[3];
    ElementType eigvec[9];
    const BLAS_INT dim = 3;
    BLAS_INT info = 0;
    ElementType work[9]; // According to Lapack's docs for dspev 
    // ---------------------------------------------------------------
    // The gradients matrix:
    BufferType lSi = new ElementType[G+1];
    BufferType lSiw = new ElementType[G+1];
    ElementType lSii[7];
    BufferType A = new ElementType[(G+1)*7];
    BufferType A2 = new ElementType[(G+1)*7];
    mataux::setValueMxArray( A, G+1, 7, 0.0f );
    for( IndexType g=0; g<(IndexType)G; ++g ){
        A[g]         = -(io->bi[g]) * (io->gi[g]) * (io->gi[g]);
        A[(G+1)+g]   = -2.0f * (io->bi[g]) * (io->gi[g]) * (io->gi[G+g]);
        A[(G+1)*2+g] = -2.0f * (io->bi[g]) * (io->gi[g]) * (io->gi[G*2+g]);
        A[(G+1)*3+g] = -(io->bi[g]) * (io->gi[G+g]) * (io->gi[G+g]);
        A[(G+1)*4+g] = -2.0f * (io->bi[g]) * (io->gi[G+g]) * (io->gi[G*2+g]);
        A[(G+1)*5+g] = -(io->bi[g]) * (io->gi[G*2+g]) * (io->gi[G*2+g]);
    }
    // ---
    ElementType scale = 0.0f;
    for( IndexType g=0; g<(IndexType)G; ++g ){
        for( unsigned int c=0; c<6; ++c )
            scale += A[(G+1)*c+g] * A[(G+1)*c+g];
    }
    scale /= (ElementType)(G*6);
    scale  = sqrt(scale);
    // ---
    for( IndexType g=0; g<(IndexType)G+1; ++g )
        A[(G+1)*6+g] = scale;
    // ---------------------------------------------------------------
    // Compute the fixed fitting matrix for ordinary LS:
    ElementType ATA[7*7];
    ElementType ATA2[7*7];
    ElementType ATA3[7*7];
    mataux::transposeMultiplyMxArray( A, (BufferType)ATA, G+1, 7 );
    // ---
    IndexType pivot1[7];
    IndexType pivot2[7];
    BLAS_INT flag = 1;
    BLAS_INT K = 7;
    BLAS_INT lname = 6; // Length of "dgetri"
    BLAS_INT largs = 0; // Length of ""
#ifdef OCTAVE_BUILD
    // ilaenv seems not to be present in liblapack. Use BS=4 as a
    // "one size fits all" thing:
    BLAS_INT BS = 4;
#else
    BLAS_INT BS  = LAPACKCALLFCN(ilaenv)(
            &flag, "dgetri", "", 
            &K, &K, &K, &K, lname, largs );
    BS = ( BS>4 ? BS : 4 );
#endif
    SizeType lwork = BS*7;
    BufferType work2 = new ElementType[lwork];
    mataux::checkAndInvertMxArray( (BufferType)ATA, 7, mxGetEps()*10, (IndexBuffer)pivot1, (IndexBuffer)pivot2, lwork, work2 );
    // ---
    BufferType AT  = new ElementType[7*(G+1)];
    mataux::transposeMxArray( A, AT, G+1, 7 );
    // ---
    BufferType piA = new ElementType[7*(G+1)];
    mataux::multiplyMxArrays( ATA, AT, piA, 7, 7, G+1 );
    // ---------------------------------------------------------------
    BufferType wi = new ElementType[G+1];
    // ---------------------------------------------------------------
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
            // Regardless of the working mode, we will always begin with an OLS
            // approximation to the solution, which reduces to the product of the
            // pre-compute pseudo-inverse of A, piA, with the acquired signal:
            memcpy( lSi, &(io->signal[i*G]), G*sizeof(ElementType) );
            for( IndexType g=0; g<(IndexType)G; ++g )
                lSi[g] = log( lSi[g]>mxGetEps() ? lSi[g] : mxGetEps() );
            lSi[G] = 0.0f;
            mataux::multiplyMxArrays( piA, lSi, (BufferType)x, 7, G+1, 1 );
            //---------------------------------------------------------------------
            if( params->mode=='w' || params->mode=='p' ){
                // In any of these cases we need to iterate to compute the
                // WLS solution
                for( unsigned int n=0; n<params->wlsit; ++n ){
                    // 1- Compute the weights wi at this iteration:
                    mataux::multiplyMxArrays( A, (BufferType)x, wi, G+1, 7, 1 );
                    ElementType maxw = -1.0f;
                    bool anynan = false;
                    for( IndexType g=0; g<(IndexType)G+1; ++g ){
                        wi[g] = exp(2.0f*wi[g]);
                        wi[g] = ( wi[g]>1.0 ? 1.0 : wi[g] );
                        wi[g] = ( wi[g]<0.0 ? 0.0 : wi[g] );
                        maxw = ( wi[g]>maxw ? wi[g] : maxw );
                        anynan |= isnan(wi[g]);
                    }
                    if(anynan)
                        break;
                    for( IndexType g=0; g<(IndexType)G+1; ++g ){
                        if( wi[g] < maxw*(params->wsc) )
                            wi[g] = maxw*(params->wsc);
                    }
                    wi[G] = maxw;
                    // 2- Create the WLS problem:
                    memcpy( A2, A, (G+1)*7*sizeof(ElementType) );
                    for( IndexType g=0; g<(IndexType)G+1; ++g ){
                        lSiw[g] = lSi[g] * wi[g];
                        for( unsigned int c=0; c<7; ++c )
                            A2[c*(G+1)+g] *= wi[g];
                    }
                    mataux::multiplyMxArrays( AT, A2, ATA, 7, G+1, 7 );
                    mataux::multiplyMxArrays( AT, lSiw, lSii, 7, G+1, 1 );
                    // 3- Solve the WLS problem:
                    IndexType result = mataux::checkAndInvertMxArray( (BufferType)ATA, 7, params->rcondth,
                                                                     (IndexBuffer)pivot1, (IndexBuffer)pivot2, lwork, work2 );
                    if( result == 0 ){
                        mataux::multiplyMxArrays( (BufferType)ATA, lSii, (BufferType)x2, 7, 7, 1 );
                        for( unsigned int c=0; c<7; ++c )
                            anynan |= isnan(x2[c]);
                        if(anynan)
                            break;
                        else
                            memcpy( x, x2, 7*sizeof(ElementType) );
                    }
                    else
                        break;
                }
            }
            //---------------------------------------------------------------------
            // Either non-negative eigenvalues correction is explicitly asked for
            // or positive, non-linear estimation will be used, it is necessary to
            // compute the eigenvalues and eigenvectors (in the latter case, to
            // compute a feasible initial iteration).
            if( params->mode=='p' || params->fixmode!='n' ){
                memcpy( (BufferType)dti, (BufferType)x, 6*sizeof(ElementType) );
                LAPACKCALLFCN(dspev)( "V", "L", &dim, (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, &dim, (BufferType)work, &info
#ifdef LAPACK_FORTRAN_STRLEN_END
                , 1, 1
#endif
                );
                if(info==0){
                    if( params->fixmode=='a' ){
                        eigval[0] = abs( eigval[0] );
                        eigval[1] = abs( eigval[1] );
                        eigval[2] = abs( eigval[2] );
                    }
                    else{
                        eigval[0] = ( eigval[0]>0 ? eigval[0] : 0.0f );
                        eigval[1] = ( eigval[1]>0 ? eigval[1] : 0.0f );
                        eigval[2] = ( eigval[2]>0 ? eigval[2] : 0.0f );
                    }
                    // Reconstruct the (unique components of the) diffusion
                    // tensor from its eigenvalues and eigenvectors:
                    DTIFromEigs( eigval[0], eigval[1], eigval[2], eigvec, x );
                }
            }
            //---------------------------------------------------------------------
            if( params->mode=='p' && info==0 ){
                // NOTE: if info!=0, the eigenvalue computation failed and
                // we cannot assure there is a feasible first iteration
                DTIFromEigs( sqrt(eigval[0]), sqrt(eigval[1]), sqrt(eigval[2]), eigvec, xr );
                xr[6] = x[6];
                IndexType success = dtiLMFit( G, &(io->signal[i*G]), lSi, lSiw,
                                             A, A2, (BufferType)ATA2, (BufferType)ATA3,
                                             (BufferType)xr, (BufferType)xr2, (BufferType)x2,
                                             (IndexBuffer)pivot1, (IndexBuffer)pivot2, lwork, work2,
                                             params );
                if(success==0){
                    squareDTI( (BufferType)xr, (BufferType)x );
                    x[6] = xr[6];
                }
            }
            //---------------------------------------------------------------------
            // At this point, we can copy the outputs to the proper place:
            memcpy( &(io->dti[6*i]), x, 6*sizeof(ElementType) );
            io->s0[i] = exp(scale*x[6]);
        }
    }
    while( start < args->getN() );
    
    blas_num_threads_thread(blas_threads);
    
    // Free memory previously allocated
    delete[] lSi;
    delete[] lSiw;
    delete[] A;
    delete[] A2;
    delete[] work2;
    delete[] AT;
    delete[] piA;
    delete[] wi;
    
    return (THFCNRET)NULL;
}

IndexType dtiLMFit( const SizeType G, const BufferType Si, BufferType Si2, BufferType Si3,
                    const BufferType A, BufferType A2, BufferType ATA, BufferType ATA2,
                    BufferType xr, BufferType xr2, BufferType x,
                    IndexBuffer pivot1, IndexBuffer pivot2, SizeType lwork, BufferType work,
                    const DTIParameters* params )
{
    ElementType rho = 0.001f; // Damping factor
    IndexType success = 1;
    bool computeJacobian = true;
    ElementType trace = 0.0f;
    ElementType cost;
    // -----------------------------
    // Compute the first value of the reconstructed signal
    // and the inital cost
    ElementType cost0 = 0.0f;
    // Square the current DT:
    squareDTI( xr, x );
    x[6] = xr[6]; // The logarithm of the baseline
    // Compute the current approx. to the exponential:
    mataux::multiplyMxArrays( A, x, Si2, G+1, 7, 1 );
    for( IndexType g=0; g<(IndexType)G; ++g ){
        Si2[g] = exp(Si2[g]);
        cost0 += (Si2[g]-Si[g])*(Si2[g]-Si[g]);
    }
    Si2[G] = exp(Si2[G]);
    cost0 += (Si2[G]-1.0f)*(Si2[G]-1.0f);
    // -----------------------------
    for( unsigned int it=0; it<params->maxiters; ++it ){ // Levenberg-Marquardt iterations
        // ---------------
        // 1- Compute the Jacobian if necessary
        if( computeJacobian ){ // Otherwise, A2 already contains the Jacobian and ATA J^TJ
            memcpy( A2, A, (G+1)*7*sizeof(ElementType) );
            squaredDTIDerivatives( xr, ATA ); // ATA has size 7*7, but here we use just 6*6
            // I can use just the first 6 columns of A due to the
            // way memory is laid out:
            mataux::multiplyMxArrays( A, ATA, A2, G+1, 6, 6 );
            // Compute the final form of the Jacobian and the current cost
            for( IndexType g=0; g<(IndexType)G+1; ++g ){
                for( unsigned int c=0; c<7; ++c )
                    A2[c*(G+1)+g] *= Si2[g];
            }
            mataux::transposeMultiplyMxArray( A2, ATA, G+1, 7 );
            trace = 0.0f;
            for( unsigned int c=0; c<7; ++c )
                trace += ATA[c*7+c];
        }
        // ---------------
        // 2- Try to pseudo-invert the jacobian matrix
        memcpy( ATA2, ATA, 7*7*sizeof(ElementType) );        
        mataux::addIdentityMxArray( ATA2, trace*rho, 7 );
        if( mataux::checkAndInvertMxArray( ATA2, 7, params->rcondth, pivot1, pivot2, lwork, work ) != 0 ){
            // The matrix is not invertible. We have to increase the damping
            // factor and try again
            rho = 2.0f*rho;
            computeJacobian = false;
            continue;
        }
        // ---------------
        // 3- Compute the step and update:
        for( IndexType g=0; g<(IndexType)G; ++g )
            Si3[g] = Si[g] - Si2[g];
        Si3[G] = 1.0f - Si2[G];
        mataux::multiplyMxArrays( Si3, A2, x, 1, G+1, 7 );
        mataux::multiplyMxArrays( ATA2, x, xr2, 7, 7, 1 );
        for( unsigned int c=0; c<7; ++c )
            xr2[c] = xr[c] + xr2[c];
        // ---------------
        // 4- Check if the cost actually decreased, otherwise
        //    increase the damping factor and continue
        squareDTI( xr2, x );
        x[6] = xr2[6];
        mataux::multiplyMxArrays( A, x, Si3, G+1, 7, 1 );
        cost = 0.0f;
        for( IndexType g=0; g<(IndexType)G; ++g ){
            Si3[g] = exp(Si3[g]);
            cost += (Si3[g]-Si[g])*(Si3[g]-Si[g]);
        }
        Si3[G] = exp(Si3[G]);
        cost += (Si3[G]-1.0f)*(Si3[G]-1.0f);
        if( cost>cost0 ){
            rho = 2.0f*rho;
            computeJacobian = false;
            continue;
        }
        // ---------------
        // 5- The iteration was successful, so that we can update
        memcpy( Si2, Si3, (G+1)*sizeof(ElementType) );
        cost0 = cost;
        rho = rho/2.0f;
        computeJacobian = true;
        // ---------------
        // 6- It only remains to check the stopping condition
        ElementType diff = 0.0f;
        ElementType norm = 1000*mxGetEps();
        for( unsigned int c=0; c<7; ++c ){
            diff += (xr[c]-xr2[c])*(xr[c]-xr2[c]);
            norm += xr[c]*xr[c];
        }
        diff = sqrt(diff/norm);
        memcpy( xr, xr2, 7*sizeof(ElementType) );
        if( diff<=params->tol ){
            success = 0;
            break;
        }
        // ---------------
    }
    return success;
}

void squareDTI( const BufferType xr, BufferType x )
{
    x[0] = xr[0]*xr[0] + xr[1]*xr[1] + xr[2]*xr[2];
    x[1] = xr[0]*xr[1] + xr[1]*xr[3] + xr[2]*xr[4];
    x[2] = xr[0]*xr[2] + xr[1]*xr[4] + xr[2]*xr[5];
    x[3] = xr[1]*xr[1] + xr[3]*xr[3] + xr[4]*xr[4];
    x[4] = xr[1]*xr[2] + xr[3]*xr[4] + xr[4]*xr[5];
    x[5] = xr[2]*xr[2] + xr[4]*xr[4] + xr[5]*xr[5];
    return;
}

void DTIFromEigs( const ElementType l0, const ElementType l1, const ElementType l2, const BufferType eigvec, BufferType dti )
{
    unsigned int pos = 0;
    for( unsigned int r=0; r<3; ++r ){
        for( unsigned int c=r; c<3; ++c ){
            dti[pos]  = 0.0f;
            dti[pos] += l0 * eigvec[r]   * eigvec[c];
            dti[pos] += l1 * eigvec[3+r] * eigvec[3+c];
            dti[pos] += l2 * eigvec[6+r] * eigvec[6+c];
            ++pos;
        }
    }
    return;
}

void squaredDTIDerivatives( const BufferType xr, BufferType ATA )
{
    mataux::setValueMxArray( ATA, 6, 6, 0.0f );
    // ---
    ATA[6*0+0] = 2.0f*xr[0];  // d T^2(0) / d T(0)
    ATA[6*0+1] = xr[1];       // d T^2(1) / d T(0)
    ATA[6*0+2] = xr[2];       // d T^2(2) / d T(0)
    // All others are ALWAYS 0
    // ---
    ATA[6*1+0] = 2.0f*xr[1];  // d T^2(0) / d T(1)
    ATA[6*1+1] = xr[0]+xr[3]; // d T^2(1) / d T(1)
    ATA[6*1+2] = xr[4];       // d T^2(2) / d T(1)
    ATA[6*1+3] = 2.0f*xr[1];  // d T^2(3) / d T(1)
    ATA[6*1+4] = xr[2];       // d T^2(4) / d T(1)
    // All others are ALWAYS 0
    // ---
    ATA[6*2+0] = 2.0f*xr[2];  // d T^2(0) / d T(2)
    ATA[6*2+1] = xr[4];       // d T^2(1) / d T(2)
    ATA[6*2+2] = xr[0]+xr[5]; // d T^2(2) / d T(2)
    ATA[6*2+4] = xr[1];       // d T^2(4) / d T(2)
    ATA[6*2+5] = 2.0f*xr[2];  // d T^2(5) / d T(2)
    // All others are ALWAYS 0
    // ---
    ATA[6*3+1] = xr[1];       // d T^2(1) / d T(3)
    ATA[6*3+3] = 2.0f*xr[3];  // d T^2(3) / d T(3)
    ATA[6*3+4] = xr[4];       // d T^2(4) / d T(3)
    // All others are ALWAYS 0
    // ---
    ATA[6*4+1] = xr[2];       // d T^2(1) / d T(4)
    ATA[6*4+2] = xr[1];       // d T^2(2) / d T(4)
    ATA[6*4+3] = 2.0f*xr[4];  // d T^2(3) / d T(4)
    ATA[6*4+4] = xr[3]+xr[5]; // d T^2(4) / d T(4)
    ATA[6*4+5] = 2.0f*xr[4];  // d T^2(5) / d T(4)
    // All others are ALWAYS 0
    // ---
    ATA[6*5+2] = xr[2];       // d T^2(2) / d T(5)
    ATA[6*5+4] = xr[4];       // d T^2(4) / d T(5)
    ATA[6*5+5] = 2.0f*xr[5];  // d T^2(5) / d T(5)
    // All others are ALWAYS 0
    // ---
    return;
}
