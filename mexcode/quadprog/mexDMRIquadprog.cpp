/*==========================================================
 * mexDMRIquadprog.cpp
 *
 * This is a test for dmriquadprog
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2022 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "../mathsmex/mexToMathsTypes.h"
#include "dmriquadprog.h"

int cholesky_invert( double* Q, double* Qi, const unsigned int N );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // ========================================================================================
    /**
     * Inputs:
     *    prhs[0]: Q, N x N
     *    prhs[1]: f, N x 1
     *    prhs[2]: A, CI x N
     *    prhs[3]: b, CI x 1
     *    prhs[4]: Aeq, CE x N
     *    prhs[5]: beq, CE x 1
     *    prhs[6]: lb, N x 1
     *    prhs[7]: ub, N x 1
     * Optional:
     *    prhs[8]: maxiters, 1 x 1
     *    prhs[9]: step, 1 x 1
     * Outputs:
     *    plhs[0]: x, N x 1
     *    plhs[1]: status, 1 x 1
     *    plhs[2]: iters, 1 x 1
     *    plhs[3]: cost, 1 x 1
     *    plhs[4]: lambda, (CI+CE+CL+CU) x 1
     */
    // ========================================================================================
    if( nlhs==0 )
        return;
    // ========================================================================================
    /* check for proper number of arguments */
    if(nrhs<8) 
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:nrhs","At least 8 inputs required");
    // ========================================================================================
    const unsigned int N = mxGetM(prhs[0]);
    if( mxGetN(prhs[0])!=N )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 0 must be square");
    // ========================================================================================
    if( (mxGetM(prhs[1])!=N) || (mxGetN(prhs[1])!=1) )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 1 must be N x 1");
    // ========================================================================================
    unsigned int CI = 0;
    if( mxGetNumberOfElements(prhs[2])>0 ){
        CI = mxGetM(prhs[2]);
        if( mxGetN(prhs[2])!=N )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 2 must be CI x N");
    }
    // ========================================================================================
    if( mxGetNumberOfElements(prhs[3])>0 ){
        if( (mxGetM(prhs[3])!=CI) || (mxGetN(prhs[3])!=1) )
            mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 3 must be CI x 1");
    }
    else if( CI>0 )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 3 must be CI x 1, but it is empty");
    // ========================================================================================
    unsigned int CE = 0;
    if( mxGetNumberOfElements(prhs[4])>0 ){
        CE = mxGetM(prhs[4]);
        if( mxGetN(prhs[4])!=N )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 4 must be CE x N");
    }
    // ========================================================================================
    if( mxGetNumberOfElements(prhs[5])>0 ){
        if( (mxGetM(prhs[5])!=CE) || (mxGetN(prhs[5])!=1) )
            mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 5 must be CE x 1");
    }
    else if( CE>0 )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 5 must be CE x 1, but it is empty");
    // ========================================================================================
    unsigned int CL = 0;
    if( mxGetNumberOfElements(prhs[6])>0 ){
        CL = mxGetM(prhs[6]);
        if( (mxGetM(prhs[6])!=N) || (mxGetN(prhs[6])!=1) )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 6 must be N x 1");
    }
    // ========================================================================================
    unsigned int CU = 0;
    if( mxGetNumberOfElements(prhs[7])>0 ){
        CU = mxGetM(prhs[7]);
        if( (mxGetM(prhs[7])!=N) || (mxGetN(prhs[7])!=1) )
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:args","Input 7 must be N x 1");
    }
    // ========================================================================================
    unsigned int maxiters = 100;
    if(nrhs>8)
        maxiters = (unsigned int)(mxGetDoubles(prhs[8])[0]);
    // ========================================================================================
    double step = 10.0;
    if(nrhs>9)
        step = mxGetDoubles(prhs[9])[0];
    // ========================================================================================
    dmriqpp::QPProblem problem;
    dmriqpp::QPAuxiliar auxiliar;
    dmriqpp::QPParams params;
    params.maxiters = maxiters;
    
    dmriqpp::allocateQPProblem( problem, N, CE, CI, CL, CU );
    
    problem.step = step;
    
    dmriqpp::allocateQPAuxiliar( problem, params, auxiliar );
    
    dmriqpp::setf( problem, mxGetDoubles(prhs[1]) );
    if(CI>0){
        dmriqpp::setA( problem, mxGetDoubles(prhs[2]) );
        dmriqpp::setb( problem, mxGetDoubles(prhs[3]) );
    }
    if(CE>0){
        dmriqpp::setAeq( problem, mxGetDoubles(prhs[4]) );
        dmriqpp::setbeq( problem, mxGetDoubles(prhs[5]) );
    }
    if(CL>0){
        dmriqpp::setlb( problem, mxGetDoubles(prhs[6]) );
    }
    if(CU>0){
        dmriqpp::setub( problem, mxGetDoubles(prhs[7]) );
    }
    // ========================================================================================
    double* Q = new double[N*N];
    memcpy( Q, mxGetDoubles(prhs[0]), N*N*sizeof(double) );
    int info = cholesky_invert( Q, problem.Qi, N );
    delete[] Q;
    if(info!=0){
        dmriqpp::freeQPAuxiliar( auxiliar );
        dmriqpp::freeQPProblem( problem );
        mexErrMsgIdAndTxt("MyToolbox:mexDMRIquadprog:invert","Matrix is not Cholesky-invertible");
    }
    // ========================================================================================
    int result = solveQuadraticProgram( problem, auxiliar, params );
    dmriqpp::freeQPAuxiliar( auxiliar );
    double cost = problem.cost;
    unsigned int iters = problem.iters;
    // ========================================================================================
    plhs[0] = mxCreateDoubleMatrix( N, 1, mxREAL );
    memcpy( mxGetDoubles(plhs[0]), problem.x, N*sizeof(double) );
    // ========================================================================================
    if( nlhs>1 ){
        plhs[1] = mxCreateDoubleMatrix( 1, 1, mxREAL );
        mxGetDoubles(plhs[1])[0] = (double)result;
    }
    // ========================================================================================
    if( nlhs>2 ){
        plhs[2] = mxCreateDoubleMatrix( 1, 1, mxREAL );
        mxGetDoubles(plhs[2])[0] = (double)iters;
    }
    // ========================================================================================
    if( nlhs>3 ){
        plhs[3] = mxCreateDoubleMatrix( 1, 1, mxREAL );
        mxGetDoubles(plhs[3])[0] = (double)cost;
    }
    // ========================================================================================
    if( nlhs>4 ){
        plhs[4] = mxCreateDoubleMatrix( CI+CE+CL+CU, 1, mxREAL );
        double* lambda = mxGetDoubles(plhs[4]);
        unsigned int pos = 0;
        if(CI>0){
            memcpy( (&lambda[pos]), problem.l, CI*sizeof(double) );
            pos += CI;
        }
        if(CE>0){
            memcpy( (&lambda[pos]), problem.leq, CE*sizeof(double) );
            pos += CE;
        }
        if(CL>0){
            memcpy( (&lambda[pos]), problem.mu, CL*sizeof(double) );
            pos += CL;
        }
        if(CU>0){
            memcpy( (&lambda[pos]), problem.eta, CU*sizeof(double) );
            pos += CU;
        }
    }
    // ========================================================================================
    dmriqpp::freeQPProblem( problem );
    return;
}

int cholesky_invert( double* Q, double* Qi, const unsigned int N )
{
    ptrdiff_t nrhs = N;
    ptrdiff_t NR = N;
    ptrdiff_t info = 0;
    for( unsigned long p=0; p<N*N; ++p )
        Qi[p] = 0.0;
    for( unsigned int n=0; n<N; ++n )
        Qi[n+N*n] = 1.0;
    
    dposv( "L", &NR, &nrhs, Q, &NR, Qi, &NR, &info );
    
    return (int)info;
}
