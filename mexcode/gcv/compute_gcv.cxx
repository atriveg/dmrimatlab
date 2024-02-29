/*==========================================================
 * compute_gcv_.c
 *
 * Implements generalized cross-validation in regularized
 * linear least squares problems. It is a means to find the
 * optimal regularization parameter as the one minimizing
 * the leave-one-out fitting error in the linear model
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2023 - Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"

#include "compute_gcv.h"
#include "../mathsmex/matrixCalculus.h"


namespace gcv
{

/** This function intiallized the default values for the
 * parameters structure used by GCV and allocates intermediate
 * buffers used by LAPACK as needed */
void allocateGCVMemory( const SizeType neqs, const SizeType NR, GCVParams* params)
{
    params->lambda0 = 0.1f;
    params->lambdastep = exp(-log(10)/5.0f);
    params->maxvals = 20;
    params->Pinv  = new ElementType[NR*NR];
    params->Pinv0 = new ElementType[NR*NR];
    params->R = new ElementType[neqs*NR];
    params->S = new ElementType[neqs*neqs];
    params->e = new ElementType[neqs];
    return;
}

void freeGCVMemory( GCVParams* params)
{
    delete[] params->Pinv;
    delete[] params->Pinv0;
    delete[] params->R;
    delete[] params->S;
    delete[] params->e;
    return;
}

/**
 * This the function does the job itself. The I/O parameters are as
 * follows:
 *
 *    Phi: a buffer with the "A matrix" of the linear problem, i.e. the
 *         dictionary where the signal is represented. It is assumed to
 *         have dimensions neqs x NR
 *    Phisq: a buffer containing the value of Phi' * Phi, hence assuming
 *         size NR x NR
 *    Lsq: a buffer containing L' * L, with L the linear operator defining
 *         the regularization. It has size NR x NR
 *    y: a buffer containing the vector of values to fit, with size neqs x 1
 *    Pinv: a buffer, sized NR x NR, with the final value of (Phi'*Phi+lambda*L'*L)^(-1)
 *         corresponding to the optimal value found for lambda
 *    lambda: the optimal value found for lambda
 *    cost: the leave-one-out cost at the optimal value of lambda
 *    neqs: the first dimension of Phi (number of equations)
 *    NR: the second dimension of Phi (number of unknowns)
 *    params: a structure containing algorithm parameters and auxiliar
 *         buffers used by LAPACK. Must be allocated with a call to
 *         "allocateGCVMemory(...)"
 *
 * Returns 0 upon success, nonzero otherwise
 * */
int computeGCV( const BufferType Phi,
                const BufferType Phisq,
                const BufferType Lsq,
                const BufferType y,
                BufferType Pinv,
                double& lambda,
                double& cost,
                const SizeType neqs,
                const SizeType NR,
                GCVParams* params
               )
{
    double lambda0 = params->lambda0;
    double cost0   = mxGetInf();
    lambda = lambda0;
    cost   = 0.0f;
    int result = 1;
    for( unsigned int i=0; i<params->maxvals; ++i ){
        //-------------------------------------------------------------------
        // Compute first the matrix to invert:
        //          Pinv <- (Phi'*Phi + lambda*L'*L)
        memcpy( Pinv, Lsq, NR*NR*sizeof(ElementType) );
        mataux::scalaropMxArray( Pinv, NR, NR, lambda, mataux::MULTIPLY );
        mataux::addMxArrays( Phisq, Pinv, Pinv, NR, NR );
        //-------------------------------------------------------------------
        // Invert the matrix using LAPACK
        ptrdiff_t NR_  = NR;
        ptrdiff_t info = 0;
        // Prepare right hand side (params->Pinv) by constructing
        // an identity matrix
        setIdentityMatrix( NR, params->Pinv );
        // Call LAPACK's implmentation
        dposv( "L", &NR_, &NR_, Pinv, &NR_, params->Pinv, &NR_, &info );
        // If dposv worked appropriately, params->Pinv contains the
        // inverse of the original Pinv matrix, i.e.
        //      params->Pinv <- (Phi'*Phi + lambda*L'*L)^(-1)
        // The value stored in Pinv is no longer valid (dpsov overwrites it)
        if(info!=0)
            return -1;
        //-------------------------------------------------------------------
        // Create the big matrix used to test the GCV cost
        mataux::multiplyMxArrays( Phi, params->Pinv, params->R, neqs, NR, NR ); // neqs x NR
        mataux::multiplyMxArraysTranspose( params->R, Phi, params->S, neqs, NR, neqs ); // neqs x neqs
        for( IndexType r=0; r<neqs; ++r ) // Subtract the identity matrix; neqs x neqs
            params->S[r+r*neqs] -= 1;
        //-------------------------------------------------------------------
        // Compute the GCV cost:
        mataux::multiplyMxArrays( params->S, y, params->e, neqs, neqs, 1 ); // neqs x 1
        ElementType trace = mataux::traceMxArray( params->S, neqs );
        cost = 0.0f;
        for( IndexType r=0; r<neqs; ++r )
            cost += (params->e[r])*(params->e[r]);
        cost /= (trace*trace);
        //-------------------------------------------------------------------
        // If the cost has increased since the last value was tried,
        // we must exit:
        if( cost>cost0 ){
            lambda = lambda0;
            cost   = cost0;
            result = 0;
            break;
        }
        else{
            lambda0 = lambda;
            cost0   = cost;
            memcpy( params->Pinv0, params->Pinv, NR*NR*sizeof(ElementType) );
            lambda = lambda*(params->lambdastep);
        }
    }
    // If the loop yields always decreasing costs and exits
    // upon reaching the maximum number of iterations,
    // the current value of lambda has been unproperly updated
    // as: lambda = lambda*(params->lambdastep), which
    // corresponds to a value of lambda that has not been
    // actually tested. We should revert this:
    if(result>0)
        lambda = lambda0;
    // In any case, the best value for the pseudo-inverse is
    // stored in params->Pinv0
    memcpy( Pinv, params->Pinv0, NR*NR*sizeof(ElementType) );
    return result;
}

/**
 * This function builds an identity matrix inside buffer,
 * which is assumed to be allocated with size NR x NR
 * */
void setIdentityMatrix( const SizeType NR, BufferType buffer )
{
    for( IndexType r=0; r<NR; ++r ){
        for( IndexType c=0; c<NR; ++c )
            buffer[r+c*NR] = 0.0f;
        buffer[r+r*NR] = 1;
    }
}

} // end namespace gcv
