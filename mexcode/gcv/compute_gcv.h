/** 
 * compute_gcv.h
 *
 * Headers for the GCV functions in compute_gcv.cxx
 *
 * Copyright 2023 - Antonio Tristán Vega
 */

#ifndef _compute_gcv_h_
#define _compute_gcv_h_

#include "../mathsmex/mexToMathsTypes.h"

namespace gcv
{

typedef struct GCVParams{
    double lambda0;        // The initial value of the reg. param.
    double lambdastep;     // The reg. param. is updated as lambda <- lambdastep*lambda
    unsigned int maxvals;  // Maximum number of reg. param. values to try
#ifdef OCTAVE_BUILD
    BufferType Phi_T;      // This will store Phi transposed
#endif
    BufferType Pinv0;      // This is an intermediate buffer for the psuedoinverse (previous step)
    BufferType Pinv;       // This is an intermediate buffer for the psuedoinverse (current step)
    BufferType R;          // This is Phi*( Phi'*Phi + lambda*L'*L )^(-1)RR
    BufferType S;          // This is the overall matrix, Phi*( Phi'*Phi + lambda*L'*L )^(-1)*Phi'- Id
    BufferType e;          // Buffer to store the vector with the errors
} GCVParams;

void allocateGCVMemory( const SizeType, const SizeType, GCVParams* );

void freeGCVMemory( GCVParams* );

int computeGCV( const BufferType, const BufferType, const BufferType, const BufferType,
                BufferType, double&, double&, const SizeType, const SizeType, GCVParams*);

void setIdentityMatrix( const SizeType, BufferType );

} // end namespace gcv

#endif
