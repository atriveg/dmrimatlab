/** 
 * This header file implements the conversion from SH representations
 * of functions defined over the unit sphere to Higher Order Tensor
 * (HOT) representations, back and forth, since SH functions up to 
 * order LMAX span the same functional space as HOT of degree L. 
 * For a reference on this, see:
 * 
 *   M. Descoteaux, E. Angelino, S. Fitzgibbons, and R. Deriche. 
 *   "Apparent Diffusion Profile estimation from High Angular 
 *   Resolution Diffusion Images: estimation and applications."
 *   Magnetic Resonance in Medicine, 56(2):395–410, 2006.
 * 
 * The evaluation of HOT at desired orientations is also implemented.
 */

#ifndef _sh2hot_h_
#define _sh2hot_h_

#include "mexToMathsTypes.h"

namespace sh2hot
{
    /**
     * The conversion matrix from HOT to SH can be analytically
     * computed by integrating the SH basis functions times the
     * terms x^nx*y^ny*z^nz, with (nx+ny+nz)=LMAX. This is done
     * here using Mathematica to compute the analytical value of
     * these integrals and evaluating them with arbitrary
     * precision (see sh2hot.mathematica). These values are 
     * hard-coded within sh2hothardcodes.cxx, where these headers
     * are implemented; the first one can be used to allocate
     * a buffer of doubles to store the coefficients of the
     * conversion matrix, and the second will compute these 
     * coefficients into this buffer.
     */
    unsigned long sh2hotHardcodesDim( const unsigned int );
    void sh2hotHardcodes( const unsigned int, BufferType );
    
    /**
     * The following functions implement the evaluation of
     * HOT within arbitrary values of the unit sphere.
     */
    void computeHOTPowers( const unsigned int, unsigned int*, unsigned int*, unsigned int* );
    void computeHOTMultiplicity( const unsigned int, unsigned long*, const unsigned int*, const unsigned int*, const unsigned int* );
    void evaluateHOT( const unsigned int, const BufferType, const BufferType, const BufferType, BufferType, const unsigned int L, 
                      const unsigned int*, const unsigned int*, const unsigned int*, const unsigned long*, const BufferType );

} // end namespace sh2hot

#endif
