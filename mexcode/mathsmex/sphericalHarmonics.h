/** Utilities to compute Spherical-Harmonics (SH)-related functions. Take a look to the
 corresponding implementation (.cxx) file for details on the meaning and usage of each
 of the functions in this file.
 
 IMPORTANT NOTE: Part of the implementations of the functions in this headers file is
 in the file sh2hot.cxx. Hence, it is mandatory to include sh2hot.cxx in any project
 including "sphericalHarmonics.h"
 */

#ifndef _sphericalHarmonics_h_
#define _sphericalHarmonics_h_

#include "mexToMathsTypes.h"

namespace shmaths
{
    
#ifndef PI
#define PI 3.14159265358979323846
#endif
    
#ifndef GAMMAEULER
#define GAMMAEULER 0.5772156649015328606065
#endif
    
    double cumprod( const unsigned long, const unsigned long );
    
    void cumprod( const unsigned long, double* );
    
    void computeAssociatedLegendrePolynomials( const double, const unsigned int, double* );
    
    void computeAssociatedLegendrePolynomialsL( const double, const unsigned int, double* );
    
    void computeLegendrePolynomials( const double, const unsigned int, double* );
    
    double* allocateBufferForAssociatedLegendrePolynomials( const unsigned int );
    
    double* allocateBufferForLegendrePolynomials( const unsigned int );
    
    void computeEvenAssociatedLegendrePolynomials( const double, const unsigned int, double* );

    void computeEvenAssociatedLegendrePolynomials( const double, const unsigned int, double*, double* );
    
    void computeEvenLegendrePolynomials( const double, const unsigned int, double* );
    
    double* allocateBufferForEvenAssociatedLegendrePolynomials( const unsigned int );
    
    double* allocateBufferForEvenLegendrePolynomials( const unsigned int );
    
    unsigned int getNumberOfLegendrePolynomials( const unsigned int );
    
    unsigned int getNumberOfAssociatedLegendrePolynomials( const unsigned int );
    
    unsigned int getNumberOfEvenLegendrePolynomials( const unsigned int );
    
    unsigned int getNumberOfEvenAssociatedLegendrePolynomials( const unsigned int );
    
    double computeP_l( const double, const unsigned int );
    
    double computeP_l_m( const double, const unsigned int, const int );
    
    double computeP_l( const unsigned int, const double* );
    
    double computeP_l_m( const unsigned int, const int, const double* );
    
    double computeY_l_m( const double, const double, const unsigned int, const int );
    
    void computeSHMatrix( const unsigned int, const double*, const double*, const unsigned int, BufferType );
    
    void computeSHMatrixSymmetric( const unsigned int, const double*, const double*, const unsigned int, BufferType );

    void computeSHMatrixSymmetric( const unsigned int, const double*, const double*, const unsigned int, BufferType, double*, double* );

    void computeSHEigMatrix( const unsigned int, BufferType );
    
    void computeSHEigMatrixSymmetric( const unsigned int, BufferType );
    
    void computeSHFRTMatrix( const unsigned int, BufferType );
    
    void computeSHFRTMatrixSymmetric( const unsigned int, BufferType );
        
    double expint( const double );
    
    double Ein( const double );
    
    void computeSphericalCoordsFromCartesian( const double, const double, const double, double&, double& );
    
    void computeSphericalCoordsFromCartesian( const double*, const double*, const double*, double*, double*, const unsigned int );
    
    double computeWigner3j( const unsigned int, const unsigned int, const unsigned int, const int, const int, const int, bool&, bool&, double* );
    
    double computeTripleComplexSHProd( const unsigned int, const unsigned int, const unsigned int, const int, const int, const int, bool&, bool&, double* );
    
    double computeTripleSHProd( const unsigned int, const unsigned int, const unsigned int, const int, const int, const int, bool&, bool&, double* );
    
    SizeType computeNumberOfNonnullWignerSymbols( const unsigned int );
    
    void computeNonnullWignerSymbols( const unsigned int, BufferType, unsigned int*, unsigned int*, unsigned int*, int*, int*, int* );
    
    SizeType computeNumberOfSquaredSHFactors( const unsigned int );
    
    void computeSquaredSHFactors( const unsigned int, BufferType, unsigned int*, unsigned int*, unsigned int*, int*, int*, int* );
    
    void unrollEvenSHIndices( const SizeType, const unsigned int*, const int*, IndexBuffer );
    
    void rollEvenSHIndices( const SizeType, const IndexBuffer, unsigned int*, int* );
    
    bool isEven( const int );
    
} // End namespace shmaths

#endif //_sphericalHarmonics_h_
