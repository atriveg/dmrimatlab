/**
 * This header file implements the basic maths to work with the MAPL
 * approach, i.e. the logic to compute and handle Hermite functions,
 * defined as exp(-x^2/2)H{n}(x), with H{n}(x) thw "physics Hermite
 * polynomials". Hence, the functions here will strongly depend
 * on those described in the header file hermitePols.cxx
 */

#ifndef _maplMaths_h_
#define _maplMaths_h_

namespace mapl
{

#ifndef PI
#define PI 3.14159265358979323846
#endif
    
    unsigned long numBasisFunctions( const unsigned int maxL );
    
    void allocateDictionaryBuffers( double*& x, double*& y, double*& z,
                          double*& xPols, double*& yPols, double*& zPols,
                          double*& Psi,
                          const unsigned long N, const unsigned int maxL );
    
    void destroyDictionaryBuffers( double*& x, double*& y, double*& z,
                          double*& xPols, double*& yPols, double*& zPols,
                          double*& Psi );
    
    void allocateRegularizationMatrix( double*& U,
                                       double*& Snm, double*& Tnm, double*& Unm,
                                       unsigned int*& nx, unsigned int*& ny, unsigned int*& nz,
                                       const unsigned int maxL );
    
    void destroyRegularizationMatrix( double*& U,
                                      double*& Snm, double*& Tnm, double*& Unm,
                                      unsigned int*& nx, unsigned int*& ny, unsigned int*& nz );
    
    void computePsiDictionary( const double* x, const double* y, const double* z,
                               double* xPols, double* yPols, double* zPols,
                               double* Psi, const double& ux, const double& uy, const double& uz,
                               const unsigned long N, const unsigned int maxL );
    
    void computePhiDictionary( const double* qx, const double* qy, const double* qz,
                               double* qxPols, double* qyPols, double* qzPols,
                               double* Phi,
                               const unsigned long Q, const unsigned int maxL );

    void computeODFDictionary( const double* dx, const double* dy, const double* dz,
                               const double* pCoeffs,
                               double* pCoeffsx, double* pCoeffsy, double* pCoeffz,
                               double* pConv1, double* pConv2,
                               double* Psi, const double& ux, const double& uy, const double& uz,
                               const unsigned long N, const unsigned int maxL, const double& contrast=2.0 );
    
    void computeRegularizationMatrix( double* U,
                                      const double* Snm, const double* Tnm, const double* Unm,
                                      const unsigned int* nx, const unsigned int* ny, const unsigned int* nz,
                                      const double& ux, const double& uy, const double& uz, const unsigned int maxL );
    
    void computenxyz( unsigned int* nx, unsigned int* ny, unsigned int* nz, const unsigned int maxL );

    void computeSTUnm( double* Snm, double* Tnm, double* Unm, const unsigned int maxL );
    
    void computePhi0Values( double* Bnnn, const unsigned int maxL );
    
    void getIndicesFromAtom( const unsigned long atom, unsigned int& L, unsigned int& n1, unsigned int& n2, unsigned int& n3 );
    
    void getAtomFromIndices( unsigned long& atom, const unsigned int L, const unsigned int n1, const unsigned int n2 );


} // namespace mapl

#endif // #ifndef _maplMaths_h_
