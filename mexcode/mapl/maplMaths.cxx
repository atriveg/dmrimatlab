#ifndef _maplMaths_cxx_
#define _maplMaths_cxx_

#include "maplMaths.h"
#include "hermitePols.h"
#include <cmath>
#include <cstring> // memcpy
#include <stdio.h>
#include <iostream>

namespace mapl
{
    
    /**
     * Implements eq. (29) in Ozarslan's paper
     * 
     * Compute the number of basis functions used by MAPL if a
     * maximum order maxL is used.
     * The expansion uses even orders L = 0, 2, ..., maxL, since
     * these orders will result in real-valued, antipodal symmetric
     * basis functions phi(x,y,z). At each L, one has to consider
     * all possible indices n1, n2, n3 >= 0 such as n1+n2+n3=L.
     * It is easy to check that there are exactly (L+1)(L+2)/2
     * such combinations, so that the total number of basis
     * functions sum up to:
     *   sum( (L+1)(L+2)/2, L=0 (L even) ... maxL )
     *      = (L+2)(L+4)(2*L+3)/24
     * NOTE: It is assumed (the function will not check it) that
     * maxL is an even integer. Otherwise, the output of the 
     * function makes no sense.
     */
    unsigned long numBasisFunctions( const unsigned int maxL )
    {
        return ((unsigned long)maxL+2)*((unsigned long)maxL+4)*((unsigned long)maxL*2+3)/24;
    }
    
    /**
     * This is a utility function to allocate the buffers required
     * by the functions that compute the dictionaries in both the q-
     * and the R-space:
     *    x (or qx): 1st coordinate, [size: N]
     *    y (or qy): 2nd coordinate, [size: N]
     *    z (or qz): 3rd coordinate, [size: N]
     *    xPols (or qxPols): the values of all Hermite polynomials up to order maxL at all x (or qx), [size: N*(maxL+1)]
     *    yPols (or qyPols): the values of all Hermite polynomials up to order maxL at all y (or qy), [size: N*(maxL+1)]
     *    zPols (or qzPols): the values of all Hermite polynomials up to order maxL at all z (or qz), [size: N*(maxL+1)]
     *    Psi (or Phi): the dictionary, [size: N*M, where M=numBasisFunctions(maxL) is the number of dictionary atoms]
     */
    void allocateDictionaryBuffers( double*& x, double*& y, double*& z,
                          double*& xPols, double*& yPols, double*& zPols,
                          double*& Psi,
                          const unsigned long N, const unsigned int maxL )
    {
        x = new double[N];
        y = new double[N];
        z = new double[N];
        xPols = new double[N*(maxL+1)];
        yPols = new double[N*(maxL+1)];
        zPols = new double[N*(maxL+1)];
        Psi = new double[N*numBasisFunctions(maxL)];
        return;
    }
    
    /**
     * This is a utility function to destroy the memory previously allocated
     * with allocateDictionaryBuffers()
     */
    void destroyDictionaryBuffers( double*& x, double*& y, double*& z,
                          double*& xPols, double*& yPols, double*& zPols,
                          double*& Psi )
    {
        if(x!=NULL){ delete[] x; x=NULL; }
        if(y!=NULL){ delete[] y; y=NULL;  }
        if(z!=NULL){ delete[] z; z=NULL; }
        if(xPols!=NULL){ delete[] xPols; xPols=NULL;  }
        if(yPols!=NULL){ delete[] yPols; yPols=NULL; }
        if(zPols!=NULL){ delete[] zPols; zPols=NULL; }
        if(Psi!=NULL){ delete[] Psi; Psi=NULL; }
        return;
    }
    
    /**
     * This is a utility function to allocate the buffer required
     * to stored the regularization matrix, U, in Fick's MAPL
     * paper, as well as auxiliary buffers required by (and
     * described in) computeRegularizationMatrix() and
     * computeSTUmn().
     * U is real and symmetric, such that c^T*U*c is the energy
     * of the Laplacian of a signal represented by coeffcients c.
     * As such, U whould be positive defintie.
     */
     void allocateRegularizationMatrix( double*& U,
                                       double*& Snm, double*& Tnm, double*& Unm,
                                       unsigned int*& nx, unsigned int*& ny, unsigned int*& nz,
                                       const unsigned int maxL )
    {
        unsigned long nfuncs = numBasisFunctions( maxL );
        U  = new double[ (unsigned long long)nfuncs * (unsigned long long)nfuncs ];
        nx = new unsigned int[nfuncs];
        ny = new unsigned int[nfuncs];
        nz = new unsigned int[nfuncs];
        Snm = new double[(maxL+1)*(maxL+1)];
        Tnm = new double[(maxL+1)*(maxL+1)];
        Unm = new double[(maxL+1)*(maxL+1)];
        return;
    }
    
    /**
     * This is a utility function to destroy the memory previously allocated
     * with allocateDictionaryBuffers()
     */
    void destroyRegularizationMatrix( double*& U,
                                      double*& Snm, double*& Tnm, double*& Unm,
                                      unsigned int*& nx, unsigned int*& ny, unsigned int*& nz )
    {
        if( U!=NULL ){
            delete[] U;
            U = NULL;
        }
        if( nx!=NULL ){
            delete[] nx;
            nx = NULL;
        }
        if( ny!=NULL ){
            delete[] ny;
            ny = NULL;
        }
        if( nz!=NULL ){
            delete[] nz;
            nz = NULL;
        }
        if( Snm!=NULL ){
            delete[] Snm;
            Snm = NULL;
        }
        if( Tnm!=NULL ){
            delete[] Tnm;
            Tnm = NULL;
        }
        if( Unm!=NULL ){
            delete[] Unm;
            Unm = NULL;
        }
        return;
    }
    
    /**
     * Implements eqs. (10) and (22) in Ozarslan's paper
     * 
     * Compute a complete dictionary Psi relating the EAP values at the N
     * (normalized) coordinates {x,y,z} with the M = numBasisFunctions(maxL)
     * MAPL coefficients assuming antipodal symmetry, i.e:
     *    [EAP(ux*x,uy*y,uz*z)]_{N} = [Psi]_{NxM} * [coeffs]_{M}
     * NOTE: It is assumed that x, y, z have been already (externally) divided 
     * by ux, uy, uz. Accordingly, the (external) evaluations of
     * Hermite polynomials have also been done at x/ux, y/uy, z/uz
     * NOTE: xPols, yPols, and zPols will be modified in-place, and will 
     * eventually become multiplied by exp(-x*x/2), exp(-y*y/2), exp(-z*z/2).
     * 
     * All buffers are externally maintained:
     *  - Inputs to the algorthm:
     *   x     [size N]:          the normalized (and rotated) 1st coordinate in the EAP domain
     *   y     [size N]:          the normalized (and rotated) 2nd coordinate in the EAP domain
     *   z     [size N]:          the normalized (and rotated) 3rd coordinate in the EAP domain
     *  - Previously computed with hermpols::evaluateAllPols(), modified on exit so that they
     *    become multiplied by exp(-x*x/2), exp(-y*y/2), exp(-z*z/2):
     *   xPols [size N*(maxL+1)]: the values of all Hermite polynomials up to order maxL evaluated at x (rotated, normalized)
     *   yPols [size N*(maxL+1)]: the values of all Hermite polynomials up to order maxL evaluated at y (rotated, normalized)
     *   zPols [size N*(maxL+1)]: the values of all Hermite polynomials up to order maxL evaluated at z (rotated, normalized)
     *  - Output:
     *   Psi   [size N*M]:        the dictionary (note M can be computed as M = numBasisFunctions(maxL))
     * NOTE: It is assumed (the function will not check it) that maxL is an 
     * even integer. Otherwise, the computation makes no sense.
     */
    void computePsiDictionary( const double* x, const double* y, const double* z,
                               double* xPols, double* yPols, double* zPols,
                               double* Psi, const double& ux, const double& uy, const double& uz,
                               const unsigned long N, const unsigned int maxL )
    {
        // First of all, multiply Hermite polynomials by negative exponentials:
        for( unsigned int p=0; p<=maxL; ++p ){
            for( unsigned long n=0; n<N; ++n ){
                xPols[n+N*p] *= exp(-0.5*x[n]*x[n]);
                yPols[n+N*p] *= exp(-0.5*y[n]*y[n]);
                zPols[n+N*p] *= exp(-0.5*z[n]*z[n]);
            }
        }
        unsigned long pos = 0; // Global position of the basis function
        for( unsigned int L=0; L<=maxL; L+=2 ){ // Loop through even integers
            for( unsigned int n1=0; n1<=L; ++n1 ){
                for( unsigned int n2=0; n2<=L-n1; ++n2 ){
                    unsigned int n3 = L-n1-n2; // n1+n2+n3 must sum up to L
                    // At this point we have all three indices n1, n2 and n3
                    // that we can use to point to the proper Hermite
                    // polynomial. The value of pos indexes the overall position
                    // of the present dictionary atom.
                    // Each atom is the product of the three separable solutions,
                    // hence we sequentially process each of the x, y, and z
                    // coordinates. Compute first the normalization constants:
                    double normx = sqrt( pow(2.0,n1+1) * PI * tgamma(n1+1)  ) * ux;
                    double normy = sqrt( pow(2.0,n2+1) * PI * tgamma(n2+1)  ) * uy;
                    double normz = sqrt( pow(2.0,n3+1) * PI * tgamma(n3+1)  ) * uz;
                    double norm  = 1.0/(normx*normy*normz);
                    // And loop for all x,y,z:
                    for( unsigned long n=0; n<N; ++n )
                        Psi[n+pos*N] = norm * xPols[n+N*n1] * yPols[n+N*n2] * zPols[n+N*n3];
                    // We are done with this dictionary atom:
                    pos++;
                }
            }
        }
        return;
    }
    
    /**
     * Implements eqs. (4) and (23) in Ozarslan's paper
     * 
     * Compute a complete dictionary Phi relating the E(q) values at the Q
     * (normalized) coordinates {qx,qy,qz} with the M = numBasisFunctions(maxL)
     * MAPL coefficients assuming antipodal symmetry, i.e:
     *    [E(qx/(2*pi*ux),qy/(2*pi*uy),qz/(2*pi*uz))]_{Q} = [Phi]_{QxM} * [coeffs]_{M}
     * NOTE: It is assumed that qx, qy, qz have been already (externally) multiplied
     * by 2*pi*ux, 2*pi*uy, 2*pi*uz. Accordingly, the (external) evaluations of
     * Hermite polynomials have also been done at 2*pi*ux*qx, 2*pi*uy*qy, 2*pi*uz*qz
     * NOTE: xPols, yPols, and zPols will be modified in-place, and will 
     * eventually become multiplied by exp(-qx*qx/2), exp(-qy*qy/2), exp(-qz*qz/2).
     * 
     * All buffers are externally maintained:
     *  - Inputs to the algorthm:
     *   qx     [size Q]:          the normalized (and rotated) 1st coordinate in the E(q) domain
     *   qy     [size Q]:          the normalized (and rotated) 2nd coordinate in the E(q) domain
     *   qz     [size Q]:          the normalized (and rotated) 3rd coordinate in the E(q) domain
     *  - Previously computed with hermpols::evaluateAllPols(), modified on exit so that they
     *    become multiplied by exp(-qx*qx/2), exp(-qy*qy/2), exp(-qz*qz/2):
     *   qxPols [size Q*(maxL+1)]: the values of all Hermite polynomials up to order maxL evaluated at qx
     *   qyPols [size Q*(maxL+1)]: the values of all Hermite polynomials up to order maxL evaluated at qy
     *   qzPols [size Q*(maxL+1)]: the values of all Hermite polynomials up to order maxL evaluated at qz
     *  - Output:
     *   Phi   [size Q*M]:         the dictionary (note M can be computed as M = numBasisFunctions(maxL))
     * NOTE: It is assumed (the function will not check it) that maxL is an 
     * even integer. Otherwise, the computation makes no sense.
     */
    void computePhiDictionary( const double* qx, const double* qy, const double* qz,
                               double* qxPols, double* qyPols, double* qzPols,
                               double* Phi,
                               const unsigned long Q, const unsigned int maxL )
    {
        // First of all, multiply Hermite polynomials by negative exponentials:
        for( unsigned int p=0; p<=maxL; ++p ){
            for( unsigned long n=0; n<Q; ++n ){
                qxPols[n+Q*p] *= exp(-0.5*qx[n]*qx[n]);
                qyPols[n+Q*p] *= exp(-0.5*qy[n]*qy[n]);
                qzPols[n+Q*p] *= exp(-0.5*qz[n]*qz[n]);
            }
        }
        unsigned long pos = 0; // Global position of the basis function
        // In Ozarslan's paper, each of the three separable solutions
        // have a sign i^(-n), so that the overall sign will become
        // i^(-n1-n2-n3) = i^(-L) = (-1)^(-L/2) = (-1)^(L/2)
        // Since we loop through even integers to sweep L, it suffices
        // to alternate the sign each time we update L:
        int sign = 1;
        for( unsigned int L=0; L<=maxL; L+=2, sign*=(-1) ){ // Loop through even integers
            for( unsigned int n1=0; n1<=L; ++n1 ){
                for( unsigned int n2=0; n2<=L-n1; ++n2 ){
                    unsigned int n3 = L-n1-n2; // n1+n2+n3 must sum up to L
                    // At this point we have all three indices n1, n2 and n3
                    // that we can use to point to the proper Hermite
                    // polynomial. The value of pos indexes the overall position
                    // of the present dictionary atom.
                    //
                    // Each atom is the product of the three separable solutions,
                    // hence we sequentially process each of the x, y, and z
                    // coordinates.Compute first the normalization constants:
                    double normx = sqrt( pow(2.0,n1) * tgamma(n1+1)  );
                    double normy = sqrt( pow(2.0,n2) * tgamma(n2+1)  );
                    double normz = sqrt( pow(2.0,n3) * tgamma(n3+1)  );
                    double norm  = 1.0/(normx*normy*normz);
                    // And loop for all qx,qy,qz:
                    for( unsigned long q=0; q<Q; ++q )
                        Phi[q+Q*pos] = norm * sign * qxPols[q+Q*n1] * qyPols[q+Q*n2] * qzPols[q+Q*n3];
                    // We are done with this dictionary atom:
                    pos++;
                }
            }
        }
        return;
    }

    /**
     * This function does not implement any of the equations in Ozarslan's
     * paper. It creates a complete dictionary, say Psi, that relates
     * MAP-MRI coefficients to the values of the Orientation Ditribution
     * Function (ODF) at desired (unit norm) spatial directions. Note the
     * EAP is the product of three Hermite polynomials times a negative
     * quadratic exponential, hence its integration along the radial
     * coordinate reduces to the sum of many terms of the form:
     *
     *  Integral( a_i(n1,n2,n3;ux,uy,uz;dx,dy,dz)*r^(i+2)*exp(-K*r^2),
     *                                                        r=0..infty )
     *
     * i.e. a sum of weighted Gamma functions evaluated in semi-integer
     * numbers. After this function is called, we can compute the ODF
     * like: [ODF(dx,dy,dz)]_{N} = [Psi]_{NxM} * [coeffs]_{M}
     *
     * All buffers are externally maintained:
     *
     *  - Inputs to the algorithm:
     *   dx [size N]: the 1st coordinate of each unit direction
     *   dy [size N]: the 2nd coordinate of each unit direction
     *   dz [size N]: the 3rd coordinate of each unit direction
     *  - Previously computed with hermpols::computePols():
     *   pCoeffs [size (maxL+1)x(maxL+1)]: coefficients of the first (maxL+1) Hermite polynomials
     *  - To be filled inside the function, all of them allocated
     *    with hermpols::allocatePolynomialCoeffs():
     *   pCoeffsx [size (maxL+1)x(maxL+1)]: modified coeffients to compute H_{n1}(dx/ux*r)
     *   pCoeffsy [size (maxL+1)x(maxL+1)]: modified coeffients to compute H_{n2}(dy/uy*r)
     *   pCoeffsz [size (maxL+1)x(maxL+1)]: modified coeffients to compute H_{n3}(dz/uz*r)
     *  - Auxiliary buffers to compute vector convolutions (polynomial
     *    products); note that n1+n2+n3 equals, at the very most, maxL,
     *    so that we only need these buffers to have size maxL+1:
     *   pConv1 [size (maxL+1)]: coefficients of H_{n1}*H_{n2}
     *   pConv2 [size (maxL+1)]: coefficients of H_{n1}*H_{n2}*H_{n3}
     *  - Output:
     *   Psi [size N*M]: the dictionary (note M can be computed as M = numBasisFunctions(maxL))
     *
     * NOTE: It is assumed (the function will not check it) that maxL is an
     * even integer. Otherwise, the computation makes no sense.
     */
    void computeODFDictionary( const double* dx, const double* dy, const double* dz,
                               const double* pCoeffs,
                               double* pCoeffsx, double* pCoeffsy, double* pCoeffsz,
                               double* pConv1, double* pConv2,
                               double* Psi, const double& ux, const double& uy, const double& uz,
                               const unsigned long N, const unsigned int maxL, const double& contrast )
    {
        // Useful constants:
        const unsigned long pcNum = hermpols::polsCoeffsSize( maxL );
        // We will proceed direction-by-direction, so that we only have
        // to hard-copy and correct polynomial coefficients once:
        for( unsigned long n=0; n<N; ++n ){ // For each ODF direction
            // Compute the coefficients of the Hermite polynomials
            // evaluated, respectively, at dx/ux*r, dy/uy*r, dz/uz*r.
            // Within pCoeffs, the i-th column contains de coeffs
            // of the i-th polynomial in ascending order.
            memcpy( pCoeffsx, pCoeffs, pcNum*sizeof(double) );
            memcpy( pCoeffsy, pCoeffs, pcNum*sizeof(double) );
            memcpy( pCoeffsz, pCoeffs, pcNum*sizeof(double) );
            const double wx0 = dx[n]/ux;
            const double wy0 = dy[n]/uy;
            const double wz0 = dz[n]/uz;
            double wx = 1.0;
            double wy = 1.0;
            double wz = 1.0;
            for( unsigned int j=0; j<=maxL; ++j ){ // For each degree
                for( unsigned int i=j; i<=maxL; ++i ){ // For each polynomial
                    pCoeffsx[j+(maxL+1)*i] *= wx;
                    pCoeffsy[j+(maxL+1)*i] *= wy;
                    pCoeffsz[j+(maxL+1)*i] *= wz;
                }
                wx *= wx0;
                wy *= wy0;
                wz *= wz0;
            }
            // Compute the constant in the argument of the exponential,
            // exp(-K*r^2) = exp( -.5*(dx^2/ux^2+dy^2/uy^2+dz^2/uz^2)*r^2)
            const double K = 0.5*(  wx0*wx0 + wy0*wy0 + wz0*wz0 );
            // Now, loop through the basis functions:
            unsigned long pos = 0; // Global position of the basis function
            for( unsigned int L=0; L<=maxL; L+=2 ){ // Loop through even integers
                for( unsigned int n1=0; n1<=L; ++n1 ){
                    for( unsigned int n2=0; n2<=L-n1; ++n2 ){
                        unsigned int n3 = L-n1-n2; // n1+n2+n3 must sum up to L
                        // At this point we have all three indices n1, n2 and n3
                        // that we can use to point to the proper Hermite
                        // polynomials. The value of pos indexes the overall
                        // position of the present dictionary atom.
                        // Compute first the normalization constants:
                        double normx = sqrt( pow(2.0,n1+1) * PI * tgamma(n1+1)  ) * ux;
                        double normy = sqrt( pow(2.0,n2+1) * PI * tgamma(n2+1)  ) * uy;
                        double normz = sqrt( pow(2.0,n3+1) * PI * tgamma(n3+1)  ) * uz;
                        double norm  = 1.0/(normx*normy*normz);
                        // Compute the coefficients of the polynomial obtained
                        // from the product of H_{n1}(dx/ux*r)*H_{n2}(dy/uy*r)
                        // *H_{n3}(dz/uz*r), which is a polynomial of degree L
                        // in r. Note H_{ni}(di/ui*r) can be casted into a
                        // polynomial in r by multiplying the coefficients of
                        // H_{ni} by the successuve powers of di/ui, which is
                        // exactly what we did to compute pCoeffsx, pCoeffsy,
                        // and pCoeffsz. Hence, we only have to convolve the
                        // coefficients:
                        hermpols::multiplyPols( &pCoeffsx[(maxL+1)*n1],
                                               &pCoeffsy[(maxL+1)*n2], pConv1, n1, n2 );
                        hermpols::multiplyPols( pConv1, &pCoeffsz[(maxL+1)*n3],
                                               pConv2, n1+n2, n3 );
                        // Now, the L+1 first positions of pConv2 contain the
                        // coefficients of the product of the three polynomials
                        // in ascending order. Therefore, it only remains to
                        // compute all the integrals and sum up:
                        double sum  = 0.0;
                        double sqK  = sqrt(K);
                        // The basic implementation to compute the probabilistic
                        // ODF should be:
                        //    double sqK3 = sqK*sqK*sqK; // K^(3/2)
                        //    double g32  = sqrt(PI)/2;  // Gamma(3/2)
                        //    double g42  = 1.0;  // Gamma(4/2) = Gamma(2)
                        // Indeed, this can be generalized to compute
                        //    Integral( r^q*EAP(dx*r,dy*r,dz*r), r=0...infty )
                        // for any q>-1 (q=2 stands for the probabilistic OPDF
                        // as described by Tristán-Vega et al. NeuroImage, 2009;
                        // q=0 stands for the ODF as described in Tuch's Q-Balls):
                        double sqK3 = pow( sqK, contrast+1.0 );
                        double g32  = tgamma( 0.5*(contrast+1.0) );
                        double g42  = tgamma( 0.5*(contrast+2.0) );
                        bool even   = true;
                        for( unsigned int d=0; d<=L; ++d){ // pConv2 has degree n1+n2+n3=L
                            // Compute:
                            //    a_d * Integral( r^2 * r^d * exp(-K*r^2), r=0..infty )
                            //                   = sqrt(K^(d+3))/2 * Gamma((d+3)/2);
                            if(even){
                                sum += pConv2[d] * g32 / sqK3;
                                // For the probabilistic ODF:
                                //    g32 *= 0.5*(double)(d+3);
                                // But, in general:
                                g32 *= 0.5*( (double)d + contrast + 1.0 );
                            }
                            else{
                                sum += pConv2[d] * g42 / sqK3;
                                // For the probabilistic ODF:
                                //    g42 *= 0.5*(double)(d+3);
                                // But, in general:
                                g42 *= 0.5*( (double)d + contrast + 1.0 );
                            }
                            even  = !even;
                            sqK3 *= sqK;
                        } // for( unsigned int d=0; d<=L; ++d)
                        // It only remains to normalize the value with the factor
                        // 0.5 left within the integral and the normalization of
                        // the basis function itself:
                        Psi[n+pos*N] = 0.5 * norm * sum;
                        // We are done with this dictionary atom:
                        pos++;
                    } // for( unsigned int n2=0; n2<=L-n1; ++n2 )
                } // for( unsigned int n1=0; n1<=L; ++n1 )
            } // for( unsigned int L=0; L<=maxL; L+=2 )

        } // for( unsigned long n=0; n<N; ++n )
        return;
    }
    
    /**
     * This function implements eq. (10) of Fick's paper, i.e. it computes a
     * regularization matrix U for the MAPL fitting problem such that c^TUc
     * stands for the energy of the Laplacian of a signal described by coefficients
     * c. As such, U should be positive definite (it is also real, symmetric, and
     * mostly sparse).
     * * The energy of the Laplacian, obviously, depends on the three scale factors
     * ux, uy, and uz.
     * All buffers are externally maintained, and can be allocated/destroyed with
     * the respective utility functions allocateRegularizationMatrix() and
     * destroyRegularizationMatrix():
     *     U:   [size: M x M, with M = numBasisFunctions(maxL) the number of basis 
     *           functions to regularize]. The output of the function
     *     Smn: [size: maxL+1 x maxL+1], a matrix that directly implements eq. (11) 
     *           in Fick's paper, it is EXTERNALLY computed. You can compute it with
     *           a call to computeSTUmn()
     *     Tmn: [size: maxL+1 x maxL+1], a matrix that directly implements eq. (12)
     *           in Fick's paper, it is EXTERNALLY computed. You can compute it with
     *           a call to computeSTUmn()
     *     Umn: [size: maxL+1 x maxL+1], a matrix that directly implements eq. (13)
     *           in Fick's paper, it is EXTERNALLY computed. You can compute it with
     *           a call to computeSTUmn()
     *     nx:  [size: M], the order of the 1st basis function, the one in x, at
     *           each dictionary atom. It is EXTERNALLY computed with computenxyz().
     *     ny:  [size: M], the order of the 2nd basis function, the one in y, at
     *           each dictionary atom. It is EXTERNALLY computed with computenxyz().
     *     nz:  [size: M], the order of the 3rd basis function, the one in z, at
     *           each dictionary atom. It is EXTERNALLY computed with computenxyz().
     */
    void computeRegularizationMatrix( double* U,
                                      const double* Snm, const double* Tnm, const double* Unm,
                                      const unsigned int* nx, const unsigned int* ny, const unsigned int* nz,
                                      const double& ux, const double& uy, const double& uz, const unsigned int maxL )
    {
        const unsigned long Q = numBasisFunctions(maxL);
        // Now, build the elements of U one by one. Fortunately, it is symmetric,
        // so we can save nearly half of the effort:
        const unsigned int M = maxL+1;
        for( unsigned long i=0; i<Q; ++i ){ // For each row
            for( unsigned long k=i; k<Q; ++k ){ // For each column
                // Retrieve xi, yi, zi, xk, yk, zk as in eq. (10)
                unsigned int xi = nx[i];
                unsigned int yi = ny[i];
                unsigned int zi = nz[i];
                unsigned int xk = nx[k];
                unsigned int yk = ny[k];
                unsigned int zk = nz[k];
                // Compute each of the six terms and add them
                double asum = 0.0;
                asum += (ux*ux*ux)/(uy*uz) * Snm[xi+M*xk] * Unm[yi+M*yk] * Unm[zi+M*zk];
                asum += (uy*uy*uy)/(ux*uz) * Snm[yi+M*yk] * Unm[zi+M*zk] * Unm[xi+M*xk];
                asum += (uz*uz*uz)/(ux*uy) * Snm[zi+M*zk] * Unm[xi+M*xk] * Unm[yi+M*yk];
                asum += 2.0*(ux*uy)/uz * Tnm[xi+M*xk] * Tnm[yi+M*yk] * Unm[zi+M*zk];
                asum += 2.0*(uy*uz)/ux * Tnm[yi+M*yk] * Tnm[zi+M*zk] * Unm[xi+M*xk];
                asum += 2.0*(ux*uz)/uy * Tnm[xi+M*xk] * Tnm[zi+M*zk] * Unm[yi+M*yk];
                // Force symmetry (read below eq. (13) of Fick's paper):
                U[i+Q*k] = asum;
                U[k+Q*i] = asum;
            }
        }
        return;
    }
    
    /**
     * This is a utility function to compute the order nx, ny, nz of each of the
     * three separable Hermite polynomials at each of the atoms in the dictionary.
     * All buffers nx, ny, nz have size M = numBasisFunctions(maxL) and are externally
     * maintained (can be allocated with allocateRegularizationMatrix() and
     * destroyed with destroyRegularizationMatrix()). This function is intended to
     * produce the inputs requred by computeRegularizationMatrix(); since nx, ny, nz
     * are data-independent, we can do it just once
     */
    void computenxyz( unsigned int* nx, unsigned int* ny, unsigned int* nz, const unsigned int maxL )
    {
        unsigned long pos = 0;
        for( unsigned int L=0; L<=maxL; L+=2 ){ // Loop through even integers
            for( unsigned int n1=0; n1<=L; ++n1 ){
                for( unsigned int n2=0; n2<=L-n1; ++n2 ){
                    nx[pos] = n1;
                    ny[pos] = n2;
                    nz[pos] = L-n1-n2;
                    ++pos;
                }
            }
        }
        return;
    }

    /**
     * This function implements eqs. (11-13) in Fick's paper, and should be called
     * right before calling computeRegularizationMatrix(), which implements eq. (10).
     * The idea is to pre-compute eqs. (11-13) only ONCE for all the data set and
     * store the results in these matrices to avoid duplicated computations (since
     * the part of the computations done here depends only on the indices of the
     * basis functions but not on the data).
     * All buffers are externally maintained, and can be allocated/destroyed with
     * the respective utility functions allocateRegularizationMatrix() and
     * destroyRegularizationMatrix():
     *     Smn: [size: (maxL+1) x (maxL+1)], a matrix that stores the outputs of eq. (11)
     *     Tmn: [size: (maxL+1) x (maxL+1)], a matrix that stores the outputs of eq. (12)
     *     Umn: [size: (maxL+1) x (maxL+1)], a matrix that stores the outputs of eq. (13)
     */
    void computeSTUnm( double* Snm, double* Tnm, double* Unm, const unsigned int maxL )
    {
        const double pi1 = sqrt(PI);
        const double pi3 = pi1*pi1*pi1;
        const double pi7 = pi3*pi3*pi1;
        int sign = 1; // For convenience
        // These matrices are massively sparse
        unsigned long N = maxL+1;
        for( unsigned long n=0; n<N; ++n ){
            double cnst;
            unsigned long m;
            for( m=0; m<N; ++m ){ // Initialize the three buffers
                Snm[n+N*m] = 0.0;
                Tnm[n+N*m] = 0.0;
                Unm[n+N*m] = 0.0;
            }
            // Compute Snm
            cnst = 2.0*pi7*sign;
            Snm[n+N*n] = cnst * 3.0*(2*n*n+2*n+1); // The term in delta[n,m]
            if(n>=2){  // The term in delta[n,m+2]
                m = n-2;
                Snm[n+N*m] = cnst * (6+4*m)*exp(0.5*(lgamma(n+1)-lgamma(m+1)));
            }
            if(n>=4){  // The term in delta[n,m+4]
                m = n-4;
                Snm[n+N*m] = cnst * exp(0.5*(lgamma(n+1)-lgamma(m+1)));
            }
            if((long)n<(long)N-2){  // The term in delta[n+2,m]
                m = n+2;
                Snm[n+N*m] = cnst * (6+4*n)*exp(0.5*(lgamma(m+1)-lgamma(n+1)));
            }
            if((long)n<(long)N-4){  // The term in delta[n+4,m]
                m = n+4;
                Snm[n+N*m] = cnst * exp(0.5*(lgamma(m+1)-lgamma(n+1)));
            }
            // Compute Tnm
            cnst = pi3*(-sign);
            Tnm[n+N*n] = cnst * (1+2*n); // The term in delta[m,n]
            if(n>=2){ // The term in delta[n,m+2]
                m = n-2;
                Tnm[n+N*m] = cnst * sqrt( (double)n*(n-1) );
            }
            if((long)n<(long)N-2){ // The term in delta[n+2,m]
                m = n+2;
                Tnm[n+N*m] = cnst * sqrt( (double)(m)*(m-1) );
            }
            // Compute Unm
            Unm[n+N*n] = 0.5/pi1 * sign;
            // Change the sign:
            sign = -sign;
        }
        return;
    }
    
    /**
     * Implements eq. (27) in Ozarslan's paper, i.e.: given a maximum order
     * for the basis functions (the even integer maxL), computes the value
     * of all basis fucntions at qx=qy=qz=0. Note that the number of basis
     * functions for each maxL, and hence the size of the externally
     * maintained buffer Bnnn, can be computed as numBasisFunctions(maxL)
     * NOTE: It is assumed (the function will not check it) that maxL is an 
     * even integer. Otherwise, the computation makes no sense.
     */
    void computePhi0Values( double* Bnnn, const unsigned int maxL )
    {
        double l2 = log(2.0);
        unsigned long pos = 0; // Global position of the basis function
        for( unsigned int L=0; L<=maxL; L+=2 ){ // Loop through even integers
            for( unsigned int n1=0; n1<=L; ++n1 ){
                for( unsigned int n2=0; n2<=L-n1; ++n2 ){
                    unsigned int n3 = L-n1-n2; // n1+n2+n3 must sum up to L
                    bool oddity = (2*(n1/2)!=n1) || (2*(n2/2)!=n2) || (2*(n3/2)!=n3); // A tribute to David Bowie
                    if(oddity)
                        Bnnn[pos] = 0.0;
                    else{
                        // The equation to implement is:
                        //   (n1!n2!n3!)^(1/2) / (n1!!n2!!n3!!)
                        // Remember n1, n2, n3 are all even, hence
                        // ni!! = ni(ni-2)(ni-4)...2
                        //      = (ni/2)(ni/2-1)(ni/2-2)...1 * 2^(ni/2)
                        //      = (ni/2)! * 2^(ni/2)
                        Bnnn[pos]  = 0.5 * ( lgamma(n1+1) + lgamma(n2+1) + lgamma(n3+1) );
                        Bnnn[pos] -= ( lgamma(n1/2+1) + lgamma(n2/2+1) + lgamma(n3/2+1) );
                        Bnnn[pos] -= ( (n1/2)*l2 + (n2/2)*l2 + (n3/2)*l2 );
                        Bnnn[pos]  = exp(Bnnn[pos]);
                    }
                    pos++;
                }
            }
        }
        return;
    }
    
    /**
     * This is a utility function to get the four indices describing the basis function
     * (L,n1,n2,n3) from its atom position within a dictionary computed with either
     * computePsiDictionary or computePhiDictionary.
     * Note this is not particularly efficient, so it should never be repeatedly called
     * from the core algorithm
     */
    void getIndicesFromAtom( const unsigned long atom, unsigned int& L, unsigned int& n1, unsigned int& n2, unsigned int& n3 )
    {
        unsigned int maxL = 0;
        unsigned long csum=0;
        while(csum<=atom)
            csum += (maxL+=2);
        unsigned long pos = 0; // Global position of the basis function
        for( L=0; L<=maxL; L+=2 ){ // Loop through even integers
            for( n1=0; n1<=L; ++n1 ){
                for( n2=0; n2<=L-n1; ++n2 ){
                    n3 = L-n1-n2; // n1+n2+n3 must sum up to L
                    if(pos==atom)
                        return;
                    pos++;
                }
            }
        }
        return;
    }
    
    /**
     * This is a utility function to get the overall position of a basis function,
     * described by its four indices (L,n1,n2,n3), inside a dictionary computed with
     * either computePsiDictionary or computePhiDictionary from its atom number.
     * Note this is not particularly efficient, so it should never be repeatedly called
     * from the core algorithm
     * NOTE: no error checking over the indices L, n1, n2, n3 is performed, so it is
     * assumed that L is an even integer and 0<=n1<=L, 0<=n2<=L, 0<=n1+n2<=L
     * NOTE: n3 is not passed since n1+n2+n3=L
     */
    void getAtomFromIndices( unsigned long& atom, const unsigned int L, const unsigned int n1, const unsigned int n2 )
    {
        if( (n1>L) || (n2>L) || (n1+n2)>L ){ return; }
        atom = 0;
        for( unsigned int l=0; l<L; l+=2 )
            atom += (l+1)*(l+2)/2;
        for( unsigned int _n1=0; _n1<=n1; ++_n1 ){
            for( unsigned int _n2=0; _n2<=L-_n1; ++_n2 ){
                if( (n1==_n1) && (n2==_n2) )
                    return;
                atom++;
            }
        }
        return;
    }
} // namespace mapl


#endif // #ifndef _maplMaths_cxx_
