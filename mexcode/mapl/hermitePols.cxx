#ifndef _hermitePols_cxx_
#define _hermitePols_cxx_

#include "hermitePols.h"
#include <stdio.h>
#include <iostream>

namespace hermpols
{

    /**
     * Return the size of the buffer of doubles required
     * to store all the coefficients of all Hermite polynomials
     * up to degree maxdeg
     */
    unsigned long polsCoeffsSize( const unsigned int maxdeg )
    {
        return (maxdeg+1)*(maxdeg+1);
    }
    
    /**
     * The next four are utility functions to allocate/destroy
     * memory used to comopute polynomial coefficients and values
     */
    void allocatePolynomialCoeffs( double*& coeffs, const unsigned int maxdeg )
    {
        coeffs = new double[ polsCoeffsSize(maxdeg) ];
        return;
    }
    
    void destroyPolynomialCoeffs( double*& coeffs )
    {
        if(coeffs!=NULL){
            delete[] coeffs;
            coeffs = NULL;
        }
        return;
    }
    
    void allocatePolynomialValues( double*& values, const unsigned int maxdeg, const unsigned long N )
    {
        values = new double[ N*(maxdeg+1) ];
        return;
    }
    
    void destroyPolynomialValues( double*& values )
    {
        if(values!=NULL){
            delete[] values;
            values = NULL;
        }
        return;
    }

    /**
     * Compute the coefficients of all Hermite polynomials up
     * to degree maxdeg, and store in a buffer coeffs externally
     * maintained (with size returned by polsCoeffsSize, allocated
     * with allocatePolynomialCoeffs() )
     * The i-th column contains de coefficients of the i-th polynomial
     * in ascending order. Since the i-th polynomial has degree i,
     * coeffs results in un upper-triangular matrix (i. e. all
     * elements below the diagonal are null)
     */
    void computePols( double* coeffs, const unsigned int maxdeg )
    {
        // We directly implement the recursive formula:
        //     H{n+1}(x) = 2xH{n}(x) - 2nH{n-1}(x)
        for(unsigned long k=0; k<(maxdeg+1)*(maxdeg+1); ++k )
            coeffs[k] = 0.0; // Initialize
        // First polynomial, H0(x) = 1.0:
        coeffs[0] = 1.0;
        // Second polynomial, if requested, H1(x) = 2*x
        if(maxdeg>0)
            coeffs[1+(maxdeg+1)] = 2.0;
        // For all the other polynomials, we use the two previous
        // polynomials:
        for( unsigned int n=2; n<=maxdeg; ++n ){ // Fopr each polynomial
            // Degree 0 is special, because it depends only on H{n-2}
            coeffs[(maxdeg+1)*n] = -2.0*(n-1)*coeffs[(maxdeg+1)*(n-2)];
            for( unsigned int d=1; d<=n; ++d ){ // For each degree >=1:
                coeffs[d+(maxdeg+1)*n] = 
                    2.0*coeffs[(d-1)+(maxdeg+1)*(n-1)] -  // 2*x*H{n}(x)
                    2.0*(n-1)*coeffs[d+(maxdeg+1)*(n-2)]; // 2*n*H{n-1}(x)
            }
        }
        return;
    }

    /**
     * Evaluate a polynomial with degree deg, whose coefficients in
     * ascending order are stored in coeffs, at the N abscissas stored
     * in x. The N evaluations are stored in the buffer f. All buffers
     * are externally maintained:
     *    coeffs: size deg+1
     *         x: size N
     *         f: size N
     */
    void evaluatePol( const double* coeffs, const double* x, double* f, const unsigned int deg, const unsigned long N )
    {
        // Horner's algorithm
        for( unsigned long n=0; n<N; ++n )
            f[n] = coeffs[deg];
        for( int cdeg=(int)deg-1; cdeg>=0; --cdeg ){
            for( unsigned long n=0; n<N; ++n )
                f[n] = x[n]*f[n] + coeffs[cdeg];
        }
        return;
    }

    /**
     * Evaluate all Hermite polynomials up to degree maxdeg (with coefficients
     * stored in the (maxdeg+1) x (maxdeg+1) buffer coeffs, and computed with
     * computePols) at the N abscissas stored in the vector x. The (maxdeg+1) x N
     * evaluations are stored in the output buffer f: the i-th column of f are the
     * values of the i-th Hermite polynomial at the N requested abscissas. All
     * buffers are externally maintained:
     *    coeffs: size (maxdeg+1)*(maxdeg+1), allocated with allocatePolynomialCoeffs()
     *         x: size N
     *         f: size N*(maxdeg+1), allocated with allocatePolynomialValues()
     */
    void evaluateAllPols(const double* coeffs, const double* x, double* f, const unsigned int maxdeg, const unsigned long N)
    {
        for( unsigned int p=0; p<=maxdeg; ++p )
            evaluatePol( &coeffs[p*(maxdeg+1)], x, &f[p*N], p, N );
        return;
    }

    /**
     * Handy print of the first maxdeg Hermite polynomials
     * to the standard output
     */
    void printPols( const double* coeffs, const unsigned int maxdeg )
    {
        for( unsigned int n=0; n<=maxdeg; ++n ){
            std::cout << "H_" << n << "(x) = ";
            for( int d=(int)n; d>=0; --d ){
                double coeff = coeffs[d+(maxdeg+1)*n];
                if(coeff!=0.0){
                    if(d==0){
                        if(n==0)
                            std::cout << coeff;
                        else{
                            if( coeff>0.0 )
                                std::cout << " + " << coeff;
                            else
                                std::cout << " - " << -coeff;
                        }
                    }
                    else{
                        if(coeff==1.0){
                            if(d==(int)n)
                                std::cout << "x";
                            else
                                std::cout << " + x";
                        }
                        else if(coeff==-1.0){
                            if(d==(int)n)
                                std::cout << "- x";
                            else
                                std::cout << " - x";
                        }
                        else{
                            if(d==(int)n)
                                std::cout << coeff << "·x";
                            else{
                                if(coeff>0.0)
                                    std::cout << " + " << coeff << "·x";
                                else
                                    std::cout << " - " << -coeff << "·x";
                            }
                        }
                        if(d>1)
                            std::cout << "^" << d << " ";
                        else
                            std::cout << " ";
                    }
                }
            }
            std::cout << std::endl;
        }
    }

} // namespace hermpols


#endif // #ifndef _hermitePols_cxx_
