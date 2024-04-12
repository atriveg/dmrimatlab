/**
 * This header file implements the logic to compute and evaluate the so-called
 * "physical" Hermite polynomials (as opposed to the "probabilistic" Hermite
 * polynomials. The advantage of the former is that the Hermite fucntions
 * computed from then, exp(-x^2/2)H_n(x) are eigen-functions for the Fourier
 * transform, i.e. they are their own Fourier transforms (up to the product
 * with the corresponding eigen-value.
 * 
 *    The first few Hermite polynomials are:
 * 
 *        H0(x) = 1
 *        H1(x) = 2x
 *        H2(x) = 4x^2 - 2
 *        H3(x) = 8x^3 - 12x
 *        H4(x) = 16x^4 - 48x^2 + 12
 * 
 * and they satisfy a simple recurrence formula:
 * 
 *        H{n+1}(x) = 2xH{n}(x) - 2nH{n-1}(x)
 * 
 * Finally, these polynomials form an orthonormal basis w.r.t. the weight
 * function exp(-x^2):
 * 
 *        int( H{m}(x) * H{n} * exp(-x^2) dx, x=-infty..infty ) = n! 2^n sqrt(pi) delta[n-m]
 */

#ifndef _hermitePols_h_
#define _hermitePols_h_

namespace hermpols
{

    unsigned long polsCoeffsSize( const unsigned int maxdeg );

    void allocatePolynomialCoeffs( double*& coeffs, const unsigned int maxdeg );
    
    void destroyPolynomialCoeffs( double*& coeffs );
    
    void allocatePolynomialValues( double*& values, const unsigned int maxdeg, const unsigned long N );
    
    void destroyPolynomialValues( double*& values );
    
    void computePols( double* coeffs, const unsigned int maxdeg );

    void evaluatePol( const double* coeffs, const double* x, double* f, const unsigned int deg, const unsigned long N );

    void evaluateAllPols(const double* coeffs, const double* x, double* f, const unsigned int maxdeg, const unsigned long N);

    void printPols( const double* coeffs, const unsigned int maxdeg );

} // namespace hermpols

#endif // #ifndef _hermitePols_h_
