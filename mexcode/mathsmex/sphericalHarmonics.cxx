#ifndef _sphericalHarmonics_cxx
#define _sphericalHarmonics_cxx

#include "sphericalHarmonics.h"
#include "math.h"

namespace shmaths
{

/** Cumulative product. It returns factorial(end)/factorial(start) or, alternatively,
 (start+1)·(start+2)·...·(end-1)·end. Note that necessarily end>start; otherwise, it
 returns simply 1. This function can be used only with SMALL unsigned integer numbers
 since it uses integer arithmetic. Otherwise, it could produce overflow.
 
 Last version: 07/05/2010
 */
double cumprod( const unsigned long start, const unsigned long end )
{
   if( start>=end )
      return 1.0;
   else
      return( cumprod(start,end-1) * (double)end );
}

void cumprod( const unsigned long N, double* values )
{
    values[0] = 1.0;
    for( unsigned long n=1; n<=N; ++n )
        values[n] = values[n-1]*n;
}

/** Compute the value of the associated Legendre polynomials of degrees l = 0, 1, ... L
 at the desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which should
    be initialized with a call of the form:
    
    float* buffer = allocateBufferForAssociatedLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later
 on by means of the delete[] operator. Note that for each degree l (2l+1) polynomials,
 corresponding to orders m = -l, ..., 0, ... l have to be evaluated. Hence, the total
 amount of polynomials to be evaluated grow to (L+1)^2.
 
 To compute the desired values for positive m, we sequentially call the function
 computeAssociatedLegendrePolynomialsL() below.
 
 For negative values of m, we can use the relation:
 
 P_l^{-m}(x) = (-1)^m·(l-m)!/(l+m)!·P_l^m(x)
 
 The function cumprod implemented above can be used for this computation.
 
 Last version: 20/04/2022
*/

void computeAssociatedLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
    unsigned int l0 = 0;
    for( unsigned int l=0; l<=L; ++l ){
        // l0 points, for each order l, at the (l,m=0) polynomial
        // Hence, we can fill all m>=0 values with a call to
        // computeAssociatedLegendrePolynomialsL():
        computeAssociatedLegendrePolynomialsL( x, l, &buffer[l0] );
        // Now, fill the values for m<0 in terms of the respective
        // values for m>0
        double sign=-1.0f;
        for( unsigned int m=1; m<=l; ++m ){
            buffer[(int)l0-(int)m] = sign*buffer[l0+m]/cumprod((int)l-(int)m,l+m);
            sign = -sign;
        }
        // Update l0 to point to (l+1,m=0):
        l0 = (l+1)*(l+2);
    }
}

/**
 * This is a robust computation of the associated Legendre polynomials of order L for
 * m form 0 to L, which is subsequentially called by computeAssociatedLegendrePolynomials
 * to compute all values up to order L.
 * It is based on matlab's implementation, which is in turned based on a fortran implementation
 * Last version: 20/04/2022
 */
void computeAssociatedLegendrePolynomialsL( const double x, const unsigned int L, double* buffer )
{
    if(L==0){ // The trivial case
        buffer[0] = 1.0f;
        return;
    }
    double s = sqrt(1.0-x*x);
    // Check different cases for robust computations
    if( s>0.0f ){
        double twocot = -2.0*x/s;
        double sn = 1.0f;
        for( unsigned int l=0; l<L; ++l )
            sn *= (-s);
        double tol = sqrt(__realmin__);
        if(abs(sn)<=tol){ // Underflow
            // Approx. solution to x*ln(x) = y
            double v  = 9.2-log(tol)/(L*s); // Large, positive number (always >9.2)
            double w  = 1.0f/log(v);        // Positive number
            double m1 = 1 + (L*s)*v*w*(1.0058+ w*(3.819 - w*12.173));
            unsigned int mm1 = ( floor(m1)<L ? (unsigned int)(floor(m1)) : L ); // mm1 = min(L, floor(m1));
            mm1 = ( mm1>=1 ? mm1 : 1 ); // mm1 = max( mm1, 1 );
            // Values of m greater than m1 will drop down to 0:
            for(unsigned int l=m1; l<=L; ++l )
                buffer[l] = 0.0f;
            // Determine the proper sign for the recursion:
            if( x>=0.0f )
                buffer[mm1-1] = __eps__*(double)( mm1==2*(mm1/2) ? -1.0f : 1.0f );
            else
                buffer[mm1-1] = __eps__*(double)( L==2*(L/2) ? 1.0f : -1.0f );
            // Apply recursion:
            double sumsq = tol;
            for( int m=(int)mm1-2; m>=0; --m ){
                buffer[m]  = ( buffer[m+1]*twocot*(m+1) - buffer[m+2]*sqrt((double)((int)L+m+2)*(double)((int)L-m-1)) );
                buffer[m] /= sqrt((double)((int)L+m+1)*(double)((int)L-m));
                sumsq += ( buffer[m] * buffer[m] );
            }
            double scale = 1.0f/sqrt( 2.0f*sumsq - buffer[0]*buffer[0] );
            for( unsigned int m=0; m<=L; ++m )
                buffer[m] *= scale;
        }
        else{ // No underflow
            // Normalization constant for m=n:
            double c = 1.0f;
            for( unsigned int d=1; d<=L; ++d )
                c *= ( 1.0f - 1.0f/(2.0f*d) );
            // Use sn = (-s).^n  to write the m = n function
            buffer[L]   = sqrt(c)*sn;
            buffer[L-1] = buffer[L]*twocot*L/sqrt(2.0f*L); // L is strictly greater than 0
            // Recursive rule down to m=0
            for( int m=(int)L-2; m>=0; --m ){
                buffer[m]  = ( buffer[m+1]*twocot*(m+1) - buffer[m+2]*sqrt((double)((int)L+m+2)*(double)((int)L-m-1)) );
                buffer[m] /= sqrt((double)((int)L+m+1)*(double)((int)L-m));
            }
        }
    }
    else{ // Polar argument, x=+-1
        for( unsigned int l=1; l<=L; ++l )
            buffer[l] = 0.0f;
        if( 2*(L/2)==L ) // even L
            buffer[0] = 1.0f;
        else // odd L
            buffer[0] = x;
    }
    // Normalize:
    for( unsigned int l=1; l<=L; ++l )
        buffer[l] *= sqrt(cumprod((int)L-(int)l,L+l));
}

/** Compute the value of the Legendre polynomials of degrees l = 0, 1, ... L at the 
 desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which sould
 be initialized with a call of the form:
 
 float* buffer = allocateBufferForLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later 
 on by means of the delete[] operator.
 
 To compute the desired values, we used the following recursive formula:
 
   P_0(x)  = 1;
   P_1(x)  = x;
   lP_l(x) = (2l-1)xP_(l-1)(x) - (l-1)P_(l-2)(x)
 
 Last version: 07/05/2010
*/
void computeLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
   for( unsigned int l=0; l<=L; ++l ){
      if( l==0 )
         buffer[l] = 1;
      else if( l==1 )
         buffer[l] = x;
      else
         buffer[l] = ( ((double)(2*l-1))*x*buffer[l-1] - ((double)(l-1))*buffer[l-2] ) / ((double)l);
   }
}

/** Allocate the buffer necessary to store the evaluations of all associated
 Legendre polynomials of order l = 0, 1, ..., L at a given real number -1 <= x <= 1.
 For each degree l, possible orders range from m=-l to m=l, so we have to allocate
 $\sum_{l=0}^{L}(2*l+1) = (L+1)^2$ doubles.
 
 THE MEMORY ALLOCATED BY THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForAssociatedLegendrePolynomials( const unsigned int L )
{
   return (new double[(L+1)*(L+1)]);
}

/** Allocate the buffer necessary to store the evaluations of all Legendre 
 polynomials of order l = 0, 1, ..., L at a given real number -1 <= x <= 1.
 We have to allocate $\sum_{l=0}^{L}1 = (L+1)$ doubles.
 
 THE MEMORY ALLOCATE WITH THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForLegendrePolynomials( const unsigned int L )
{
   return (new double[L+1]);
}

/** Compute the value of the associated Legendre polynomials of even degrees l = 0, 2, ... L
 at the desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which sould
 be initialized with a call of the form:
 
 float* buffer = allocateBufferForEvenAssociatedLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later
 on by means of the delete[] operator. Note that for each degree l (2l+1) polynomials,
 corresponding to orders m = -l, ..., 0, ... l have to be evaluated. Hence, the total
 amount of polynomials to be evaluated grow to (L+1)·(L/2+1) (only even degrees are
 stored). L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED TO THIS RESPECT.
 
 This function is based on computeAssociatedLegendrePolynomials; in fact, the values for all
 degrees (both even and odd) are computed, and only those corresponding to even orders are
 stored in the final buffer. This is because we need a recursive formula to compute the value
 of the polynomials, involving even and odd degree polyomials. Note that this function is
 provided only for convenience, since it is very useful in problems showing spherical symmetry.
 
 Last version: 07/05/2010
 */
void computeEvenAssociatedLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
   // Compute all the associated Legendre polynomials, both for even and odd orders:
   double* aux = allocateBufferForAssociatedLegendrePolynomials( L );
   computeAssociatedLegendrePolynomials( x, L, aux );
   // We have to keep only the polynomials associated to even orders:
   for( unsigned int l=0; l<=L/2; ++l ){
      // The position of degree 2l, m=0 in the actual indexing of the auxiliar buffer (auxiliar l index):
      unsigned int ali = (2*l)*(2*l) + (2*l);
      // The position of degree 2l, m=0 in the actual indexing of the final buffer (final buffer l index):
      unsigned int fli = 2*l*l + l;
      // Now, we can copy from the auxiliar buffer to the final buffer for each order m:
      for( int m=-(int)(2*l); m<=(int)(2*l); ++m )
         buffer[(int)fli+m] = aux[(int)ali+m];
   }
   // Delete the buffer with all the associated polynomials:
   delete[] aux;
}

/**
 * This function does exactly the same as the previous one, but it uses a externally maintained buffer "aux" for
 * internal computations. This may be useful in case this function has to be repeatedly called from outisde.
 * NOTE: aux should have size (L+1)*(L+1), as allocated by allocateBufferForAssociatedLegendrePolynomials(L)
 */
void computeEvenAssociatedLegendrePolynomials( const double x, const unsigned int L, double* buffer, double* aux )
{
   // Compute all the associated Legendre polynomials, both for even and odd orders:
   computeAssociatedLegendrePolynomials( x, L, aux );
   // We have to keep only the polynomials associated to even orders:
   for( unsigned int l=0; l<=L/2; ++l ){
      // The position of degree 2l, m=0 in the actual indexing of the auxiliar buffer (auxiliar l index):
      unsigned int ali = (2*l)*(2*l) + (2*l);
      // The position of degree 2l, m=0 in the actual indexing of the final buffer (final buffer l index):
      unsigned int fli = 2*l*l + l;
      // Now, we can copy from the auxiliar buffer to the final buffer for each order m:
      for( int m=-(int)(2*l); m<=(int)(2*l); ++m )
         buffer[(int)fli+m] = aux[(int)ali+m];
   }
}

/** Compute the value of the Legendre polynomials of even degrees l = 0, 2, 4, ... L
 at the desired point x. Such point must satisfy the condition -1 <= x <= 1, and NO CHECKING
 IS PERFORMED TO THIS RESPECT. The values computed are stored in "buffer" which sould
 be initialized with a call of the form:
 
 float* buffer = allocateBufferForEvenLegendrePolynomials( L );
 
 for the same value of L as used in this function. This buffer has to be erased later 
 on by means of the delete[] operator. L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED
 TO THIS RESPECT.
 
 This function is based on computeLegendrePolynomials; in fact, the values for all degrees
 (both even and odd) are computed, and only those corresponding to even orders are stored
 in the final buffer. This is because we need a recursive formula to compute the value of the 
 polynomials, involving even and odd degree polyomials. Note that this function is provided
 only for convenience, since it is very useful in problems showing spherical symmetry.
 
 Last version: 07/05/2010
 */
   
void computeEvenLegendrePolynomials( const double x, const unsigned int L, double* buffer )
{
   // Use the general purpose function to compute the values for all orders:
   double* aux = allocateBufferForLegendrePolynomials( L );
   computeLegendrePolynomials( x, L, aux );
   // Store the values only for even degrees:
   for( unsigned int l=0; l<=L/2; ++l )
      buffer[l] = aux[2*l];
   // Delete auxiliar buffer:
   delete[] aux;
}

/** Allocate the buffer necessary to store the evaluations of all even degree associated
 Legendre polynomials of order l = 0, 2, ..., L at a given real number -1 <= x <= 1.
 For each degree l, possible orders range from m=-l to m=l, so we have to allocate
 $\sum_{l=0}^{L/2}(2*(2l)+1) = (L+1)·(L/2+1)$ doubles.
 
 L IS ASSUMED TO BE ODD, AND NO OTHER CHECKING IS DONE
 
 THE MEMORY ALLOCATE WITH THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForEvenAssociatedLegendrePolynomials( const unsigned int L )
{
   return (new double[(L+1)*(L/2+1)]);
}

/** Allocate the buffer necessary to store the evaluations of all Legendre 
 polynomials of even orders l = 0, 2, ..., L at a given real number -1 <= x <= 1.
 We have to allocate $\sum_{l=0}^{L/2}1 = (L/2+1)$ doubles.
 
 L IS ASSUMED TO BE ODD, AND NO OTHER CHECKING IS DONE
 
 THE MEMORY ALLOCATE WITH THIS FUNCTION MUST BE DELETED EXTERNALLY WITH THE
 delete[] OPERATOR
 
 Last version: 07/05/2010
 */
double* allocateBufferForEvenLegendrePolynomials( const unsigned int L )
{
   return (new double[L/2+1]);
}
   
/** Return the number of Legendre polynomials up to degree L */
unsigned int getNumberOfLegendrePolynomials( const unsigned int L )
{
   return (L+1);
}

/** Return the number of associated Legendre polynomials up degree order L */
unsigned int getNumberOfAssociatedLegendrePolynomials( const unsigned int L )
{
   return ( (L+1)*(L+1) );
}

/** Return the number of Legendre polynomials up to degree L (only even degrees).
 L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED*/
unsigned int getNumberOfEvenLegendrePolynomials( const unsigned int L )
{
   return (L/2+1);
}

/** Return the number of associated Legendre polynomials up to degree L (only even degrees).
 L IS ASSUMED TO BE EVEN, AND NO CHECKING IS PERFORMED*/
unsigned int getNumberOfEvenAssociatedLegendrePolynomials( const unsigned int L )
{
   return ( (L+1)*(L/2+1) );
}
   
/** Compute the value of Legendre polynomial of degree l at the point x, which
 is assumed to lay in the range [-1,1] (NO CHECKING IS PERFORMED). Note that this
 function internally uses computeLegendrePolynomials, so the values for all degrees
 r<l have to be computed each time the function is called. Consider using instead:
 
 double computeP_l( const unsigned int l, const double* buffer );
 
 where buffer has to be computed with the aforementioned computeLegendrePolynomials
 for the same value of x. For succesive calls with the same argument x, precomputing
 the buffer can considerably accelerate the execution.
*/
double computeP_l( const double x, const unsigned int l )
{
   double*  buffer = allocateBufferForLegendrePolynomials( l );
   computeLegendrePolynomials( x, l, buffer );
   double   value  = computeP_l( l, buffer );
   delete[] buffer;
   return (value);
}

/** Compute the value of associated Legendre polynomial of degree l ans order m
 at the point x, which is assumed to lay in the range [-1,1] (NO CHECKING IS
 PERFORMED). Note that this function internally uses
 computeAssociatedLegendrePolynomials, so the values for all degrees r<l have
 to be computed each time the function is called. Consider using instead:
 
 double computeP_l( const unsigned int l, const int m, double* buffer );

 where buffer has to be computed with the aforementioned
 computeAssociatedLegendrePolynomials for the same value of x. For succesive
 calls with the same argument x, precomputing the buffer can considerably accelerate
 the execution.
 
 Note that m is assumed to lay in the range [-l,l], and no bound checking is
 performed.
 */
double computeP_l_m( const double x, const unsigned int l, const int m )
{
   double*  buffer = allocateBufferForAssociatedLegendrePolynomials( l );
   computeAssociatedLegendrePolynomials( x, l, buffer );
   double   value = computeP_l_m( l, m, buffer );
   delete[] buffer;
   return (value);
}

/** Compute the value of Legendre polynomial of degree l at the point x, which
 is assumed to lay in the range [-1,1] (NO CHECKING IS PERFORMED). The buffer
 argument can be obtained via the computeLegendrePolynomials function.
 THE USER IS RESPONSIBLE OF EXTERNALLY ALLOCATING, COMPUTING AND DELETING
 THE BUFFER
*/
double computeP_l( const unsigned int l, const double* buffer )
{
   return (buffer[l]);
}

/** Compute the value of associated Legendre polynomial of degree l ans order m
 at the point x, which is assumed to lay in the range [-1,1] (NO CHECKING IS
 PERFORMED). The buffer argument can be obtained via the
 computeAssociatedLegendrePolynomials function. THE USER IS RESPONSIBLE OF
 EXTERNALLY ALLOCATING, COMPUTING AND DELETING THE BUFFER
 
 Note that m is assumed to lay in the range [-l,l], and no bound checking
 is performed.
 */
double computeP_l_m( const unsigned int l, const int m, const double* buffer )
{
   return (   buffer[ (int)(l*(l+1)) + m ]   );
}

/** Compute the real-valued Spherical Harmonics basis function of degree l and
 order m for a point given, in spherical coordinates, as:
 
   (r=1,theta,phi)
 
 Note that physics convention is used, so that the unit vectors on the
 canonycal cartesian basis are given by:
 
  (1,0,0) -> ( r=1, theta=pi/2, phi=0         )
  (0,1,0) -> ( r=1, theta=pi/2, phi=pi/2      )
  (0,0,1) -> ( r=1, theta=0,    phi=undefined )
 
 This is a modified, real-valued basis; instead of complex exponentials,
 cosines (for m>0) and sines (for m<0) are used.
 */
double computeY_l_m( const double phi, const double theta, const unsigned int l, const int m )
{
   if( m==0 )
      return (   sqrt( (2*l+1)/(4*PI) ) * computeP_l( cos(theta), l )   );
   else if( m>0 )
      return (   sqrt( (2*l+1)/(2*PI)/cumprod((int)l-m,(int)l+m) ) * computeP_l_m( cos(theta), l, m )   ) * sin( m*phi );
   else
      return (   sqrt( (2*l+1)/(2*PI)*cumprod((int)l+m,(int)l-m) ) * computeP_l_m( cos(theta), l, m )   ) * cos( m*phi );      
}

/** Compute the whole matrix of spherical harmonics for least squares fitting.
 The resulting matrix has size NxH. Hence, N represents the number of points over
 the unit-radius sphere for which the real SH-basis is evaluated (therefore, both
 theta and phi have to be length-N vectors). Each row of the resulting matrix sh
 represents the evaluation of the first H basis functions for the point defined
 by spherical coordinates (r=1, theta, phi) with physics convention. The first
 column represents (l=0,m=0). The second one is (l=1,m=-1); then (l=1,m=0), 
 (l=1,m=1),(l=2,m=-2), and so on (spherical harmonics are provided in the same
 order as the associated Legendre polynomials).
 
 NOTE: This function is meant to populated the buffer of a mxArray of doubles; 
 it is assumed that a memory block with the proper size is allocated and freed
 from the calling function.
 */
void computeSHMatrix( const unsigned int N, const double* theta, const double* phi, const unsigned int L, BufferType shdata )
{
   // Set the proper size for the resulting matrix:
   unsigned int K = getNumberOfAssociatedLegendrePolynomials(L);
   // Allocate the buffer to precompute the associated Legendre polynomials:
   double* buffer = allocateBufferForAssociatedLegendrePolynomials( L );
   // Now, for each point in the unit sphere (and hence, for each row of the matrix):
   for( unsigned int n=0; n<N; ++n ){
      // Get the point to evaluate the polynomials:
      double x0 = cos( theta[n] );
      // Precompute the buffer of associated Legendre polynomials:
      computeAssociatedLegendrePolynomials( x0, L, buffer );
      // Now, compute the SH values and store them in the same buffer:
      unsigned int pos=0; // Auxiliar position in the buffer
      for( unsigned int l=0; l<=L; ++l ){  // For each degree
         for( int m=-(int)l; m<0; ++m ){  // For each negative order
            buffer[pos] = sqrt( (2*l+1)/(2*PI)*cumprod((int)l+m,(int)l-m) ) * buffer[pos] * cos( m * phi[n] );
            ++pos;
         }
         // For m=0:
         buffer[pos] = sqrt( (2*l+1)/(4*PI) ) * buffer[pos];
         ++pos;
         for( int m=1; m<=(int)l; ++m ){  // For each positive order
            buffer[pos] = sqrt( (2*l+1)/(2*PI)/cumprod((int)l-m,(int)l+m) ) * buffer[pos] * sin( m * phi[n] );
            ++pos;
         }
      }
      // Now the buffer contains the whole n-th row of the resulting matrix. We can insert it directly:
      for( unsigned int k=0; k<K; ++k )
          shdata[k*N+n] = buffer[k];
   }
   // Delete the buffer:
   delete[] buffer;
}
   
/** Compute the whole matrix of spherical harmonics for least squares fitting.
 The resulting matrix has size NxH. Hence, N represents the number of points over
 the unit-radius sphere for which the real SH-basis is evaluated (therefore, both
 theta and phi have to be length-N vectors). Each row of the resulting matrix sh
 represents the evaluation of the first H basis functions for the point defined
 by spherical coordinates (r=1, theta, phi) with physics convention. The first
 column represents (l=0,m=0). The second one is (l=1,m=-1); then (l=1,m=0), 
 (l=1,m=1),(l=2,m=-2), and so on (spherical harmonics are provided in the same
 order as the associated Legendre polynomials).
 
 As opposed to computeSHMatrix, this function considers only even degrees of
 the spherical harmonics. This implies that only functions with radial symmetry
 can be represented.
 
 NOTE: This function is meant to populated the buffer of a mxArray of doubles; 
 it is assumed that a memory block with the proper size is allocated and freed
 from the calling function.
 */
void computeSHMatrixSymmetric( const unsigned int N, const double* theta, const double* phi, const unsigned int L, BufferType shdata )
{
   // Set the proper size for the resulting matrix:
   unsigned int K = getNumberOfEvenAssociatedLegendrePolynomials(L);
   // Allocate the buffer to precompute the associanted Legendre polynomials:
   double* buffer = allocateBufferForEvenAssociatedLegendrePolynomials( L );
   // Now, for each point in the unit sphere (and hence, for each row of the matrix):
   for( unsigned int n=0; n<N; ++n ){
      // Get the point to evaluate the polynomials:
      double x0 = cos( theta[n] );
      // Precompute the buffer of associated Legendre polynomials:
      computeEvenAssociatedLegendrePolynomials( x0, L, buffer );
      // Now, compute the SH values and store them in the same buffer:
      unsigned int pos=0; // Auxiliar position in the buffer
      for( unsigned int l=0; l<=L/2; ++l ){  // For each even degree
         for( int m=-(int)(2*l); m<0; ++m ){  // For each negative order
            buffer[pos] = sqrt( (4*l+1)/(2*PI)*cumprod((int)(2*l)+m,(int)(2*l)-m) ) * buffer[pos] * cos( m * phi[n] );
            ++pos;
         }
         // For m=0:
         buffer[pos] = sqrt( (4*l+1)/(4*PI) ) * buffer[pos];
         ++pos;
         for( int m=1; m<=(int)(2*l); ++m ){  // For each positive order
            buffer[pos] = sqrt( (4*l+1)/(2*PI)/cumprod((int)(2*l)-m,(int)(2*l)+m) ) * buffer[pos] * sin( m * phi[n] );
            ++pos;
         }
      }
      // Now the buffer contains the whole n-th row of the resulting matrix. We can insert it directly:
      for( unsigned int k=0; k<K; ++k )
          shdata[k*N+n] = buffer[k];
   }
   // Delete the buffer:
   delete[] buffer;
}

/**
 * This function does exactly the same as the previous one, but it uses externally maintained buffers for
 * internal computations. This may be useful in case this function has to be repeatedly called from outisde.
 * NOTE: buffer should have size (L+1)*(L+2)/2, as allocated with allocateBufferForEvenAssociatedLegendrePolynomials(L)
 *       buffer2 should have size (L+1)*(L+1), as allocated by allocateBufferForAssociatedLegendrePolynomials(L)
 */
void computeSHMatrixSymmetric( const unsigned int N, const double* theta, const double* phi, const unsigned int L, BufferType shdata, double* buffer, double* buffer2 )
{
   // Set the proper size for the resulting matrix:
   unsigned int K = getNumberOfEvenAssociatedLegendrePolynomials(L);
   // Now, for each point in the unit sphere (and hence, for each row of the matrix):
   for( unsigned int n=0; n<N; ++n ){
      // Get the point to evaluate the polynomials:
      double x0 = cos( theta[n] );
      // Precompute the buffer of associated Legendre polynomials:
      computeEvenAssociatedLegendrePolynomials( x0, L, buffer, buffer2 );
      // Now, compute the SH values and store them in the same buffer:
      unsigned int pos=0; // Auxiliar position in the buffer
      for( unsigned int l=0; l<=L/2; ++l ){  // For each even degree
         for( int m=-(int)(2*l); m<0; ++m ){  // For each negative order
            buffer[pos] = sqrt( (4*l+1)/(2*PI)*cumprod((int)(2*l)+m,(int)(2*l)-m) ) * buffer[pos] * cos( m * phi[n] );
            ++pos;
         }
         // For m=0:
         buffer[pos] = sqrt( (4*l+1)/(4*PI) ) * buffer[pos];
         ++pos;
         for( int m=1; m<=(int)(2*l); ++m ){  // For each positive order
            buffer[pos] = sqrt( (4*l+1)/(2*PI)/cumprod((int)(2*l)-m,(int)(2*l)+m) ) * buffer[pos] * sin( m * phi[n] );
            ++pos;
         }
      }
      // Now the buffer contains the whole n-th row of the resulting matrix. We can insert it directly:
      for( unsigned int k=0; k<K; ++k )
          shdata[k*N+n] = buffer[k];
   }
}
   
/** This function generates a diagonal matrix (eig) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Laplace-Beltrami operator (i.e., the part of the Laplacian depending only on
 the angular coordinates (theta,phi). The eigenvalues for the SH functions up
 to degree L are computed and stored in the diagonal of eig with the usual 
 order: (l=0,m=0), (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ...
 (l=L,m=-L), ...(l=L,m=L).
 
 Interestingly, the associated eigenvalue for each l is the same independently
 on the value of m, and it is computed as: \lambda_{l,m} = -l*(l+1).
 
 NOTE: This function is meant to populated the buffer of a mxArray of doubles; 
 it is assumed that a memory block with the proper size is allocated and freed
 from the calling function.
 */
void computeSHEigMatrix( const unsigned int L, BufferType eigdata )
{
    unsigned int K = getNumberOfAssociatedLegendrePolynomials(L);
    // Set the proper values on the diagonal of the matrix:
    unsigned int pos = 0; // Auxiliar position counter
    for( unsigned long p=0; p<K*K; ++p )
        eigdata[p] = 0.0;
    for( unsigned int l=0; l<=L; ++l ){ // For each degree
       for( int m=-(int)l; m<=(int)l; ++m, ++pos ){
          eigdata[pos*K+pos] = -(double)( l*(l+1) );
       }
    }
}

/** This function generates a diagonal matrix (eig) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Laplace-Beltrami operator (i.e., the part of the Laplacian depending only on
 the angular coordinates (theta,phi). The eigenvalues for the SH functions up
 to degree L are computed and stored in the diagonal of eig with the usual
 order: (l=0,m=0), (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ...
 (l=L,m=-L), ...(l=L,m=L).
 
 Interestingly, the associated eigenvalue for each l is the same independently
 on the value of m, and it is computed as: \lambda_{l,m} = -l*(l+1)
 
 As opposed to computeSHEigMatrix, this function is restricted to even degrees
 of the SH basis, and hence it is only usuful for fucntions defined over the
 unit sphere showing radial symmetry. L IS ASSUMED TO BE EVEN AND NO CHECKING
 IS DONE TO THIS RESPECT.
 
 NOTE: This function is meant to populated the buffer of a mxArray of doubles; 
 it is assumed that a memory block with the proper size is allocated and freed
 from the calling function.
 */
void computeSHEigMatrixSymmetric( const unsigned int L, BufferType eigdata )
{
    unsigned int K = getNumberOfEvenAssociatedLegendrePolynomials(L);
    // Set the proper values on the diagonal of the matrix:
    unsigned int pos = 0; // Auxiliar position counter
    for( unsigned long p=0; p<K*K; ++p )
        eigdata[p] = 0.0;
    for( unsigned int l=0; l<=L/2; ++l ){ // For each degree
       for( int m=-(int)(2*l); m<=(int)(2*l); ++m, ++pos ){
          eigdata[pos*K+pos] = -(double)( (2*l)*(2*l+1) );
       }
    }
}
   
/** This function generates a diagonal matrix (frt) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Funk-Radon transform (FRT) operator (SH are also eigenfunctions with respect
 to this linear operator. The eigenvalues for the SH functions up to degree L
 are computed and stored in the diagonal of eig with the usual order: (l=0,m=0),
 (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ... (l=L,m=-L), 
 ...(l=L,m=L).

 Once again, the corresponding eigenvalue in each case depends only on the
 degree l of the SH basis function and does not depend on the order m. It
 can be computed in terms of the value of Legendre polynomials at x=0:
 
   \lambda_{l,m} = 2·\pi·P_l(0)
 
 NOTE: The eigenvalue associated to odd degrees l is always 0 (the SH
 correspond to antisymmetric functions, and hence the integration in a 
 whole equator yields 0). For even degrees, an alternative expression can
 be obtained:
 
   \lambda_{l,m} = 2·\pi·(-1)^(l/2) (l-1)!!/l!!,
 
 where !! denotes de bouble factorial: (l-1)!! = (l-1)·(l-3)·...·3·1, and
 l!! = l·(l-2)·(l-4)·...·4·2, for even l. Anyway, we made use of the
 recursive computation of Legendre polynomials instead.
 
 NOTE: This function is meant to populated the buffer of a mxArray of doubles; 
 it is assumed that a memory block with the proper size is allocated and freed
 from the calling function.
 */
void computeSHFRTMatrix( const unsigned int L, BufferType frtdata )
{
    unsigned int K = getNumberOfAssociatedLegendrePolynomials(L);
    // Compute the values of the Legendre polynomials as a block:
    double* buffer = allocateBufferForLegendrePolynomials(L);
    computeLegendrePolynomials( 0.0f, L, buffer );
    // And place the corresponding values:
    for( unsigned long p=0; p<K*K; ++p )
        frtdata[p] = 0.0;
    unsigned int pos = 0; // Auxiliar position counter
    for( unsigned int l=0; l<=L; ++l ){ // For each degree
       for( int m=-(int)l; m<=(int)l; ++m, ++pos ){
          frtdata[pos*K+pos] = (2*PI) * buffer[l];
       }
    }
    delete[] buffer;
}

/** This function generates a diagonal matrix (frt) whose diagonal entries are
 the eigenvalues of the Spherical Harmonics basis functions associated to the
 Funk-Radon transform (FRT) operator (SH are also eigenfunctions with respect
 to this linear operator. The eigenvalues for the SH functions up to degree L
 are computed and stored in the diagonal of eig with the usual order: (l=0,m=0),
 (l=1,m=-1), ..., (l=1,m=1), (l=2,m=-2), ..., (l=2,m=2), ... (l=L,m=-L),
 ...(l=L,m=L).
 
 Once again, the corresponding eigenvalue in each case depends only on the
 degree l of the SH basis function and does not depend on the order m. It
 can be computed in terms of the value of Legendre polynomials at x=0:
 
   \lambda_{l,m} = 2·\pi·P_l(0)
 NOTE: The eigenvalue associated to odd degrees l is always 0 (the SH
 correspond to antisymmetric functions, and hence the integration in a 
 whole equator yields 0). For even degrees, an alternative expression can
 be obtained:
 
 \lambda_{l,m} = 2·\pi·(-1)^(l/2) (l-1)!!/l!!,
 
 where !! denotes de bouble factorial: (l-1)!! = (l-1)·(l-3)·...·3·1, and
 l!! = l·(l-2)·(l-4)·...·4·2, for even l. Anyway, we made use of the
 recursive computation of Legendre polynomials instead.
 
 As opposed to computeSHFRTMatrix, this function is restricted to even degrees
 of the SH basis, and hence it is only usuful for fucntions defined over the
 unit sphere showing radial symmetry. L IS ASSUME DTO BE EVEN AND NO CHECKING
 IS DONE TO THIS RESPECT.
 
 NOTE: This function is meant to populated the buffer of a mxArray of doubles; 
 it is assumed that a memory block with the proper size is allocated and freed
 from the calling function.
 */
void computeSHFRTMatrixSymmetric( const unsigned int L, BufferType frtdata )
{
    unsigned int K = getNumberOfEvenAssociatedLegendrePolynomials(L);
    // Compute the values of the Legendre polynomials as a block:
    double* buffer = allocateBufferForEvenLegendrePolynomials(L);
    computeEvenLegendrePolynomials( 0.0f, L, buffer );
    // And place the corresponding values:
    for( unsigned long p=0; p<K*K; ++p )
        frtdata[p] = 0.0;
    unsigned int pos = 0; // Auxiliar position counter
    for( unsigned int l=0; l<=L/2; ++l ){ // For each degree
      for( int m=-(int)(2*l); m<=(int)(2*l); ++m, ++pos ){
         frtdata[pos*K+pos] = (2*PI) * buffer[l];
      }
    }
    delete[] buffer;
}

/** The exponential integral. Ei(x) is defined as \int_x^\infty \exp(-t)/t dt,
 and hence shows a singularity at x=0. It may be defined in terms of a power
 series expansion (regular part) and the logarithm function (singularity
 at x=0). Although this power series converge as fast as that for the
 exponential function or even faster, it is not an efficient implementation.
 
 Instead, we use the implementation provided by ALGLIB under the GPL license:
 
                  http://www.alglib.net/
 
 This implementation is very similar to that in matlab. Compared to the original
 in ALGLIB, we have stablished much less conservative thresholds, since we
 need a fast performance. Even so, it is accurate enough. If more precission is
 required, the places marked with the word "Tolerance" should be revised.
 
 NOTE: It is assumed that the input, x, is positive. Otherwise, the
 logarithm of x is returned (NaN).
 
 NOTE(2): We have dropped the argument "n" in the original ALGLIB routine, 
 since we need to compute only E_1.
 */
double expint( const double x )
{
    double result;
    double r;
    double t;
    double yk;
    double xk;
    double pk;
    double pkm1;
    double pkm2;
    double qk;
    double qkm1;
    double qkm2;
    double psi;
    double z;
   int k;
    double big;
   unsigned int cont = 0;
    
    big = 1.44115188075855872e17;
    if( x<=1e-3f ){ // asymptotic approximation -- Tolerance
      // Since x is quite small, the Taylor series expansion can be computed with
      // very few coefficients (we keep the first three terms; the modulus
      // of the next term, k=4, will be less than (1e-3)^4/(4·4!) ~ 1e-14:
      return ( -::log(x) - GAMMAEULER + x - x*x/4 + x*x*x/18 );
   }
   if( x>=30.0f ) // For greater values, it is practically zero (below 3e-15) -- Tolerance
      return 0.0f;
   if( x<=1.0f ){
      psi = - GAMMAEULER - ::log(x);
      z   = -x;
        xk  = 0;
        yk  = 1;
        pk  = 0;
        result = 0.0;
        do{
         xk      = xk + 1;
         yk      = yk * z / xk;
         pk     += 1;
         result += yk/pk;
            if( ::fabs(result)>1e-9 ) // Tolerance
            t = fabs(yk/result);
            else
            t = 1;
        }
        while( t>=1e-6 && cont++<100 ); //Tolerance
        result = psi-result;
        return result;
    }
    else{
      k      = 1;
      pkm2   = 1;
        qkm2   = x;
        pkm1   = 1.0;
        qkm1   = x + 1;
        result = pkm1/qkm1;
      do{
         if( ++k%2 == 1 ){
                yk = 1;
                xk = 1 + double(k-1)/double(2);
            }
            else{
            yk = x;
            xk = double(k)/double(2);
         }
         pk = pkm1*yk + pkm2*xk;
         qk = qkm1*yk + qkm2*xk;
            if( ::fabs(qk)>1e-9 ){ // Tolerance
            r = pk/qk;
            t = fabs((result-r)/r);
            result = r;
         }
            else
            t = 1;
         pkm2 = pkm1;
            pkm1 = pk;
            qkm2 = qkm1;
            qkm1 = qk;
            if( fabs(pk)>=big ){
                pkm2 = pkm2/big;
                pkm1 = pkm1/big;
                qkm2 = qkm2/big;
                qkm1 = qkm1/big;
            }
        }
        while( t>=1e-6 && cont++<100 ); //Tolerance
        result = result*exp(-x);
    }
    return result;
}
   
/** The non-singular exponential integral, defined as the power series
 in the definition of the exponential integral without the logarithmic
 singularity and the Euler constant. Although the convergence of this 
 power series expansion is quite fast, it is preferable to use the
 expint routine (for rasonable values) instead.
 
 The relative error in the computation of the Ein function is at most
 1.5e-4 for x=6. The absolut error is nearly 3.6e-4 for x=6.
 
 This errors can be decreased (or incrased, with the advantage
 of faster performance) modifying the thresholds used to decied wether
 to use expint or not.
 */
double Ein( const double x )
{
   if( x<0.1f ){
      // x is very small, so the Taylor series expansion with 3 terms is the most efficient
      // implementation; note that the modulus of the first neglected term (k=4) is at most
      // (1e-1)^4/(4·4!) ~ 1e-6. This means that the maximum error is, at most, in the order
      // of 1e-5.
      return( -x + x*x/4 - x*x*x/18 );
   }
   else if( x>6.0f ){
      // For larger values of x, the expint function is negligible compared to the logarithm,
      // so it makes no sense to compute its value; the relative error with this approximation
      // is at most in the order of 5e-4
      return( -::log(x) - GAMMAEULER );
   }
   else // The normal case. The computation is based on the expint function
      return ( -expint(x) - ::log(x) - GAMMAEULER );
}
   
/** This is quite a specialized function to compute the angular coordinates theta and phi of
 a spherical coordinates system from the cartesian coordinates (x,y,z). Physics convention is
 used, so that the 'z' axis correspond to theta=0, the 'x' axis to (theta=0, phi=0), and the
 'y' axis to (theta=0, phi=PI/2).
 
 Regardless on the actual 'r' coordinate, the vector [x,y,z] is normalized so that it has
 norm 1. This because this function is designed for functions defined over the unit sphere
 (which are those we are able to represent by means of SH expansions).
 
 Besides, this function is specifically designed to work with functions showing radial 
 symmetry (only EVEN degrees of the SH basis functions/associated Legendre polynomials).
 For convenience, all points are translated to the hemisphere defined by 0 <= phi < pi.
 Each point lying in the opposite hemisphere is projected onto its radially symmetric
 point.
 
 Last version: 20/04/2022
 */
void computeSphericalCoordsFromCartesian( const double x, const double y, const double z, double& theta, double& phi )
{
   //===============================================================================================================
   // First of all, eliminate the radial coordinate and work with normalized coordinates:
   double r  = ::sqrt( x*x + y*y + z*z );
   if( r < __eps__ ){ // The vector has norm zero
      theta = phi = 0.0f;
      return;
   }
   //===============================================================================================================
   // The computation of the theta coordinate is quite simple, since it only depends on z0:
   theta     = ::acos( z/r );
   // If the theta coordinate corresponds to the 'z' or '-z' axis, the phi coordinate is undefined:
   if( (::abs(x)<r*__eps__) && (::abs(y)<r*__eps__) )
       phi = 0.0f; // Avoid domain errors in calls to atan2
   else
       phi = ::atan2( y, x );
   //===============================================================================================================
   // The last step is to project the points in the hemisphere y<0 to points in the hemisphere y>=0. The '-z'
   // axis is also projected onto its radially symmetric position, 'z'.
   if( (PI-theta)<__eps__ ) // This is the '-z' axis
      theta = phi = 0.0f;
   else if( phi<0.0f ){   // Wrong semi-plane. Reflect...
      theta = PI - theta;
      phi   = phi + PI; // phi-PI+2*PI
   }
   else if( (PI-phi)<__eps__ ){ // '-x'--'z' plane
      theta = PI - theta;
      phi   = phi - PI;
   }
}
   
/** This function sequentially uses:
 
 ComputeShericalCoordsFromCartesian( const double x, const double y, const double z, double& theta, double& phi )
 
 to obtain the angular spherical coordinates theta and phi of a set of points given in Cartesian coordinates
 by the vectors x, y, and z with length N (physics conventions). See the documentation above for details.
 */
void computeSphericalCoordsFromCartesian( const double* x, const double* y, const double* z, double* theta, double* phi, const unsigned int N )
{
   for( unsigned int n=0; n<N; ++n )
      computeSphericalCoordsFromCartesian( x[n], y[n], z[n], theta[n], phi[n] );
}

/** This is an implementation of Wigner's 3j symbols, related to Clebsh Gordan symbols in quantum mechanics. 
 * Basically, they serve to analytically determine the integral of the triple product of 3 SH complex basis
 * functions (l1,m1), (l2,m2), (l3,m3), which is non-null in very few cases (the nonnull flag can be used to
 * determine it). Based on:
 * 
 *    Weisstein, Eric W. "Wigner 3j-Symbol." From MathWorld--A Wolfram Web Resource
 *    http://mathworld.wolfram.com/Wigner3j-Symbol.html
 * 
 * we use the Racah formula, as described in:
 * 
 *     David Terr (2022). Wigner3j.m
 *     https://www.mathworks.com/matlabcentral/fileexchange/5275-wigner3j-m)
 *     MATLAB Central File Exchange. Retrieved February 18, 2022.
 * 
 * NOTE: the factorials are precomputed (outside) and passed in the double* buffer passed as the last input
 * to avoid repeated (and very slow) calls to an ad-hoc function. If the buffer points to a NULL, then the 
 * function just checks if the corresponding Wigner symbol is nonnull.
 */
double computeWigner3j( const unsigned int l1, const unsigned int l2, const unsigned int l3, const int m1, const int m2, const int m3, bool& nonnull, bool& resnan, double* factorials )
{
    resnan  = false;   
    nonnull = true;
    if( (m1>(int)l1) || (m1<-(int)l1) || (m2>(int)l2) || (m2<-(int)l2) || (m3>(int)l3) || (m3<-(int)l3) ){
        resnan = true;
        return NAN;
    }
    if( m1+m2+m3 != 0 ){
        nonnull = false;
        return 0.0f;
    }
    unsigned int al12 = l1+l2;
    unsigned int sl12 = ( l1>l2 ? l1-l2 : l2-l1 );
    if( (l3<sl12) || (l3>al12) ){
        nonnull = false;
        return 0.0f;
    }
    /** NOTE: we do not check the integer perimeter rule, l1+l2+l3 integer,
     * since this function is already specialized for integer arguments */
    if(factorials==(double*)NULL)
        return 0.0f; // No need to make computations if we just need to know if this coefficient is nonnull
    // THIRD PART OF THE FORMULA: the combined factorials
    double result  = factorials[ (unsigned int)((int)l1+m1) ];
    result        *= factorials[ (unsigned int)((int)l1-m1) ];
    result        *= factorials[ (unsigned int)((int)l2+m2) ];
    result        *= factorials[ (unsigned int)((int)l2-m2) ];
    result        *= factorials[ (unsigned int)((int)l3+m3) ];
    result        *= factorials[ (unsigned int)((int)l3-m3) ];
    result         = ::sqrt(result);
    // FIRST PART OF THE FORMULA: the sign
    if(   !isEven( (int)l1 - (int)l2 - m3 )   )
        result = -result;
    // FOURTH PART OF THE FORMULA: the sum in "t":
    int t1  = (int)l2 - m1 - (int)l3;
    int t2  = (int)l1 + m2 - (int)l3;
    int t3  = (int)l1 + (int)l2 - (int)l3;
    int t4  = (int)l1 - m1;
    int t5  = (int)l2 + m2;
    int t12 = ( t1>t2 ? t1 : t2 );
    int t45 = ( t4<t5 ? t4 : t5 );    
    unsigned int mint = ( t12>0 ? (unsigned int)t12 : 0 );
    unsigned int maxt = ( t3<t45 ? (unsigned int)t3 : (unsigned int)t45 );
    double sum = 0.0f;
    int    sign  = ( isEven((long)mint) ? 1 : -1 );
    for( unsigned int t=mint; t<=maxt; ++t ){
        double x  = factorials[t];
        x        *= factorials[ (unsigned int)( (int)t  - (int)t1 ) ];
        x        *= factorials[ (unsigned int)( (int)t  - (int)t2 ) ];
        x        *= factorials[ (unsigned int)( (int)t3 - (int)t  ) ];
        x        *= factorials[ (unsigned int)( (int)t4 - (int)t  ) ];
        x        *= factorials[ (unsigned int)( (int)t5 - (int)t  ) ];
        sum  += ((double)sign)/x;
        sign *= -1;
    }
    result *= sum;
    // SECOND PART OF THE FORMULA: the triangle coefficient
    double triangle  = factorials[ (unsigned int)(  (int)l1 + (int)l2 - (int)l3 ) ];
    triangle        *= factorials[ (unsigned int)(  (int)l1 - (int)l2 + (int)l3 ) ];
    triangle        *= factorials[ (unsigned int)( -(int)l1 + (int)l2 + (int)l3 ) ];
    triangle        /= factorials[ l1+l2+l3+1 ];
    // EXIT
    return (result*::sqrt(triangle));
}

/** Compute the integral of the triple product of 3 COMPLEX SH basis functions (l1,m1), (l2,m2), (l3,m3)
 * based on Wigner's 3j symbols. The resnan flag is set true if mi>|li|, for which the SH are not defined.
 * The nonnull flag is set true only if an actual computation has to be done (note the Wigner's 3j symbols
 * are very sparse, so that most of the {l1,m1,l2,m2,l3,m3} combinations lead to null output values.
 * 
 * NOTE: a double* buffer with the factorials required by computeWigner3j must be externally allocated
 *       and calculated. The function cumprod can be used to this end. If a NULL is passed, this function
 *       will simply check if the corresponding integral is nonnull.
 */
double computeTripleComplexSHProd( const unsigned int l1, const unsigned int l2, const unsigned int l3, const int m1, const int m2, const int m3, bool& nonnull, bool& resnan, double* factorials )
{
    nonnull = true;
    resnan  = false;
    bool nonnull1, nonnull2, resnan1, resnan2;
    double res1 = computeWigner3j(l1,l2,l3,m1,m2,m3,nonnull1,resnan1,factorials);
    double res2 = computeWigner3j(l1,l2,l3,0,0,0,nonnull2,resnan2,factorials);
    if( resnan1 || resnan2 ){
        resnan = true;
        return NAN;
    }
    if( !nonnull1 || !nonnull2 ){
        nonnull = false;
        return 0.0f;
    }
    double res  = res1*res2;
    res        *= ::sqrt( (double)(2*l1+1)*(double)(2*l2+1)*(double)(2*l3+1)/ (4*PI) );
    return res;
}

/** Compute the integral of the triple product of 3 REAL SH basis functions (l1,m1), (l2,m2), (l3,m3)
 * based on Wigner's 3j symbols. The resnan flag is set true if mi>|li|, for which the SH are not defined.
 * The nonnull flag is set true only if an actual computation has to be done (note the Wigner's 3j symbols
 * are very sparse, so that most of the {l1,m1,l2,m2,l3,m3} combinations lead to null output values.
 * 
 * NOTE: a double* buffer with the factorials required by computeWigner3j must be externally allocated
 *       and calculated. The function cumprod can be used to this end. If a NULL is passed, this function
 *       will simply check if the corresponding integral is nonnull.
 */
double computeTripleSHProd( const unsigned int l1, const unsigned int l2, const unsigned int l3, const int m1, const int m2, const int m3, bool& nonnull, bool& resnan, double* factorials )
{
    // In the real SH basis, those functions where mi<=0 are chosen as the (scaled) real part of the corresponding
    // complex function, meanwhile mi>0 corresponds to the (scaled) imaginary part. The real part has
    // cos(mi*phi), and the imaginary part sin(mi*phi). For the integral in phi to be nonnull, it is necessary that
    // the integrand is an even function in phi, meaning either all mi are <=0 or only one mi is <=0, i.e.
    // we need m1*m2*m3 <=0
    nonnull = true;
    resnan  = false;
    if( (m1>(int)l1) || (m1<-(int)l1) || (m2>(int)l2) || (m2<-(int)l2) || (m3>(int)l3) || (m3<-(int)l3) ){
        resnan = true;
        return NAN;
    }
    // Any of the following cases means the integrand in phi is an odd function,
    // and there is no need to perform any actual computation:
    if(   (m1*m2*m3>0)
        || ( (m1==0) && (m2*m3<0) ) || ( (m2==0) && (m1*m3<0) ) || ( (m3==0) && (m1*m2<0) )
        || ( (m1==0) && (m2==0) && (m3>0) ) || ( (m1==0) && (m3==0) && (m2>0) ) || ( (m2==0) && (m3==0) && (m1>0) )   )
    {
        nonnull = false;
        return 0.0f;
    }
    double res = 0.0f; // Returned value
    int wm1 = m1; // +-m1, to express the sine/cosine in Euler's notation
    int wm2 = m2; // +-m2, to express the sine/cosine in Euler's notation
    int wm3 = m3; // +-m3, to express the sine/cosine in Euler's notation
    int sg1[2] = {1,1}; sg1[1] = ( m1<=0 ? 1 : -1 ); // Real part: z+conj(z), imag part: z-conj(z)
    int sg2[2] = {1,1}; sg2[1] = ( m2<=0 ? 1 : -1 ); // Real part: z+conj(z), imag part: z-conj(z)
    int sg3[2] = {1,1}; sg3[1] = ( m3<=0 ? 1 : -1 ); // Real part: z+conj(z), imag part: z-conj(z)
    // Now, take into account that for SH: conj(Y_l^m) = (-1)^m Y_l^{-m}
    sg1[1] = ( m1==2*(m1/2) ? sg1[1] : -sg1[1] );
    sg2[1] = ( m2==2*(m2/2) ? sg2[1] : -sg2[1] );
    sg3[1] = ( m3==2*(m3/2) ? sg3[1] : -sg3[1] );
    double w1 = ( m1==0 ? 0.5f : 0.5f*::sqrt(2) ); //  if m1=0, the real SH is not scaled; otherwise, it is scaled a factor *sqrt(2)
    double w2 = ( m2==0 ? 0.5f : 0.5f*::sqrt(2) ); //  if m2=0, the real SH is not scaled; otherwise, it is scaled a factor *sqrt(2)
    double w3 = ( m3==0 ? 0.5f : 0.5f*::sqrt(2) ); //  if m3=0, the real SH is not scaled; otherwise, it is scaled a factor *sqrt(2)
    double weight = w1*w2*w3;
    nonnull = false; // At least one of the 8 addends has to be non null
    for( unsigned int i1=0; i1<2; ++i1 ){
        for( unsigned int i2=0; i2<2; ++i2 ){
            for( unsigned int i3=0; i3<2; ++i3 ){
                bool tmpnn, tmprn;
                // ----------------
                double tmp = (double)(sg1[i1]*sg2[i2]*sg3[i3])
                    *computeTripleComplexSHProd( l1, l2, l3, wm1, wm2, wm3, tmpnn, tmprn, factorials );
                // ----------------
                resnan  = ( resnan | tmprn );  // nan if any addend is nan
                nonnull = ( nonnull | tmpnn ); // non null if any addend is non null
                // ----------------
                // The integrand is not an odd function only if either
                // all mi are negative, so that it is the product of three
                // even functions, or there are two positive mi and one
                // negative mi, so that we have even*even*odd=even. In the
                // latter case, we have the product of two imaginary parts,
                // so that we have to place a minus sign to equal the j*j
                // term in the denominator:
                if( (m1>0) || (m2>0) || (m3>0) )
                    res -= weight*tmp; // weight stands for the computation of real/imag parts and real SH scaling
                else
                    res += weight*tmp; // weight stands for the computation of real/imag parts and real SH scaling
                // ----------------
                wm3 = -wm3;
            }
            wm2 = -wm2;
        }
        wm1 = -wm1;
    }
    if(resnan){
        nonnull = true;
        return NAN;
    }
    return res;
}

/** For a given non-negative, even SH order L, compute the number of nonnull
 * triple real SH products, i.e. the number of integrals:
 * 
 *    int_{S} Y_l1^m1 Y_l2^m2 Y_l3^m3 dS
 * 
 *    for l1=0,2,...,L, m1=-l1..l1; l2=0,2,...,L, m2=-l2..l2; l3=0,2,...L,L+1,...,2L, m3=-l3..l3;
 * 
 * that have a value different from 0. Note the set of all possible indices {(l1,m1),(l2,m2),(l3,m3)}
 * has cardinal (L+1)(L+2)/2 x (L+1)(L+2)/2 x (2L+1)(2L+2)/2, but only a small subset of them yields
 * to nonnull integrals
 * 
 * NOTE: It is not internally checked if L is actually even
 */
SizeType computeNumberOfNonnullWignerSymbols( const unsigned int L )
{
    bool nonnull,resnan;
    SizeType pos = 0;
    for( unsigned int l3=0; l3<=2*L; l3+=2 ){
        for( int m3=-(int)l3; m3<=(int)l3; ++m3 ){
            for( unsigned int l2=0; l2<=L; l2+=2 ){
                for( int m2=-(int)l2; m2<=(int)l2; ++m2 ){
                    for( unsigned int l1=0; l1<=L; l1+=2 ){
                        for( int m1=-(int)l1; m1<=(int)l1; ++m1 ){
                            computeWigner3j(l1,l2,l3,m1,m2,m3,nonnull,resnan,(double*)NULL);
                            if(nonnull)
                                pos++;
                        }
                    }
                }
            }
        }
    }
    return pos;
}

/** For a given non-negative, even SH order L, compute all possible values of
 * triple real SH products, i.e. the values of the integrals:
 * 
 *    int_{S} Y_l1^m1 Y_l2^m2 Y_l3^m3 dS
 * 
 *    for l1=0,2,...,L, m1=-l1..l1; l2=0,2,...,L, m2=-l2..l2; l3=0,2,...L,L+1,...,2L, m3=-l3..l3;
 * 
 * that have a value different from 0. Note the set of all possible indices {(l1,m1),(l2,m2),(l3,m3)}
 * has cardinal (L+1)(L+2)/2 x (L+1)(L+2)/2 x (2L+1)(2L+2)/2, but only a small subset of them yields
 * to nonnull integrals
 * 
 * Note the indices for which these values are obtained are also returned.
 * 
 * NOTE: This function will not allocate memory for any of the buffers; the calling function must do
 *       so. Use computeNumberOfNonnullWignerSymbols() to get the proper size for these buffers
 * NOTE: It is not internally checked if L is actually even
 */
void computeNonnullWignerSymbols( const unsigned int L, BufferType wigner, unsigned int* ll1, unsigned int* ll2, unsigned int* ll3, int* mm1, int* mm2, int* mm3 )
{
    unsigned long NF = 4*L+1;
    double* factorials = new double[NF+1];
    cumprod( NF, factorials );
    bool nonnull,resnan;
    SizeType pos = 0;
    for( unsigned int l3=0; l3<=2*L; l3+=2 ){
        for( int m3=-(int)l3; m3<=(int)l3; ++m3 ){
            for( unsigned int l2=0; l2<=L; l2+=2 ){
                for( int m2=-(int)l2; m2<=(int)l2; ++m2 ){
                    for( unsigned int l1=0; l1<=L; l1+=2 ){
                        for( int m1=-(int)l1; m1<=(int)l1; ++m1 ){
                            double res = computeWigner3j(l1,l2,l3,m1,m2,m3,nonnull,resnan,factorials);
                            if(nonnull){
                                wigner[pos] = res;
                                ll1[pos]    = l1;
                                ll2[pos]    = l2;
                                ll3[pos]    = l3;
                                mm1[pos]    = m1;
                                mm2[pos]    = m2;
                                mm3[pos]    = m3;
                                pos++;
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] factorials;
}

/** When computing the SH coefficients of the square of a real, antipodal symmetric function
 * from its SH coefficients, each of the former can be written in terms of a symmetric quadratic
 * form on the latter. These quadratic forms are vastly sparse, so that there is no need to 
 * compute a whole 3-D array with all its entries, but instead a sparse matrix. This function 
 * computes the number of nonzeros in such matrix. Use this function to allocate memory before 
 * calling computeSquaredSHFactors.
 */
SizeType computeNumberOfSquaredSHFactors( const unsigned int L )
{
    bool nonnull,resnan;
    SizeType pos = 0;
    for( unsigned int l3=0; l3<=2*L; l3+=2 ){
        for( int m3=-(int)l3; m3<=(int)l3; ++m3 ){
            for( unsigned int l2=0; l2<=L; l2+=2 ){
                for( int m2=-(int)l2; m2<=(int)l2; ++m2 ){
                    for( unsigned int l1=0; l1<=L; l1+=2 ){
                        for( int m1=-(int)l1; m1<=(int)l1; ++m1 ){
                            double res = computeTripleSHProd(l1,l2,l3,m1,m2,m3,nonnull,resnan,(double*)NULL);
                            if(nonnull)
                                pos++;
                        }
                    }
                }
            }
        }
    }
    return pos;
}


/** When computing the SH coefficients of the square of a real, antipodal symmetric function
 * from its SH coefficients, each of the former can be written in terms of a symmetric quadratic
 * form on the latter. These quadratic forms are vastly sparse, so that there is no need to 
 * compute a whole 3-D array with all its entries, but instead a sparse matrix. This function 
 * computes the nonzero entries in such matrix. NOTE this function will neither allocate nor
 * free any memory, and it is the calling function which must do so. Use
 * computeNumberOfSquaredSHFactors() to compute the number of ElmentTypes to allocate.
 * 
 */
void computeSquaredSHFactors( const unsigned int L, BufferType factors, unsigned int* ll1, unsigned int* ll2, unsigned int* ll3, int* mm1, int* mm2, int* mm3 )
{
    unsigned long NF = 4*L+1;
    double* factorials = new double[NF+1];
    cumprod( NF, factorials );
    bool nonnull,resnan;
    SizeType pos = 0;
    for( unsigned int l3=0; l3<=2*L; l3+=2 ){
        for( int m3=-(int)l3; m3<=(int)l3; ++m3 ){
            for( unsigned int l2=0; l2<=L; l2+=2 ){
                for( int m2=-(int)l2; m2<=(int)l2; ++m2 ){
                    for( unsigned int l1=0; l1<=L; l1+=2 ){
                        for( int m1=-(int)l1; m1<=(int)l1; ++m1 ){
                            double res = computeTripleSHProd(l1,l2,l3,m1,m2,m3,nonnull,resnan,factorials);
                            if(nonnull){
                                factors[pos] = res;
                                ll1[pos]     = l1;
                                ll2[pos]     = l2;
                                ll3[pos]     = l3;
                                mm1[pos]     = m1;
                                mm2[pos]     = m2;
                                mm3[pos]     = m3;
                                pos++;
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] factorials;
}

/** This function translates two buffers of double indexed SH functions, as a duple (l,m), into
 * a single indexed buffer corresponding to the order of corresponding SH functions in a stack
 * of unrolled functions (0,0), (2,-2), (2,-1), ..., (2,2), (4,-4), ..., (4,4) will translate
 * into 0, 1, 2, ..., 5, 6, ...., 14.
 * NOTE: No sanity checks are performed on either l[n] or m[n]
 */
void unrollEvenSHIndices( const SizeType N, const unsigned int* l, const int* m, IndexBuffer k )
{
    for( SizeType n=0; n<N; ++n )
        k[n] = ((IndexType)l[n]+1)*(IndexType)l[n]/2 + (IndexType)m[n]; 
    return;
}

/** This is the reverse translation of that implemented in unrollEvenSHIndices: takes one input
 * buffer, which comprises ordinals in the even, real SH basis, and outputs two buffers with the
 * corresponding l and m indices.
 */
void rollEvenSHIndices( const SizeType N, const IndexBuffer k, unsigned int* l, int* m )
{
    for( SizeType n=0; n<N; ++n ){
        double ld = ( 1 + ::sqrt((double)(1+8*k[n])) )/4 + 100*mxGetEps();
        l[n] = 2*(unsigned int)ld;
        m[n] = (int)( k[n] - ((IndexType)l[n]+1)*(IndexType)l[n]/2 );
    }
    return;
}

bool isEven( const int L )
{
    return ( L==2*(L/2) );
}

} // End namespace shmaths

#endif // #ifndef _sphericalHarmonics_cxx
