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

#ifndef _sh2hot_cxx_
#define _sh2hot_cxx_

#include "sphericalHarmonics.h"

namespace sh2hot
{
    /** This function is intended to compute all possible combinations of powers of each component of the
     vector [x,y,z]^t to be included in each term of the tensor product. For a tensor of order L, the tensor
     sumation is defined as terms of the form:
     
              p)          q)         r)             p)            q)            r)
     T_{1,1, ..., 1,2,2, ...,2,3,3, ..., 3} · x·x· ... ·x · y·y· ... ·y · z·z· ... ·z = T_{i_1,i_2,...,i_L} · x^p · y^q · z^r
     
     where p+q+r=L. Since we deal only with hypersymmetric tensors, the term T_{i_1,i_2,...,i_L} is 
     invariant under any permutation of the set of indices {i_1,i_2,...,i_L}:
     
     T_{1,1,2,2,2,3,3} = T_{3,1,2,1,2,3,2} = T_{3,3,2,2,2,1,1} = ...
     
     It is easy to prove that there are exactly (L+1)(L+2)/2 independent terms for a symmetric tensor
     of order L. The purpose of this function is to find the powers of x, y, and z for each corresponding
     term in the summation.
     
     The memory to store the indices has to be allocated externally. Note that HOT form a basis for the
     same functional space as even-order SH, so the memory required can be allocated via:
     
     unsigned int* nx = new unsigned int[ getNumberOfEvenAssociatedLegendrePolynomials(L) ];
     */
    void computeHOTPowers( const unsigned int L, unsigned int* nx, unsigned int* ny, unsigned int* nz )
    {
        // Fill the buffers. The methodology is quite simple; since we need that nx+ny+nz=l for
        // each term, we may pick up 'nx' any of 0, 1, ... l. Once 'nx' has been fixed, we may pick
        // up 'ny' any of: 0, 1, ... (l-nx). Once 'nx' and 'ny' have been fixed, 'nz' can take only
        // the value nz=l-nx-ny.
        // Note we go backwards from L to 0 with both nx and ny to get a consistent representation
        // with usual rank-2 tensors
        unsigned int pos = 0; // Auxiliar absolute position in the buffer.
        for( int x=(int)L; x>=0; --x ){
            for( int y=(int)L-x; y>=0; --y ){
                nx[pos] = (unsigned int)x;
                ny[pos] = (unsigned int)y;
                nz[pos] = (unsigned int)( (int)L - x - y );
                pos++;
            }
        }
    }
    
    /** This function is intended to compute the multiplicities of each component of a HOT. For a tensor
     of order L, the tensor sumation is defined as terms of the form:
     
              p)          q)         r)             p)            q)            r)
     T_{1,1, ..., 1,2,2, ...,2,3,3, ..., 3} · x·x· ... ·x · y·y· ... ·y · z·z· ... ·z = T_{i_1,i_2,...,i_L} · x^p · y^q · z^r
     
     where p+q+r=L. Since we deal only with hypersymmetric tensors, the term T_{i_1,i_2,...,i_L} is 
     invariant under any permutation of the set of indices {i_1,i_2,...,i_L}:
     
     T_{1,1,2,2,2,3,3} = T_{3,1,2,1,2,3,2} = T_{3,3,2,2,2,1,1} = ...
     
     It is easy to prove that there are exactly (L+1)(L+2)/2 independent terms for a symmetric tensor
     of order L, each of them appearing with a certain multiplicity in the overall sumation. The
     purpose of this method is precisely to compute such factors.
     
     The memory to store the indices has to be allocated externally. Note that HOT form a basis for the
     same functional space as even-order SH, so the memory required can be allocated via:
     
     unsigned int* mu = new unsigned int[ getNumberOfEvenAssociatedLegendrePolynomials(L) ];
     
     The arguments nx, ny, and nz can be computed via the computeHOTPowers() function.
     */
    void computeHOTMultiplicity( const unsigned int L, unsigned long* mu, const unsigned int* nx, const unsigned int* ny, const unsigned int* nz )
    {
        // Get the total number of independent terms of the tensor:
        unsigned int l = (L+1)*(L+2)/2;
        // The multiplicity of a term with powers nx, ny, and nz is easily
        // proven to be: l!/nx!/ny!/nz!. Instead of explicitly computing the
        // factorials, we use the function cumprod:
        for( unsigned int pos=0; pos<l; ++pos ){
            // Find the greater degree among nx, ny, and nz:
            if( nx[pos]>=ny[pos] && nx[pos]>=nz[pos] ){
                mu[pos] = (unsigned long)( shmaths::cumprod( nx[pos], L ) );
                for( unsigned int d=ny[pos]; d>1; --d )
                    mu[pos] /= d;
                for( unsigned int d=nz[pos]; d>1; --d )
                    mu[pos] /= d;
            }
            else if( ny[pos]>=nx[pos] && ny[pos]>=nz[pos] ){
                mu[pos] = (unsigned long)( shmaths::cumprod( ny[pos], L ) );
                for( unsigned int d=nx[pos]; d>1; --d )
                    mu[pos] /= d;
                for( unsigned int d=nz[pos]; d>1; --d )
                    mu[pos] /= d;
            }
            else{
                mu[pos] = (unsigned long)( shmaths::cumprod( nz[pos], L ) );
                for( unsigned int d=nx[pos]; d>1; --d )
                    mu[pos] /= d;
                for( unsigned int d=ny[pos]; d>1; --d )
                    mu[pos] /= d;
            }
        }
    }
    
    /**
     * This function is intended to evaluate a scalar function defined over the unit sphere given by
     * a hypersymmetric higher order tensor. An L-order tensor is evaluated (for 3-D vectors) as:
     *
     *         3      3           3
     *   sum                             A                    x    · x   · ... · x
     *       i_1=1, i_2=1, ..., i_L=1     i_1, i_2, ..., i_L    i_1   i_2         i_L
     *
     * However, since we deal only with hypersymmetric tensors, this summation of 3^L elements can be
     * reduced to a small subset of (L+1)(L+2)/2 independent terms with given muliplicities:
     *
     *         L                 p_j    q_j    r_j
     *    sum      T_j · mu_j · x    · y    · z
     *        j=1
     *
     * The multiplicities mu_j and the powers of x, y, and z for each independent term are obtained
     * via the computeHOTMultiplicity() and the computeHOTPowers() functions, respectively. 
     *
     * The meaning of the parameters are as follows:
     *
     *    N:              The number of points where the HOT has to be evaluated
     *    x, y, and z:    Each of them is a vector of length N. For a given position j=0,1,...N-1,
     *                    the HOT is evaluated for [x[j],y[j],z[j]]^T, which necessarily has
     *                    norm 1.
     *    val:            The vector with the evaluations of the HOT (output).
     *    L:              The order of the HOT (only even degrees 0, 2, 4, 6, 8, ... are considered)
     *    nx, ny, and nz: The degrees of x, y, and z, respectively, for each free component of the HOT.
     *                    Compute them with the computeHOTPowers() function.
     *    mu:             The multiplicity of each free component of the tensor. Compute it via the
     *                    computeHOTMultiplicity() function.
     *    T:              The actual independent components of the HOT, T_j.
     */
    void evaluateHOT( const unsigned int N, const BufferType x, const BufferType y, const BufferType z, 
                                     BufferType val,
                                     const unsigned int L, const unsigned int* nx, const unsigned int* ny, const unsigned int* nz,
                                     const unsigned long* mu,
                                     const BufferType T )
    {
        // Compute the number of terms of the summation depending on the order
        // of the tensor:
        unsigned int K = (L+1)*(L+2)/2;
        for( unsigned int n=0; n<N; ++n ){ // At each point of the unit sphere
            // Initiallize the output:
            val[n] = 0.0f;
            // Implement the summation:
            for( unsigned int j=0; j<K; ++j ){ // For each term
                // Initiallize the independent term, x^p·y^q·z^r:
                double cum = 1.0f;
                // Compute the whole independent term:
                for( unsigned int p=0; p<nx[j]; ++p ) // Implements x^p
                    cum *= x[n];
                for( unsigned int q=0; q<ny[j]; ++q ) // Implements y^q
                    cum *= y[n];
                for( unsigned int r=0; r<nz[j]; ++r ) // Implements z^r
                    cum *= z[n];
                // Sum the result with its corresponding multiplicity:
                val[n] += T[j] * mu[j] * cum;
            }
        }
    }

} // end namespace sh2hot
#endif
