/**
 *
 * DSIInterp.cxx
 *
 * This is auxiliar code to perform interpolation of the EAP
 * computed in a regular Cartesian grid (i.e. the output of
 * HYDI-DSI) by using either nearest neoghbor, or linear or
 * cubic B-spline interpolation
 */

#ifndef _DSIInterp_cxx_
#define _DSIInterp_cxx_

#include "DSIInterp.h"

namespace dsiinterp{

SizeType allocate_extended_grid( BufferType& grid, const SizeType& Nx,  const SizeType& Ny,  const SizeType& Nz )
{
    SizeType M = (2*Ny+5)*(2*Nx+5)*(2*Nz+5);
    grid = new ElementType[M];
    return M;
}

void fill_extended_grid( const BufferType eap, BufferType grid, const SizeType& Nx,  const SizeType& Ny,  const SizeType& Nz )
{
    // From the "reduced" form of the EAP as returned by the estimation function,
    // fill the (antipodal symmetric) grid with its full size. Note the original
    // grid with size (2*Nx+1)(2*Ny+1)(2*Nz+1) is extended to (2*Nx+5)(2*Ny+5)(2*Nz+5)
    // by inserting 0-valued nodes at the beginning and the end of each dimension
    IndexType pos = 0;
    for( IndexType nz=0; nz<=(IndexType)Nz; ++nz ){
        // Note the weird indexing nz -> nx -> ny is consistent with
        // Matlab's meshgrid behavior, which is externally used to
        // call these functions:
        IndexType nx = ( nz==0 ? 0 : -(IndexType)Nx );
        for( ; nx<=(IndexType)Nx; ++nx ){
            IndexType ny = ( (nx==0 && nz==0) ? 0 : -(IndexType)Ny );
            for( ; ny<=(IndexType)Ny; ++ny ){
                // Main grid point:
                IndexType nxe = nx + Nx + 2;
                IndexType nye = ny + Ny + 2;
                IndexType nze = nz + Nz + 2;
                // Linear index:
                IndexType li = nye + (2*Ny+5)*nxe + (2*Ny+5)*(2*Nx+5)*nze;
                grid[li] = eap[pos];
                // Antipode:
                nxe = -nx + Nx + 2;
                nye = -ny + Ny + 2;
                nze = -nz + Nz + 2;
                // Linear index:
                li = nye + (2*Ny+5)*nxe + (2*Ny+5)*(2*Nx+5)*nze;
                grid[li] = eap[pos];
                // Next one:
                ++pos;
            }
        }
    }
    return;
}

bool interpolate_NN( const BufferType grid,
                                  const SizeType& Nx, const SizeType& Ny, const SizeType& Nz,
                                  const ElementType& x, const ElementType& y, const ElementType& z, ElementType& result )
{
    // Perform nearest negihbor interpolation within a grid with size (2Nx+5) x (2Ny+5)  x (2Nz+5)
    // The EAP is estimated in a mesh (-Nx..Nx) x (-Ny..Ny) x (-Nz..Nz). We extend each dimension
    // with a 0 value at each end: (-Nx-1..Nx+1) x (-Ny-1..Ny+1) x (-Nz-1..Nz+1). The additional point
    // along each dimension is added for consistency with B-splines interpolation:
    result = 0.0;
    if( (x<=-(ElementType)Nx-1) || (x>=(ElementType)Nx+1) || (y<=-(ElementType)Ny-1) || (y>=(ElementType)Ny+1) || (z<=-(ElementType)Nz-1) || (z>=(ElementType)Nz+1) )
        return false;
    // Otherwise, we are inbounds and we can actually interpolate. Compute
    // the associated grid coordinates for each dimension:
    IndexType xi = (IndexType)( round(x) + (Nx+2.0) ); // 1 < xi < 2*Nx+3
    IndexType yi = (IndexType)( round(y) + (Ny+2.0) ); // 1 < yi < 2*Ny+3
    IndexType zi = (IndexType)( round(z) + (Nz+2.0) ); // 1 < zi < 2*Nz+3
    // Translate to a linear index within the buffer:
    IndexType lidx  = yi + xi*(2*Ny+5) + zi*(2*Nx+5)*(2*Ny+5);
    // Interpolate:
    result = grid[lidx];
    return true;
}

bool interpolate_LN( const BufferType grid,
                                  const SizeType& Nx, const SizeType& Ny, const SizeType& Nz,
                                  const ElementType& x, const ElementType& y, const ElementType& z, ElementType& result )
{
    // Perform trilinear interpolation within a grid with size (2Nx+5) x (2Ny+5)  x (2Nz+5)
    // The EAP is estimated in a mesh (-Nx..Nx) x (-Ny..Ny) x (-Nz..Nz). We extend each dimension
    // with a 0 value at each end: (-Nx-1..Nx+1) x (-Ny-1..Ny+1) x (-Nz-1..Nz+1). The additional point
    // along each dimension is added for consistency with B-splines interpolation:
    result = 0.0;
    if( (x<=-(ElementType)Nx-1) || (x>=(ElementType)Nx+1) || (y<=-(ElementType)Ny-1) || (y>=(ElementType)Ny+1) || (z<=-(ElementType)Nz-1) || (z>=(ElementType)Nz+1) )
        return false;
    // Otherwise, we are inbounds and we can actually interpolate. Compute
    // the associated grid coordinates for each dimension:
    IndexType xi = (IndexType)( floor(x) + (Nx+2.0) ); // 1 < xi < 2*Nx+3
    IndexType yi = (IndexType)( floor(y) + (Ny+2.0) ); // 1 < yi < 2*Ny+3
    IndexType zi = (IndexType)( floor(z) + (Nz+2.0) ); // 1 < zi < 2*Nz+3
    ElementType xg = x - floor(x); // 0 <= xg <= 1
    ElementType yg = y - floor(y); // 0 <= yg <= 1
    ElementType zg = z - floor(z); // 0 <= zg <= 1
    // Trilinear interpolations will use the 2 x 2 x 2 nearest neighbors of the grid point
    // with coordinates [xi,yi,zi]. Let's compute:
    for( unsigned short k=0; k<2; ++k ){
        for( unsigned short i=0; i<2; ++i ){
            for( unsigned short j=0; j<2; ++j ){
                // The value at each grid point is multiplied by the
                // tensor product of the 4 kernels:
                // result += (grid[xi+i][yi+j][zi+k])
                //     * wx * wy * wz;
                // BUT: the grid is stored in a linear format, so we have to
                // translate the triple indexing [xi,yi,zi] to a unique linear
                // index
                IndexType lidx = (yi+j) + (xi+i)*(2*Ny+5) + (zi+k)*(2*Nx+5)*(2*Ny+5);
                ElementType wx = i*xg + (1-i)*(1.0-xg);
                ElementType wy = j*yg + (1-j)*(1.0-yg);
                ElementType wz = k*zg + (1-k)*(1.0-zg);
                result += (grid[lidx]) * wx * wy * wz;
            }
        }
    }
    return true;
}

bool interpolate_BS( const BufferType grid,
                                  const SizeType& Nx, const SizeType& Ny, const SizeType& Nz,
                                  const ElementType& x, const ElementType& y, const ElementType& z, ElementType& result )
{
    // Perform B-spline interpolation within a grid with size (2Nx+5) x (2Ny+5)  x (2Nz+5)
    // The EAP is estimated in a mesh (-Nx..Nx) x (-Ny..Ny) x (-Nz..Nz). We extend each dimension
    // with a 0 value at each end: (-Nx-1..Nx+1) x (-Ny-1..Ny+1) x (-Nz-1..Nz+1). Then, we
    // add yet and additional node at each end because B-splines require it.
    // Beyond the additional 0-valued node, the EAP is assumed to be 0:
    result = 0.0;
    if( (x<=-(ElementType)Nx-1) || (x>=(ElementType)Nx+1) || (y<=-(ElementType)Ny-1) || (y>=(ElementType)Ny+1) || (z<=-(ElementType)Nz-1) || (z>=(ElementType)Nz+1) )
        return false;
    // Otherwise, we are inbounds and we can actually interpolate. Compute
    // the associated grid coordinates for each dimension:
    IndexType xi = (IndexType)( floor(x) + (Nx+2.0) ); // 1 < xi < 2*Nx+3
    IndexType yi = (IndexType)( floor(y) + (Ny+2.0) ); // 1 < yi < 2*Ny+3
    IndexType zi = (IndexType)( floor(z) + (Nz+2.0) ); // 1 < zi < 2*Nz+3
    ElementType xg = x - floor(x); // 0 <= xg <= 1
    ElementType yg = y - floor(y); // 0 <= yg <= 1
    ElementType zg = z - floor(z); // 0 <= zg <= 1

    // Evaluate all the kernel functions for B-spline interpolation:
    ElementType Bsplines[3][4];
    splinesBi( xg, Bsplines[0][0], Bsplines[0][1], Bsplines[0][2], Bsplines[0][3] );
    splinesBi( yg, Bsplines[1][0], Bsplines[1][1], Bsplines[1][2], Bsplines[1][3] );
    splinesBi( zg, Bsplines[2][0], Bsplines[2][1], Bsplines[2][2], Bsplines[2][3] );
    // B-spline interpolations will use the 4 x 4 x 4 nearest neighbors of the grid point
    // with coordinates [xi,yi,zi]. Let's compute:
    for( unsigned short k=0; k<4; ++k ){
        for( unsigned short i=0; i<4; ++i ){
            for( unsigned short j=0; j<4; ++j ){
                // The value at each grid point is multiplied by the
                // tensor product of the 4 kernels:
                // result += (grid[xi+i-1][yi+j-1][zi+k-1])
                //     * Bsplines[0][i] * Bsplines[1][j] * Bsplines[2][k];
                // BUT: the grid is stored in a linear format, so we have to
                // translate the triple indexing [xi,yi,zi] to a unique linear
                // index
                IndexType lidx = ((IndexType)(yi+j)-1) + ((IndexType)(xi+i)-1)*(2*Ny+5) + ((IndexType)(zi+k)-1)*(2*Nx+5)*(2*Ny+5);
                result += (grid[lidx])
                     * Bsplines[0][i] * Bsplines[1][j] * Bsplines[2][k];
            }
        }
    }
    return true;
}

void splinesBi( const ElementType& x, ElementType& b0, ElementType& b1, ElementType& b2, ElementType& b3 )
{
    // Evaluate the 4 cubic B-spline functions using Horner's algotihm
    // Note x is assumed (but not checked) to be in [0,1]
    // Cubic terms:
    b0 = -1.0;
    b1 = 3.0;
    b2 = -3.0;
    b3 = 1.0;
    // Cuadratic terms:
    b0 = b0*x + 3.0;
    b1 = b1*x - 6.0;
    b2 = b2*x + 3.0;
    b3 = b3*x; // Just 0!
    // Linear terms:
    b0 = b0*x - 3.0;
    b1 = b1*x; // Just 0!
    b2 = b2*x + 3.0;
    b3 = b3*x; // Just 0!
    // Independent terms:
    b0 = b0*x + 1.0;
    b1 = b1*x + 4.0;
    b2 = b2*x + 1.0;
    b3 = b3*x; // Just 0!
    // Factor 6.0:
    b0 /= 6;
    b1 /= 6;
    b2 /= 6;
    b3 /= 6;
    return;
}

/** Solve a block-tridiagonal system by means of the Thomas algorithm */

// Allocate memory for auxiliar buffers:
void allocate_thomas3D( BufferType*& gamma, BufferType*& theta, const SizeBuffer sizes, const unsigned int order )
{
    gamma = new BufferType[order+1];
    theta = new BufferType[order+1];
    for( unsigned int d=0; d<=order; ++d ){
        const SizeType Nd = sizes[d]+1;
        SizeType Md = 1;
        for( unsigned int k=0; k<d; ++k ){ Md *= sizes[k]; }
        gamma[d] = new ElementType[Nd];
        theta[d] = new ElementType[Md*Nd];
    }
    return;
}

// Free memory for auxiliar buffers:
void free_thomas3D( BufferType*& gamma, BufferType*& theta, const unsigned int order )
{
    for( unsigned int d=0; d<=order; ++d ){
        delete[] gamma[d];
        delete[] theta[d];
    }
    delete[] gamma;
    delete[] theta;
    return;
}

// Actually solve it:
void thomas3D( const BufferType values, BufferType* gamma, BufferType* theta, BufferType results, const SizeBuffer sizes, const unsigned int order )
{
    // ---------------------------
    const SizeType N = sizes[order]+1;
    SizeType M = 1;
    for( unsigned int k=0; k<order; ++k ){ M *= sizes[k]; }
    // ---------------------------
    // Forwards swept to compute the recursive multiplicative factors:
    gamma[order][0] = 1;
    for( IndexType d=0; d<M; ++d )
        theta[order][d] = 0.0;
    for( IndexType k=0; k<sizes[order]; ++k ){
        gamma[order][k+1] = -1.0/( gamma[order][k] + 4.0 );
        if( order==0 )
            theta[0][(k+1)*M] = (-gamma[0][k+1])*( 6.0*values[k] - theta[0][k*M] );
        else{
            thomas3D( &values[k*M], gamma, theta,
                                      &(theta[order][(k+1)*M]),
                                      (SizeBuffer)sizes, order-1 );
            for( IndexType d=0; d<M; ++d )
                theta[order][(k+1)*M+d] = (-gamma[order][k+1])*( 6.0*theta[order][(k+1)*M+d] - theta[order][k*M+d] );
        }
    }
    // ---------------------------
    // Backwards swept
    for( IndexType d=0; d<M; ++d )
        results[(sizes[order]-1)*M+d] = (  (theta[order][(N-1)*M+d])/(1.0-gamma[order][N-1])  )*gamma[order][N-1] + theta[order][(N-1)*M+d];
    for( IndexType k=(IndexType)sizes[order]-2; k>=0; --k ){
        for( IndexType d=0; d<M; ++d )
            results[k*M+d] = (gamma[order][k+1])*(results[(k+1)*M+d]) + (theta[order][(k+1)*M+d]);
    }
    // ---------------------------
    return;
}

}

#endif
