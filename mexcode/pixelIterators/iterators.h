/**
* This header file implements helper classes to loop through pixels
* in 3-D images represented as BufferType
*/

#ifndef _iterators_h_
#define _iterators_h_

#include "mex.h"
#include "../mathsmex/mexToMathsTypes.h"
#include <cstring>
#include <stdexcept>

#define NDIM 3

namespace dmriiters{

    class Iterator
    {
    protected:
        unsigned int ndim;   // The number of dimensions of the image
        SizeType N;          // The number of elements each pixel has (0-th dimension), to deal  with vector images
        SizeBuffer fov;      // The size of the 3-D image
        BufferType img;      // The buffer of the image itself
        IndexBuffer pos;     // The current position of the iterator within each dimension
        IndexBuffer start;   // The beginning of the region to iterate across
        SizeBuffer extent;   // The extent of the region to iterate across
    protected:
        Iterator();  // Non-public on purpose
        bool CheckRegion(void);
    public:
        // Minimum required arguments: number of dimensions, fov of the image, pixel
        // size, and buffer pointer:
        Iterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType );
        // Additionally, region to iterate across as a pair of origin, extent:
        Iterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const IndexBuffer, const SizeBuffer );
        ~Iterator();
        void Begin(void);
        void Next(void);
        bool End(void);
        IndexType AbsolutePosition(void);              // Get the absolute position within the image buffer
        IndexType AbsolutePosition(const IndexBuffer); // Get the absolute position within the image buffer
        void GetIndex(IndexBuffer);                    // Get the current pixel location
        bool SetIndex(const IndexBuffer);              // Set the current pixel location
        void GetPixel(BufferType);                     // Get the pixel at the current iteration
        bool SetPixel(const BufferType);               // Set the pixel at the current iteration
    };

    class NeighborhoodIterator : virtual public Iterator
    {
    protected:
        SizeBuffer   radii;   // Neighborhood radii
        BufferType   fillval; // Fill value for out-of-bounds indices
        SizeType     nneighs; // The total number of neighbors
        IndexBuffer* offset;  // The offsets at each neighborhood position
        IndexBuffer  neigh;   // A neighboring position
        IndexType    npos;    // The global position within the neighborhood
    protected:
        NeighborhoodIterator(); // Non-public on purpose
        void AllocateBuffers(void);
        void ComputeOffsets(void);
    public:
        // Minimum required arguments: number of dimensions, fov of the image, pixel
        // size, and buffer pointer (will use a trivial 1x1x1x...x1 neighborhood):
        NeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType );
        // Additionally, shape of the neighborhood to use, in the form of an N-D
        // radii, i.e. a radii [2x3x1] will lead to a [5x7x3] neighborhood centered
        // in the pixel being iterated
        NeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const SizeBuffer );
        // Or, instead, region to iterate across as a pair of origin, extent
        // (will use a trivial 1x1x1x...x1 neighborhood):
        NeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const IndexBuffer, const SizeBuffer );
        // Or both the shape of the neighborhood and the iteration region:
        NeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const SizeBuffer, const IndexBuffer, const SizeBuffer );
        ~NeighborhoodIterator();
        void SetFillValue(const BufferType);            // Override the default Neumann boundary condition
        SizeType GetNhoodSize(void);                    // Get the number of positions in the neighborhood
        void Rewind(void);                              // Restore the neighborhood pointer to its starting position
        bool Play(BufferType);                          // Get the next neighboor
        void GetNeighborhood(BufferType);               // Get the neighborhood of the current pixel
        void GetOffset( const IndexType, IndexBuffer ); // Get the offset of a particular neighbor
    };

    class DirectionalIterator : virtual public Iterator
    {
    protected:
        unsigned int dir;    // The current direction to iterate along
        unsigned int endDir; // A mask to check if we have reached past the last direction
        IndexBuffer step;    // The steps to move towars each dimension
        IndexBuffer init;    // The initial positions at each dimension for the present direction
        IndexBuffer end;     // The ending positions at each dimension for the present direction
    protected:
        DirectionalIterator(); // Non-public on purpose
    public:
        // Minimum required arguments: number of dimensions, fov of the image, pixel
        // size, and buffer pointer:
        DirectionalIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType );
        // Addtionally, region to iterate across as a pair of origin, extent:
        DirectionalIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const IndexBuffer, const SizeBuffer );
        ~DirectionalIterator();
        void BeginDir(); // Go to the first of the 2^NDIM directions
        void NextDir();  // Go to the next direction
        bool EndDir();   // Check if we have reached the final direction
        void Begin();    // Go to the first pixel at the current direction
        void Next();     // Get the next pixel within the current direction
        bool End();      // Check if we have reached the final pixel of the current direction
    };

    class DirectionalNeighborhoodIterator : public NeighborhoodIterator, public DirectionalIterator
    {
    protected:
        DirectionalNeighborhoodIterator(); // Non-public on purpose
        void ComputeOffsets(void);
    public:
        // Minimum required arguments: number of dimensions, fov of the image, pixel
        // size, and buffer pointer (will use a trivial 1x1x1x...x1 neighborhood):
        DirectionalNeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType );
        // Additionally, shape of the neighborhood to use, in the form of an N-D
        // radii, i.e. a radii [2x3x1] will lead to a [5x7x3] neighborhood centered
        // in the pixel being iterated
        DirectionalNeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const SizeBuffer );
        // Or, instead, region to iterate across as a pair of origin, extent
        // (will use a trivial 1x1x1x...x1 neighborhood):
        DirectionalNeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const IndexBuffer, const SizeBuffer );
        // Or both the shape of the neighborhood and the iteration region:
        DirectionalNeighborhoodIterator( const unsigned int, const SizeBuffer, const SizeType, const BufferType, const SizeBuffer, const IndexBuffer, const SizeBuffer );
    };

}

#endif // end #ifndef _iterators_h_
