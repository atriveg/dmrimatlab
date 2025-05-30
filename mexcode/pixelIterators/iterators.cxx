/**
* This header file implements helper classes to loop through pixels
* in 3-D images represented as BufferType
*/

#include "iterators.h"

namespace dmriiters{

    /** ********************************************************************************* */
    /** BASE ITERATOR CLASS                                                               */
    /** ********************************************************************************* */

    /**
     * This constructor is protected, so that it cannot be used anywhere
     * outside the methods of this class and subclasses
     */
    Iterator::Iterator()
    {
        this->ndim = 0;
        this->N = 0;
        this->fov = NULL;
        this->img = NULL;
        this->pos = NULL;
        this->start = NULL;
        this->extent = NULL;
    }

    /**
     * This constructor is used to initalize the iterator with a default
     * requested region filling the FOV
     */
    Iterator::Iterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img
    ) : Iterator()
    {
        if(ndim<1)
            throw std::invalid_argument("The N-D image must have at least 1 dimension");
        if(ndim>32)
            throw std::invalid_argument("The N-D image must have at most 32 dimensions");
        this->ndim = ndim;
        this->N = N;
        this->img = img;
        this->fov = new SizeType[ndim];
        memcpy( this->fov, fov, ndim*sizeof(SizeType) );
        this->pos = new IndexType[ndim];
        this->start = new IndexType[ndim];
        for( unsigned int d=0; d<ndim; ++d )
            this->start[d] = 0;
        this->extent = new SizeType[ndim];
        memcpy( this->extent, fov, ndim*sizeof(SizeType) );
        this->Begin();
    }

    /**
     * This constructor is used to provide a custom requested region
     * to iterate across, in the form of an origin and an extent
     */
    Iterator::Iterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const IndexBuffer start,
        const SizeBuffer extent
    ) : Iterator(ndim,fov,N,img)
    {
        memcpy( this->start, start, ndim*sizeof(IndexType) );
        memcpy( this->extent, extent, ndim*sizeof(SizeType) );
        if(!this->CheckRegion())
            throw std::invalid_argument("Requested region outside FOV");
        this->Begin();
    }

    Iterator::~Iterator()
    {
        if(this->fov!=NULL)
            delete[] this->fov;
        if(this->pos!=NULL)
            delete[] this->pos;
        if(this->start!=NULL)
            delete[] this->start;
        if(this->extent!=NULL)
            delete[] this->extent;
    }

    /**
     * Check if the requested region to iterate is in-bounds
     */
    bool Iterator::CheckRegion(void)
    {
        bool good = true;
        for( unsigned int d=0; d<this->ndim-1; ++d ){
            if( (this->start[d]<0) || (this->start[d]>=this->fov[d]) )
                good = false;
            if( this->start[d] + this->extent[d] > this->fov[d] )
                good = false;
        }
        return good;
    }

    /**
     * Move the iterator to its very first position
     */
    void Iterator::Begin(void)
    {
        for( unsigned int d=0; d<this->ndim; ++d )
            this->pos[d] = this->start[d];
        return;
    }

    /**
     * Move the iterator to the next position
     */
    void Iterator::Next(void)
    {
        for( unsigned int d=0; d<this->ndim; ++d ){
            (this->pos[d])++;
            if( (this->pos[d]) == (this->start[d]+this->extent[d]) ){
                if(d!=this->ndim-1)
                    this->pos[d] = this->start[d];
            }
            else
                break;
        }
        return;
    }

    /**
     * Check if the iterator has reached its last position
     */
    bool Iterator::End(void)
    {
        return ( this->pos[this->ndim-1]>=(this->start[this->ndim-1]+this->extent[this->ndim-1]) );
    }

    /**
     * Get the absolute position within the image buffer for
     * the current location of the iterator.
     */
    IndexType Iterator::AbsolutePosition(void)
    {
        IndexType apos = this->pos[this->ndim-1];
        for( int d=(int)(this->ndim)-2; d>=0; --d )
            apos = apos*(this->fov[d]) + (this->pos[d]);
        return( apos );
    }

    /**
     * Get the absolute position within the image buffer for
     * an arbitrary location of the iterator. BEWARE: no bounds
     * checking is performed
     */
    IndexType Iterator::AbsolutePosition(const IndexBuffer pos)
    {
        IndexType apos = pos[this->ndim-1];
        for( int d=(int)(this->ndim)-2; d>=0; --d )
            apos = apos*(this->fov[d]) + pos[d];
        return( apos );
    }

    /**
     * Retrieve the current position of the iterator in
     * an IndexBuffer
     */
    void Iterator::GetIndex(IndexBuffer index)
    {
        memcpy( index, this->pos, (this->ndim)*sizeof(IndexType) );
        return;
    }

    /**
     * Force the iterator to move to a desired position,
     * return true if succeeded, false if failed (because
     * the requested position is out of bounds for the
     * requested region)
     */
    bool Iterator::SetIndex(const IndexBuffer index)
    {
        for( unsigned int d=0; d<this->ndim; ++d ){
            if( index[d] < this->start[d] )
                return false;
            if( index[d] >= this->start[d] + (IndexType)(this->extent[d]) )
                return false;
        }
        memcpy( this->pos, index, (this->ndim)*sizeof(IndexType) );
        return true;
    }

    /**
     * Get the value of the pixel  pointed by the current position of
     * the iterator (note the pixel might be a vecotr, so that a
     * BufferType is used to store it)
     */
    void Iterator::GetPixel(BufferType pixel)
    {
        if(this->img != NULL)
            memcpy( pixel, &(this->img[(this->AbsolutePosition())*(this->N)]), (this->N)*sizeof(ElementType) );
        return;
    }

    /**
     * Set the value of the pixel currently pointed. Return true
     * if succeeded, false otherwise (because the iterator is
     * pointing out of bounds for the requested region)
     */
    bool Iterator::SetPixel(const BufferType pixel)
    {
        for( unsigned int d=0; d<this->ndim; ++d ){
            if( this->pos[d] < this->start[d] )
                return false;
            if( this->pos[d] >= this->start[d] + (IndexType)(this->extent[d]) )
                return false;
        }
        memcpy( &(this->img[(this->AbsolutePosition())*(this->N)]), pixel, (this->N)*sizeof(ElementType) );
        return true;
    }

    /** ********************************************************************************* */
    /** NEIGHBORHOOD ITERATOR CLASS                                                       */
    /** ********************************************************************************* */

    /**
     * This constructor is protected, so that it cannot be used anywhere
     * outside the methods of this class and subclasses
     */
    NeighborhoodIterator::NeighborhoodIterator() : Iterator()
    {
        this->radii = NULL;
        this->fillval = NULL;
        this->nneighs = 0;
        this->offset = NULL;
        this->neigh = NULL;
        this->npos = 0;
    }

    /**
     * Use this constructor to get a default radii 0,0,0, i.e.
     * the same thing as a non-neighbor iterator, and also a
     * default iteration region matching the entire fov of the
     * image
     */
    NeighborhoodIterator::NeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img
    ) : Iterator(ndim,fov,N,img)
    {
        this->radii = new SizeType[ndim];
        for( unsigned int d=0; d<ndim; ++d ){ this->radii[d] = 0; }
        this->AllocateBuffers();
    }

    /**
     * Use this constructor to get a default radii 0,0,0, i.e.
     * the same thing as a non-neighbor iterator, but with a custom
     * iteration region given by an start index and an extent
     */
    NeighborhoodIterator::NeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const IndexBuffer start,
        const SizeBuffer extent
    ) : Iterator(ndim,fov,N,img,start,extent)
    {
        this->radii = new SizeType[ndim];
        for( unsigned int d=0; d<ndim; ++d ){ this->radii[d] = 0; }
        this->AllocateBuffers();
    }

    /**
     * Use this constructor to get a custom neighborhood shape,
     * described by a vector of radii and centered at the pixel
     * being iterated. The iteration region matches the entire
     * fov of the image
     */
    NeighborhoodIterator::NeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const SizeBuffer radii
    ) : NeighborhoodIterator(ndim,fov,N,img)
    {
        memcpy( this->radii, radii, NDIM*sizeof(SizeType) );
        this->AllocateBuffers();
    }

    /**
     * Use this constructor to get a custom neighborhood shape,
     * described by a vector of radii and centered at the pixel
     * being iterated. Use also a custom iteration region given
     * by an start index and an extent
     */
    NeighborhoodIterator::NeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const SizeBuffer radii,
        const IndexBuffer start,
        const SizeBuffer extent
    ) : NeighborhoodIterator(ndim,fov,N,img,start,extent)
    {
        memcpy( this->radii, radii, NDIM*sizeof(SizeType) );
        this->AllocateBuffers();
    }

    /**
     * Destructor
     */
    NeighborhoodIterator::~NeighborhoodIterator()
    {
        if(this->radii != NULL)
            delete[] this->radii;
        if(this->fillval != NULL)
            delete[] this->fillval;
        if(this->neigh != NULL)
            delete[] this->neigh;
        if(this->offset != NULL){
            for(IndexType i=0; i<this->nneighs; ++i )
                delete[] this->offset[i];
            delete[] this->offset;
        }
    }

    void NeighborhoodIterator::AllocateBuffers(void)
    {
        // Compute the total number of neighbors
        // depending on the radii:
        this->nneighs = 1;
        for( unsigned int d=0; d<this->ndim; ++d )
            this->nneighs *= (2*radii[d]+1);
        // Only if the neighborhood is non-trivial
        // we actually need to allocate memory:
        if(this->nneighs!=1){
            this->neigh = new IndexType[this->nneighs];
            this->offset = new IndexBuffer[this->nneighs];
            for(IndexType i=0; i<this->nneighs; ++i )
                this->offset[i] = new IndexType[this->ndim];
            // Now, pre-compute the offsets (i.e. fill up
            // this->offset ):
            this->ComputeOffsets();
        }
        else{
            this->neigh = NULL;
            this->offset = NULL;
        }
        // The fill value is always intialized to NULL, since
        // we will be using Neumann boundary conditions by
        // default:
        this->fillval = NULL;
        this->ComputeOffsets();
        return;
    }

    void NeighborhoodIterator::ComputeOffsets(void)
    {
        if(this->nneighs==1)
            return;
        for( IndexType i=0; i<this->nneighs; ++i ){
            SizeType value   = i;
            SizeType divider = this->nneighs;
            for( int d=(int)(this->ndim-1); d>=0; --d ){
                divider /= (2*(this->radii[d])+1);
                this->offset[i][d] = (IndexType)(value/divider)-(IndexType)(this->radii[d]);
                value %= divider;
            }
        }
        return;
    }

    /**
     * Override the default Neumann boundary condition for out-of-bounds
     * positions of the iterators. If the argument passed is non-NULL, then
     * out-of-bounds will be fill with this constant value, i.e. a Dirichlet
     * boundary condition will be used
     */
    void NeighborhoodIterator::SetFillValue(const BufferType fillval)
    {
        if(fillval!=NULL){
            if(this->fillval==NULL)
                this->fillval = new ElementType[this->N];
            memcpy( this->fillval, fillval, (this->N)*sizeof(ElementType) );
        }
        else{
            if(this->fillval!=NULL){
                delete[] this->fillval;
                this->fillval = NULL;
            }
        }
        return;
    }

    /**
     * Use this method to retrieve the total number of elements in the
     * neighborhood. This can be used to allocate the buffer to store this
     * neighborhhod upon a call to ::GetNeighborhood(). However, note the
     * pixels can be vectors, so that you should allocate:
     *
     * BufferType nhood = new ElementType[ N*iterator.GetNhoodSize() ];
     *
     * with N the number of components of the vector pixel
     */
    SizeType NeighborhoodIterator::GetNhoodSize(void)
    {
        return this->nneighs;
    }

    /**
     * The following two methods are used to iterate through the neighborhood
     * and retrieve neighbors one by one, by doing something like:
     *
     *   // N is the number of components in the vector pixel:
     *   BufferType currentPixel = new ElementType[N];
     *   iterator.Rewind();
     *   while( iterator.Play(currentPixel) ){
     *      // do something with currentPixel
     *   }
     *
     * The first one simply rewinds the position of the neighborhood to its
     * initial position:
     */
    void NeighborhoodIterator::Rewind(void)
    {
        this->npos = 0;
    }

    /**
     * The second one does the actual job: gets the current neighboor (checking
     * for out-of-bounds) and moves towards the next one. It return true if
     * there are still neighboors to be retrieved, false otherwise.
     * The pixels of the neighborhoods are retrieved in the natural order,
     * i.e., if radii={1,1,1} the N components at offset (-1,-1,-1)
     * are retroeved first in currentPixel; then the N components at offset
     * (0,-1,-1), then (1,-1,-1), then (-1,0,-1) ... (1,1,1)
     */
    bool NeighborhoodIterator::Play(BufferType currentPixel)
    {
        if(this->img==NULL)
            return false; // Nothing to do
        if( this->npos >= (IndexType)(this->nneighs) )
            return false; // Passed the end of the neighborhhod
        if(this->nneighs==1)
            this->GetPixel(currentPixel);
        else{
            bool outbounds = false;
            for( unsigned int d=0; d<this->ndim; ++d ){
                this->neigh[d] = this->pos[d] + this->offset[this->npos][d];
                if( this->neigh[d]<0 ){
                    outbounds = true;
                    this->neigh[d] = 0;
                }
                if( this->neigh[d]>=this->fov[d] ){
                    outbounds = true;
                    this->neigh[d] = this->fov[d]-1;
                }
            }
            if( outbounds && (this->fillval!=NULL) ) // Dirichlet boundary condition
                memcpy( currentPixel, this->fillval, (this->N)*sizeof(ElementType) );
            else // Neumann boundary condition
                memcpy( currentPixel, &(this->img[(this->AbsolutePosition(this->neigh))*(this->N)]), (this->N)*sizeof(ElementType) );
        }
        ++this->npos;
        return ( this->npos < (IndexType)(this->nneighs) );
    }

    /**
     * Get the whole neighborhood of the currently pointed
     * pixel, stored in a buffer. This buffer sholud have size
     * prod(2*radii+1)*N, and it is externally maintained. Use:
     *
     *    BufferType nhood = new ElementType[ N*iterator.GetNhoodSize() ];
     *
     * The pixels of the neighborhoods are stored in the natural order,
     * i.e., if radii={1,1,1} the N components at offset (-1,-1,-1)
     * are stored first in nhood; then the N components at offset
     * (0,-1,-1), then (1,-1,-1), then (-1,0,-1) ... (1,1,1)
     */
    void NeighborhoodIterator::GetNeighborhood(BufferType nhood)
    {
        if(nhood==NULL)
            return;
        // Iterate through the neighborhood:
        this->Rewind();
        while(   this->Play( &(nhood[(this->npos)*(this->N)]) )   );
        return;
    }

    /**
     * Get the offset associated to a particular neighbor
     * index, which is stored in the ndim-sized buffer
     * offser.
     */
    void NeighborhoodIterator::GetOffset( const IndexType pos, IndexBuffer offset )
    {
        for( unsigned int d=0; d<this->ndim; ++d )
            offset[d] = 0;
        if( (pos<0) || (pos>=this->nneighs) )
            return;
        memcpy( offset, this->offset[pos], (this->ndim)*sizeof(IndexType) );
        return;
    }

    /** ********************************************************************************* */
    /** DIRECTIONAL ITERATOR CLASS                                                        */
    /** ********************************************************************************* */

    /**
     * This constructor is protected, so that it cannot be used anywhere
     * outside the methods of this class and subclasses
     */
    DirectionalIterator::DirectionalIterator() : Iterator()
    {
        this->dir  = 0;
        this->step = NULL;
        this->init = NULL;
        this->end  = NULL;
    }

    /**
     * This constructor is used to initalize the iterator with a default
     * requested region filling the FOV
     */
    DirectionalIterator::DirectionalIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img
    ) : Iterator(ndim,fov,N,img)
    {
        this->step = new IndexType[this->ndim];
        this->init = new IndexType[this->ndim];
        this->end  = new IndexType[this->ndim];
        this->endDir = 1;
        this->endDir = (~(this->endDir)) << ndim-1;
        this->BeginDir();
        this->Begin();
    }
    /**
     * This constructor is used to provide a custom requested region
     * to iterate across, in the form of an origin and an extent
     */
    DirectionalIterator::DirectionalIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const IndexBuffer start,
        const SizeBuffer extent
    ) : DirectionalIterator(ndim,fov,N,img)
    {
        memcpy( this->start, start, ndim*sizeof(IndexType) );
        memcpy( this->extent, extent, ndim*sizeof(SizeType) );
        this->BeginDir();
        this->Begin();
    }

    DirectionalIterator::~DirectionalIterator()
    {
        if(this->step!=NULL)
            delete[] this->step;
        if(this->init!=NULL)
            delete[] this->init;
        if(this->end!=NULL)
            delete[] this->end;
    }

    /**
     * Go to the beginning of the set of all causal and
     * anticausal directions (and move the iterator to
     * the beginning of such direction)
     */
    void DirectionalIterator::BeginDir( )
    {
        this->dir = 0;
        this->Begin();
        return;
    }

    /**
     * Go to the next direction (and move the iterator to
     * the beginning of such direction)
     */
    void DirectionalIterator::NextDir( )
    {
        this->dir++;
        this->Begin();
        return;
    }

    /**
     * Check if we have reached past the last direction
     */
    bool DirectionalIterator::EndDir( )
    {
        return ( ((this->dir)&(this->endDir)) != 0 );
    }

    /**
     * Move the iterator to its very first position
     */
    void DirectionalIterator::Begin( )
    {
        unsigned char base;
        for( unsigned int d=0, base=1; d<this->ndim; ++d, base=(base<<1) ){
            // Set the proper origin, the proper step and the proper
            // end position depending whether the direction at each
            // dimension is causal/non-causal
            if( ((this->dir)&base)==0 ){
                this->init[d] = this->start[d];
                this->end[d]  = this->start[d] + (IndexType)(this->extent[d]);
                this->step[d] = 1;
            }
            else{
                this->init[d] = this->start[d] + (IndexType)(this->extent[d]) - 1;
                this->end[d]  = this->start[d]-1;
                this->step[d] = -1;
            }
            // Actually initialize the position of the iterator:
            this->pos[d] = this->init[d];
        }
        return;
    }

    /**
     * Move the iterator to its next position
     */
    void DirectionalIterator::Next( )
    {
        for( unsigned int d=0; d<(this->ndim)-1; ++d ){
            this->pos[d] += this->step[d];
            if( this->pos[d] == this->end[d] )
                this->pos[d] = this->init[d];
            else
                return;
        }
        this->pos[this->ndim-1] += this->step[this->ndim-1];
        return;
    }

    /**
     * Check if the iterator has reached past
     * its final position
     */
    bool DirectionalIterator::End( )
    {
        return ( this->pos[this->ndim-1] == this->end[this->ndim-1] );
    }

    /** ********************************************************************************* */
    /** DIRECTIONAL NEIGHBORHOOD ITERATOR CLASS                                           */
    /** ********************************************************************************* */

    /**
     * This constructor is protected, so that it cannot be used anywhere
     * outside the methods of this class and subclasses
     */
    DirectionalNeighborhoodIterator::DirectionalNeighborhoodIterator() : Iterator(), NeighborhoodIterator(), DirectionalIterator() {}

    /**
     * Use this constructor to get a default radii 0,0,0, i.e.
     * the same thing as a non-neighbor iterator, and also a
     * default iteration region matching the entire fov of the
     * image
     */
    DirectionalNeighborhoodIterator::DirectionalNeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img
    ) : Iterator(ndim,fov,N,img), NeighborhoodIterator(ndim,fov,N,img), DirectionalIterator(ndim,fov,N,img) {}

    /**
     * Use this constructor to get a default radii 0,0,0, i.e.
     * the same thing as a non-neighbor iterator, but with a custom
     * iteration region given by an start index and an extent
     */
    DirectionalNeighborhoodIterator::DirectionalNeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const IndexBuffer start,
        const SizeBuffer extent
    ) : Iterator(ndim,fov,N,img,start,extent), NeighborhoodIterator(ndim,fov,N,img,start,extent), DirectionalIterator(ndim,fov,N,img,start,extent) {}

    /**
     * Use this constructor to get a custom neighborhood shape,
     * described by a vector of radii and centered at the pixel
     * being iterated. The iteration region matches the entire
     * fov of the image
     */
    DirectionalNeighborhoodIterator::DirectionalNeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const SizeBuffer radii
    ) : Iterator(ndim,fov,N,img), NeighborhoodIterator(ndim,fov,N,img,radii), DirectionalIterator(ndim,fov,N,img) {}

    /**
     * Use this constructor to get a custom neighborhood shape,
     * described by a vector of radii and centered at the pixel
     * being iterated. Use also a custom iteration region given
     * by an start index and an extent
     */
    DirectionalNeighborhoodIterator::DirectionalNeighborhoodIterator(
        const unsigned int ndim,
        const SizeBuffer fov,
        const SizeType N,
        const BufferType img,
        const SizeBuffer radii,
        const IndexBuffer start,
        const SizeBuffer extent
    ) : Iterator(ndim,fov,N,img,start,extent), NeighborhoodIterator(ndim,fov,N,img,radii,start,extent), DirectionalIterator(ndim,fov,N,img,start,extent) {}

} // end namespace
