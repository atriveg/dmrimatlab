/*==========================================================
 * pixelIteratorTest.cpp
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2025- Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "iterators.h"
#include <stdio.h>
#include <iostream>


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    SizeType fov[3] = {3,2,4};
    SizeType N = 3;
    BufferType img = new ElementType[fov[0]*fov[1]*fov[2]*N];
    BufferType pixel = new ElementType[N];
    unsigned long gpos=0;
    for( unsigned int z=0; z<fov[2]; ++z ){
        for( unsigned int y=0; y<fov[1]; ++y ){
            for( unsigned int x=0; x<fov[0]; ++x ){
                img[(x+y*fov[0]+z*fov[0]*fov[1])*N+0] = (ElementType)gpos - 0.5;
                img[(x+y*fov[0]+z*fov[0]*fov[1])*N+1] = (ElementType)gpos - 0.0;
                img[(x+y*fov[0]+z*fov[0]*fov[1])*N+2] = (ElementType)gpos + 0.5;
                gpos++;
            }
        }
    }

    IndexType origin[3] = {0,0,1};
    SizeType  extent[3] = {3,2,1};

    dmriiters::DirectionalIterator iterator( 3, fov, N, img, origin, extent );
    IndexType position[3];

    for( iterator.BeginDir(); !iterator.EndDir(); iterator.NextDir() ){
        std::cout << "NEW DIRECTION: " << std::endl;
        for( iterator.Begin(); !iterator.End(); iterator.Next() ){
            iterator.GetIndex( (IndexBuffer)position );
            IndexType ap1 = iterator.AbsolutePosition();
            IndexType ap2 = iterator.AbsolutePosition((IndexBuffer)position);
            iterator.GetPixel(pixel);
            std::cout << "   [" << position[0] << ", " << position[1] << ", " << position[2] << "] (" << ap1 << " / " << ap2 << ") : {" << pixel[0] << ", " << pixel[1] << ", " << pixel[2] << "}" << std::endl;
        }
    }


    std::cout << __FILE__ << ": " << __LINE__ << std::endl;

    SizeType radii[3];
    radii[0] = 3;
    radii[1] = 1;
    radii[2] = 1;

    dmriiters::DirectionalNeighborhoodIterator iterator2( 3, fov, N, img, radii, origin, extent );

    std::cout << "NON-DIRECTIONAL, NEIGHBORHOOD ITERATOR" << std::endl;

    for( iterator2.BeginDir(); !iterator2.EndDir(); iterator2.NextDir() ){
        std::cout << "NEW DIRECTION: " << std::endl;
        for( iterator2.Begin(); !iterator2.End(); iterator2.Next() ){
            iterator2.GetIndex( (IndexBuffer)position );
            IndexType ap1 = iterator2.AbsolutePosition();
            IndexType ap2 = iterator2.AbsolutePosition((IndexBuffer)position);
            iterator2.GetPixel(pixel);
            std::cout << "   [" << position[0] << ", " << position[1] << ", " << position[2] << "] (" << ap1 << " / " << ap2 << ") : {" << pixel[0] << ", " << pixel[1] << ", " << pixel[2] << "}" << std::endl;
        }
    }

    std::cout << "NEIGHBORHOOD TESTS" << std::endl;
    ElementType fillval[N] = {-0.2,0.1,0.2};
    iterator2.SetFillValue(fillval);
    BufferType nhood = new ElementType[ N*iterator2.GetNhoodSize() ];

    for( iterator2.BeginDir(); !iterator2.EndDir(); iterator2.NextDir() ){
        for( iterator2.Begin(); !iterator2.End(); iterator2.Next() ){
            iterator2.GetIndex( (IndexBuffer)position );
            if( (position[0]==1) && (position[1]==1) && (position[2]==1) ){
                iterator2.GetNeighborhood(nhood);
            }
        }
    }

    unsigned long pos = 0;
    for( int z=-(int)(radii[2]); z<=(int)(radii[2]); ++z ){
        std::cout << "[" << std::endl;
        for( int y=-(int)(radii[1]); y<=(int)(radii[1]); ++y ){
            for( int x=-(int)(radii[0]); x<=(int)(radii[0]); ++x ){
                std::cout << nhood[pos] << "   ";
                pos += N;
            }
            std::cout << std::endl;
        }
        std::cout << "]" << std::endl;
    }

    delete[] nhood;
    delete[] img;
    delete[] pixel;
    return;
}
