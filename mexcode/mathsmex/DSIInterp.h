/**
 *
 * DSIInterp.h
 *
 * This is auxiliar code to perform interpolation of the EAP
 * computed in a regular Cartesian grid (i.e. the output of
 * HYDI-DSI) by using either nearest neoghbor, or linear or
 * cubic B-spline interpolation
 */

#ifndef _DSIInterp_h_
#define _DSIInterp_h_

#include <cmath>
#include "mexToMathsTypes.h"

namespace dsiinterp{

    SizeType allocate_extended_grid( BufferType&, const SizeType&,  const SizeType&,  const SizeType& );

    void fill_extended_grid( const BufferType, BufferType, const SizeType&,  const SizeType&,  const SizeType& );

    bool interpolate_NN( const BufferType,
                                  const SizeType&, const SizeType&, const SizeType&,
                                  const ElementType&, const ElementType&, const ElementType&, ElementType& );

    bool interpolate_LN( const BufferType,
                                  const SizeType&, const SizeType&, const SizeType&,
                                  const ElementType&, const ElementType&, const ElementType&, ElementType& );

    bool interpolate_BS( const BufferType,
                                  const SizeType&, const SizeType&, const SizeType&,
                                  const ElementType&, const ElementType&, const ElementType&, ElementType& );


    void splinesBi( const ElementType&, ElementType&, ElementType&, ElementType&, ElementType& );


    void allocate_thomas3D( BufferType*&, BufferType*&, const SizeBuffer, const unsigned int = 2 );


    void free_thomas3D( BufferType*&, BufferType*&, const unsigned int = 2 );


    void thomas3D( const BufferType, BufferType*, BufferType*, BufferType, const SizeBuffer, const unsigned int = 2 );

}

#endif
