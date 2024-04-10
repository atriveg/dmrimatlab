/** Utilities to work with REAL double matrixes from mxArray data types.
 * Basically, it implements simple matrix operations (add, subtract,
 * multiply, divide) working over mxDouble* buffers meant to be obtained
 * from mxArray's describing 2-D matrixes.
 * Note that all mxDouble* (a.k.a. BufferType) buffers are assumed to be
 * allocated and freed from the calling function, so no sanity checks
 * are performed in any case.
 */

#ifndef _matrixCalculus_h_
#define _matrixCalculus_h_

#include "mexToMathsTypes.h"
#include "lapack.h"
#include "blas.h"
#include <string.h>

namespace mataux
{
    typedef enum{ADD,SUBTRACT,MULTIPLY,DIVIDE} ScalarOperatorType;
    
    void setValueMxArray( BufferType, const SizeType, const SizeType, const ElementType );
    
    void scalaropMxArray( BufferType, const SizeType, const SizeType, const ElementType, const ScalarOperatorType );
    
    void multiplyMxArrays( const BufferType, const BufferType, BufferType, const SizeType, const SizeType, const SizeType );

    void multiplyMxArraysTranspose( const BufferType, const BufferType, BufferType,
                                    const SizeType, const SizeType, const SizeType );
    
    void transposeMultiplyMxArray( const BufferType, BufferType, const SizeType, const SizeType );
    
    void addMxArrays( const BufferType, const BufferType, BufferType, const SizeType, const SizeType );
    
    void subtractMxArrays( const BufferType, const BufferType, BufferType, const SizeType, const SizeType );
    
    void transposeMxArray( const BufferType, BufferType, const SizeType, const SizeType );
    
    IndexType checkAndInvertMxArray( BufferType, const SizeType, const ElementType, IndexBuffer, IndexBuffer, SizeType, BufferType );
    
    IndexType checkAndInvertMxSymArray( BufferType, const SizeType, const ElementType, IndexBuffer, IndexBuffer, SizeType, BufferType );
    
    IndexType checkAndInvertMxPosArray( BufferType, const SizeType, const ElementType, IndexBuffer, IndexBuffer, SizeType, BufferType );
    
    IndexType divideMxArrays( const BufferType, const BufferType, BufferType, const SizeType, const SizeType );
    
    ElementType sumMxNDArray( const BufferType, const SizeType );
    
    ElementType rmsMxNDArray( const BufferType, const SizeType );
    
    SizeType nonnullsMxNDArray( const BufferType, const SizeType );
    
    bool notNull(const ElementType);
    
    ElementType traceMxArray( const BufferType, const SizeType );
    
    void addIdentityMxArray( BufferType, const ElementType, const SizeType );
    
} // End namespace mataux

#endif // _matrixCalculus_h_
