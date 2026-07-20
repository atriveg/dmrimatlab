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

#include <string.h>
#include "mexToMathsTypes.h"

#ifdef OCTAVE_BUILD
    #if defined(_SYSTEM_BLAS_BUILD_)
        #include "./include/blas.h"
        #include "./include/lapack.h"
    #elif defined(_SYSTEM_OPENBLAS_BUILD_)
        #include "./include/blas.h"
        #include "openblas/cblas.h"
        #include "openblas/lapack.h"
        #ifndef _USE_OPENBLAS_THREAD_CONTROL
            #define _USE_OPENBLAS_THREAD_CONTROL
        #endif
    #elif defined(_LOCAL_OPENBLAS_BUILD_)
        #include "openblas_config.h"
        #include "./include/blas.h"
        #include "./include/lapack.h"
        #ifndef _USE_OPENBLAS_THREAD_CONTROL
            #define _USE_OPENBLAS_THREAD_CONTROL
        #endif
    #elif defined(_LOCAL_BLIS_BUILD_)
        #include "./include/blas.h"
        #include "./include/lapack.h"
        #ifndef _USE_BLIS_THREAD_CONTROL
            #define _USE_BLIS_THREAD_CONTROL
        #endif
    #elif defined(_MKL_BLAS_BUILD_)
        #include "mkl_blas.h"
        #include "mkl_lapack.h"
        #ifndef _USE_MKL_THREAD_CONTROL
            #define _USE_MKL_THREAD_CONTROL
        #endif
    #else
        #error "Unknown BLAS implementation"
    #endif
#else
    #include "blas.h"
    #include "lapack.h"
#endif

#ifdef OCTAVE_BUILD
    #if defined(_SYSTEM_BLAS_BUILD_) || defined(_SYSTEM_OPENBLAS_BUILD_) || defined(_LOCAL_OPENBLAS_BUILD_) || defined(_LOCAL_BLIS_BUILD_)
        #define LAPACKCALLFCN(FUNC) FUNC##_
        #define BLASCALLFCN(FUNC) FUNC##_
    #elif defined(_MKL_BLAS_BUILD_)
        typedef MKL_INT BLAS_INT;
        #define LAPACKCALLFCN(FUNC) FUNC
        #define BLASCALLFCN(FUNC) FUNC
    #else
        #error "Unknown BLAS implementation"
    #endif
    #ifndef BLAS_INT
        #if defined(blasint)
            typedef blasint BLAS_INT;
        #else
            typedef int BLAS_INT;
        #endif
    #endif
#else
    typedef ptrdiff_t BLAS_INT;
    #define LAPACKCALLFCN(FUNC) FUNC
    #define BLASCALLFCN(FUNC) FUNC
#endif


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

    ElementType normMxNDArray( const BufferType, const SizeType );

    ElementType rmsMxNDArray( const BufferType, const SizeType );

    SizeType nonnullsMxNDArray( const BufferType, const SizeType );

    bool notNull(const ElementType);

    ElementType traceMxArray( const BufferType, const SizeType );

    void addIdentityMxArray( BufferType, const ElementType, const SizeType );

} // End namespace mataux

#endif // _matrixCalculus_h_
