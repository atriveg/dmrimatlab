/** Types definitions to serve as an interface from low-level C/C++
 * code to frontend mex functions
 */
#include "mex.h"
#include "float.h"

// - NOTE: for Mac, we have disabled lpthreads until we find a safe way to
//         prevent BLAS/LAPACK from using their own inner threads

/*
#if defined (__unix__) || (defined (__APPLE__) && defined (__MACH__))
#define _HAS_POSIX_THREADS
#endif
*/
#if defined (__unix__)
#define _HAS_POSIX_THREADS
#endif

#if MX_HAS_INTERLEAVED_COMPLEX
typedef mxDouble ElementType;
#else
#define mxGetDoubles mxGetPr
typedef double ElementType;
#endif

typedef ElementType*  BufferType;
typedef mwSize        SizeType;
typedef mwSignedIndex IndexType;
typedef IndexType*    IndexBuffer;
typedef mxArray*      MatrixType;

#ifndef __realmin__
#define __realmin__ DBL_MIN
#endif
#ifndef __eps__
#define __eps__ DBL_EPSILON
#endif
