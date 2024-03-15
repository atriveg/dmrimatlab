/**
* This header file implements helper types and functions
* to ease multi-threaded coding in Linux, MAC, and Windows
*/

#ifndef _threadHelper_h_
#define _threadHelper_h_

#if (defined (__APPLE__) && defined (__MACH__))
    #define _HAS_POSIX_THREADS
#elif defined (_WIN32)
    #define _HAS_MKL_BLAS
#elif defined (__unix__)
    #define _HAS_POSIX_THREADS
    #define _HAS_MKL_BLAS
#else
    #error "Unable to auto-detect SO/architecture. Cannot define threads structure"
#endif

#ifdef _HAS_POSIX_THREADS
    #include <pthread.h>
    #include <unistd.h>
#else
    #include <sysinfoapi.h>
    #include <windows.h>
    #include <process.h>
#endif

#ifdef _HAS_MKL_BLAS

    #include "omp.h"
    #if defined (__unix__)
        /**
        * This declaration is required because Matlab does not provide any
        * header file interfacing to mkl.so, so we need to manually define
        * this header. According to Intel's MKL documentation, the preferred
        * way to avoid BLAS/LAPACK using their own threads is using the "local"
        * set_num_threads function, so that it won't interfere with other 
        * processes/threads using MKL.
        *
        * NOTE: this implies it will be necessary to manually link against
        * mkl.so when compiling:
        */
        extern "C" {
            int mkl_serv_set_num_threads_local( int );   
        }
    
    #else
        /**
        * In Windows, however, it seems that MinGW is unable to dinamically
        * link against mkl.dll (no mkl.lib available), so we have to rely 
        * on the OpenMP implementation.
        */
    #endif

#else

    /**
    * In MAC OS, BLAS is based on the OpenBLAS implementation, which
    * provides the routines:
    *      void openblas_set_num_threads(int num_threads);
    *      int openblas_get_num_threads(void);
    * However, Matlab does not show the header file for these routines,
    * and we will have to define them here:
    */
    extern "C" {
        void openblas_set_num_threads( int );
        int openblas_get_num_threads( void );
    }

#endif

#ifdef _HAS_POSIX_THREADS
typedef void*     THFCNRET;
typedef pthread_t THHANDLER;
#else
typedef unsigned __stdcall  THFCNRET;
typedef HANDLE              THHANDLER;
#endif

#include <mutex>

#include "mex.h"
#include "../mathsmex/mexToMathsTypes.h"

class DMRIThreader
{
protected:
    // Make the default constructor protected so that
    // it is ensured that a problem size and chunk size
    // are fixed:
    DMRIThreader( );
public:
    DMRIThreader( const SizeType, const SizeType );
    void setThreadInfo( const unsigned int );
    void setPixelCounter( std::mutex*, IndexType* );
    void setProcessSize( const SizeType, const SizeType );
    void claimNewBlock( IndexType*, IndexType* );
    bool threadedProcess( const unsigned int, THFCNRET (*)(void*) );
    inline SizeType getN(){return this->N;}
    unsigned int getThid( void );
private:
    SizeType     N;     // The total number of voxels to process
    SizeType     chunk; // The block-size of voxels to process each time
    unsigned int nth;   // The total number of threads
    std::mutex*  mtx;   // The mutex used to control couter variable pos
    IndexType*   pos;   // The last position claimed by any thread (needs mutex-control)
    unsigned int thid;  // A thread identifier, only used when getThid() is called
};

unsigned int get_number_of_threads( const unsigned int );

unsigned int blas_num_threads(const unsigned int);

int dmriCreateThread( THHANDLER*, THFCNRET (*)(void *), void* );

int dmriJoinThread( THHANDLER );

bool dmriCloseThread( THHANDLER );


#endif // end #ifndef _threadHelper_h_
