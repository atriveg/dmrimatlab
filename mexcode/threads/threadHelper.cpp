/**
* This header file implements helper types and functions
* to ease multi-threaded coding in Linux, MAC, and Windows
*/

#include "threadHelper.h"

unsigned int get_number_of_threads( const unsigned int th )
{
    if(th==0)
        return 1;
#ifdef _HAS_POSIX_THREADS
    unsigned int maxthreads = sysconf(_SC_NPROCESSORS_CONF);
#else
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    unsigned int maxthreads = (unsigned int)(sysinfo.dwNumberOfProcessors);
#endif
    return ( th<maxthreads ? th : maxthreads );
}

unsigned int blas_num_threads(const unsigned int nth)
{
    /**
     * It is important to ensure that BLAS will not create its own threads
     * that will interfere with the threads created by this package.
     * Otherwise, the overall number of threads will blow up dramatically
     * putting down the performance. 
     */
#ifdef _HAS_MKL_BLAS

#if defined (__unix__)
    /**
     * According to Intel's MKL documenation, a call to mkl_set_max_threads_local()
     * should be the appropriate choice here. However, Matlab provides only
     * a reduced subset of MKL routines, and the aforementioned is not one 
     * of them. Seeting the environment variable "MKL_NUM_THREADS" will not
     * work either, since the variable is read when the mkl library is loaded,
     * then kept read only. Setting the OMP threads seems to do the trick
     * (though it is documented that "races" might occur upon different calls
     * to this routine.
     */
    if(nth==0)
        return (unsigned int)( omp_get_max_threads() );
    else{
        unsigned int rnth = (unsigned int)( omp_get_max_threads() );
        /**
        * According to Intel's MKL documentation, the preferred way to
        * avoid BLAS/LAPACK using their own threads is using the "local"
        * set_num_threads function, so that it won't interfere with other 
        * processes/threads using MKL.
        *
        * NOTE: this implies it will be necessarty to manually link against
        * mkl.so when compiling:
        */
        mkl_serv_set_num_threads_local( (int)nth );
        /**
        * Otherwise, we should call OpenMP's set_num_threads, but it might
        * collide with other threads using OpenMP:
        */
        // omp_set_num_threads( nth );
        return rnth;
    }
#else
    /**
    * In Windows, however, it seems that MinGW is unable to dinamically
    * link against mkl.dll (no mkl.lib available), so we have to rely 
    * on the OpenMP implementation.
    */
    if(nth==0)
        return (unsigned int)( omp_get_max_threads() );
    else{
        unsigned int rnth = (unsigned int)( omp_get_max_threads() );
        omp_set_num_threads( (int)nth );
        return rnth;
    }
#endif

#else
    /**
     * The default implementation for MAC OS is OpenBLAS, which provides 
     * (according to the documentation) specific routines to set the number
     * of computational threads:
     */
    if(nth==0)
        return (unsigned int)( openblas_get_num_threads() );
    else{
        unsigned int rnth = (unsigned int)( openblas_get_num_threads() );
        openblas_set_num_threads( nth );
        return rnth;
    }
#endif
}

int dmriCreateThread( THHANDLER* tid, THFCNRET (*funct_addr)(void *), void* args )
{
#ifdef _HAS_POSIX_THREADS
    return pthread_create( tid, NULL, funct_addr, args );
#else
    *tid = (THHANDLER)_beginthreadex( NULL, 0, funct_addr, args, 0, NULL );
    return ( (uintptr_t)(*tid)==-1L ? 1 : 0);
#endif
}

int dmriJoinThread( THHANDLER tid )
{
#ifdef _HAS_POSIX_THREADS
    return pthread_join( tid, NULL );
#else
    return (int)( WaitForSingleObject( tid, INFINITE ) );
#endif
}

bool dmriCloseThread( THHANDLER tid )
{
#ifdef _HAS_POSIX_THREADS
    return true;
#else
    return ( CloseHandle(tid) );
#endif
}

DMRIThreader::DMRIThreader( )
{
    this->N = 0;
    this->chunk= 1;
    this->setThreadInfo( 1 );
    this->setPixelCounter( NULL, NULL );
    this->thid = 0;
}

DMRIThreader::DMRIThreader(
    const SizeType N,
    const SizeType chunk
    )
{
    DMRIThreader();
    this->N = N;
    this->chunk= chunk;
}

void DMRIThreader::setThreadInfo( const unsigned int nth )
{
    this->nth = nth;
}

void DMRIThreader::setPixelCounter( std::mutex* mtx, IndexType* pos )
{
    this->mtx = mtx;
    this->pos = pos;
}

void DMRIThreader::setProcessSize( const SizeType N, const SizeType chunk )
{
    this->N     = N;
    this->chunk = chunk;
}

void DMRIThreader::claimNewBlock( IndexType* init, IndexType* end )
{
    // Make sure only one thread is accessing the counter pointed by
    // pos at each time. Reading and writing of this counter must be
    // done always during a mutex lock.
    this->mtx->lock();     // Lock the mutex
    *init = *(this->pos);
    if( *init > (IndexType)(this->N) ){
        *init        = (IndexType)(this->N);
        *end         = (IndexType)(this->N);
        *(this->pos) = (IndexType)(this->N);
    }
    else{
        *(this->pos) += chunk;
        if( *(this->pos) > (IndexType)(this->N) )
            *(this->pos) = (IndexType)(this->N);
        *end = *(this->pos);
    }
    this->mtx->unlock();
    return;
}

bool DMRIThreader::threadedProcess(
    const unsigned int nthreads,
    THFCNRET (* function)(void*)
    )
{
    IndexType  pos = 0;
    std::mutex mtx;

    this->setThreadInfo( nthreads );
    this->setPixelCounter( &mtx, &pos );

    THHANDLER* threads = new THHANDLER[nthreads];
    int*       rets    = new int[nthreads];
    
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads(1);

    for( unsigned int tid=0; tid<nthreads; ++tid )
        rets[tid] = dmriCreateThread( &(threads[tid]), function, (void*)(this) );
    for( unsigned int tid=0; tid<nthreads; ++tid )
        dmriJoinThread( threads[tid] );
    for( unsigned int tid=0; tid<nthreads; ++tid )
        dmriCloseThread( threads[tid] );

    // Revert BLAS threads usage to its default:
    blas_num_threads(blas_threads);

    delete[] threads;
    delete[] rets;

    return false;
}

unsigned int DMRIThreader::getThid( void )
{
    // Get a unique thread identifier if needed,
    // use the mutex lock.
    this->mtx->lock();     // Lock the mutex
    unsigned int _thid = this->thid;
    (this->thid)++;
    this->mtx->unlock();
    return _thid;
}
