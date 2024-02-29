#ifndef _matrixCalculus_cxx
#define _matrixCalculus_cxx

#include "matrixCalculus.h"
#include "math.h"

namespace mataux
{
    /** Fix all elements of the array to a scalar value
     * NOTE: This is done IN PLACE!
     */
    void setValueMxArray(
        BufferType in,
        const SizeType M,
        const SizeType N,
        const ElementType value )
    {
        for( unsigned long pos=0; pos<M*N; ++pos )
            in[pos] = value;
    }
    
    /** Perform element-wise operation of each element with a scalar
     * NOTE: This is done IN PLACE!
     */
    void scalaropMxArray(
        BufferType in,
        const SizeType M,
        const SizeType N,
        const ElementType op,
        const ScalarOperatorType type )
    {
        ptrdiff_t N_ = M*N;
        ptrdiff_t inc_ = 1;
        ElementType op_ = op;
        switch(type){
            case ADD:
                for( unsigned long pos=0; pos<M*N; ++pos ){ in[pos]+=op_; }
                break;
            case SUBTRACT:
                for( unsigned long pos=0; pos<M*N; ++pos ){ in[pos]-=op_; }
                break;
            case DIVIDE:
                //for( unsigned long pos=0; pos<M*N; ++pos ){ in[pos]/=op_; }
                op_ = 1.0f/op_;
                dscal( &N_, &op_, in, &inc_ );
                break;
            case MULTIPLY:
                //for( unsigned long pos=0; pos<M*N; ++pos ){ in[pos]*=op_; }
                dscal( &N_, &op_, in, &inc_ );
                break;
        }
    }
    
    /** Compute the product of two mxArrays representing double matrixes and return the result
        as another mxArray (MxN)*(NxP) = (MxP) */
    void multiplyMxArrays( 
        const BufferType in1,
        const BufferType in2,
        BufferType out,
        const SizeType M,
        const SizeType N,
        const SizeType P )
    {
        // Multiplying two matrices is no piece of
        // cake if they are large enough. A regular loop
        // like this:
        /*
        unsigned long pos = 0;
        for(unsigned long p=0; p<P; ++p ){
            for(unsigned long m=0; m<M; ++m, ++pos ){
                out[pos] = 0;
                for( unsigned long n=0; n<N; ++n )
                    out[pos] += in1[M*n+m]*in2[N*p+n];
            }
        }
        */
        // may be one order of magnitude slower than calling
        // BLAS' dgemm
        ptrdiff_t M_ = (ptrdiff_t)M;
        ptrdiff_t N_ = (ptrdiff_t)N;
        ptrdiff_t P_ = (ptrdiff_t)P;
        ElementType alpha = 1.0f;
        ElementType beta = 0.0f;
        dgemm(
            "N",    // First matrix is not transposed
            "N",    // Second matrix is not transposed
            &M_,    // First matrix (not transposed) has M rows
            &P_,    // Second matrix (not transposed) has P columns
            &N_,    // First matrix (not transoposed) has N columns
            &alpha, // scale applied to in1*in2
            in1,    // First matrix
            &M_,    // Leading dimension of the first matrix
            in2,    // Second matrix
            &N_,    // Leading dimension of the second matrix
            &beta,  // Scale applied to the output (not used)
            out,    // The output matrix
            &M_     // Leading dimension of the output matrix
        );
    }

    /** Compute the product of two mxArrays representing double matrixes (with the
     *  second one transposed) and return the result as another
     *  mxArray (MxN)*(PxN)' = (MxP)
     *   Sizes should be:
     *      in1: M x N
     *      in2: P x N
     *      out: M x P
     */
    void multiplyMxArraysTranspose(
        const BufferType in1,
        const BufferType in2,
        BufferType out,
        const SizeType M,
        const SizeType N,
        const SizeType P )
    {
        // Multiplying two matrices is no piece of
        // cake if they are large enough. A regular loop
        // like this:
        /*
        unsigned long pos = 0;
        for(unsigned long p=0; p<P; ++p ){
            for(unsigned long m=0; m<M; ++m, ++pos ){
                out[pos] = 0;
                for( unsigned long n=0; n<N; ++n )
                    out[pos] += in1[M*n+m]*in2[N*n+p];
            }
        }
        */
        // may be one order of magnitude slower than calling
        // BLAS' dgemm
        ptrdiff_t M_ = (ptrdiff_t)M;
        ptrdiff_t N_ = (ptrdiff_t)N;
        ptrdiff_t P_ = (ptrdiff_t)P;
        ElementType alpha = 1.0f;
        ElementType beta  = 0.0f;
        dgemm(
            "N",    // First matrix is not transposed
            "T",    // Second matrix IS transposed
            &M_,    // First matrix (not transposed) has M rows
            &P_,    // Second matrix (not transposed) has P columns
            &N_,    // First matrix (not transoposed) has N columns
            &alpha, // scale applied to in1*in2
            in1,    // First matrix
            &M_,    // Leading dimension of the first matrix
            in2,    // Second matrix
            &P_,    // Leading dimension of the second matrix
            &beta,  // Scale applied to the output (not used)
            out,    // The output matrix
            &M_     // Leading dimension of the output matrix
        );
    }
    
    /** Given a matrix A with size MxN, compute A^T*A with size NxN */
    void transposeMultiplyMxArray(
        const BufferType in,
        BufferType out,
        const SizeType M,
        const SizeType N
    )
    {
        // Multiplying two matrices is no piece of
        // cake if they are large enough. A regular loop
        // like this:
        /*
        for( IndexType nr=0; nr<(IndexType)N; ++nr ){
            for( IndexType nc=nr; nc<(IndexType)N; ++nc ){
                out[N*nr+nc] = 0.0;
                for( IndexType p=0; p<(IndexType)M; ++p )
                    out[N*nr+nc] += (in[M*nc+p]) * (in[M*nr+p]);
                out[N*nc+nr] = out[N*nr+nc];
            }
        }
        */
        // may be one order of magnitude slower than calling
        // BLAS' dsyrk
        ptrdiff_t M_ = (ptrdiff_t)M;
        ptrdiff_t N_ = (ptrdiff_t)N;
        ElementType alpha = 1.0f;
        ElementType beta  = 0.0f;
        dsyrk(
            "U",    // Upper part computed
            "T",    // Compute in^T*in instead of in*in^T
            &N_,    // The size of the output matrix
            &M_,    // Since second argument is "T", the number of rows of in
            &alpha, // Just 1
            in,     // Input matrix
            &M_,    // Leading dimension of the input
            &beta,  // Just 0
            out,    // Output matrix
            &N_     // Leading dimension of the output
        );
        // This is a generic function, and we don't know if the
        // lower part of out will be further used. Hence, we must
        // fill it just in case (dsyrk won't do it for us)
        for( IndexType c=0; c<(IndexType)N; ++c ){
            for( IndexType r=c+1; r<(IndexType)N; ++r )
                out[N*c+r] = out[N*r+c];
        }
    }
    
    /** Compute the sum of two mxArrays representing double matrixes and return the result
        as another mxArray (MxN)+(MxN) = (MxN) */
    void addMxArrays(
        const BufferType in1,
        const BufferType in2,
        BufferType out,
        const SizeType M,
        const SizeType N )
    {
        unsigned long pos = 0;
        for(unsigned long pos=0; pos<M*N; ++pos )
            out[pos] = in1[pos]+in2[pos];
    }
    
    /** Compute the subtraction of two mxArrays representing double matrixes and return the result
        as another mxArray (MxN)-(MxN) = (MxN) */
    void subtractMxArrays(
        const BufferType in1,
        const BufferType in2,
        BufferType out,
        const SizeType M,
        const SizeType N )
    {
        
        unsigned long pos = 0;
        for(unsigned long pos=0; pos<M*N; ++pos )
            out[pos] = in1[pos]-in2[pos];
    }
    
    /** Transpose a mxArray and return another mxArray (MxN)' = (NxM) */
    void transposeMxArray(
        const BufferType in,
        BufferType out,
        const SizeType M,
        const SizeType N )
    {
        unsigned long pos = 0;
        for(unsigned long m=0; m<M; ++m ){
            for(unsigned long n=0; n<N; ++n, ++pos )
                out[pos] = in[M*n+m];
        }
    }
    
    /**
     * The following function uses Lapack to:
     * 1- Compute the norm of the input matrix in, whose size is MxM
     * 2- Compute the LU factorization of in
     * 3- Compute an estimate of the reciprocal condition number of in
     * 4- If the condition number is larger than a given threshold, rcnth, invert in using the LU factorization
     * Three work buffers have to be provided from outside:
     *   pivot1: Mx1
     *   pivot2: Mx1
     *   work:   4*M
     * Returns true if the inversion is possible, false otherwise
     * BEWARE!!! in is modified in-place
     */
    IndexType checkAndInvertMxArray(
        BufferType in,
        const SizeType M,
        ElementType rcnth,
        IndexBuffer pivot1,
        IndexBuffer pivot2,
        const SizeType lwork,
        BufferType work
    )
    {
        IndexType info;
        ElementType anorm;
        ElementType rcond;
        ptrdiff_t MM = (ptrdiff_t)M;
        /** First, compute the norm of the input matrix */
        anorm = dlange( "1", &MM, &MM, in, &MM, work ); // in remains untouched
        /** Second, compute L-U factorization using dgetrf */
        dgetrf( &MM, &MM, in, &MM, pivot1, &info ); // in is modified here!!!
        if(info!=0)
            return (IndexType)(-1);
        /** Third, estimate the reciprocal condition number */
        dgecon( "1", &MM, in, &MM, &anorm, &rcond, work, pivot2, &info ); // in remains untouched
        if(info!=0)
            return (IndexType)(-2);
        if( rcond<rcnth )
            return (IndexType)(-3);
        /** Finally, invert the matrix */
        ptrdiff_t llwork = (ptrdiff_t)lwork;
        dgetri( &MM, in, &MM, pivot1, work, &llwork, &info );
        if(info!=0)
            return (IndexType)(-4);
        return 0;
    }
    
    /**
     * The following function uses Lapack to:
     * 1- Compute the norm of the symmetric input matrix in, whose size is MxM
     * 2- Compute the LU factorization of in
     * 3- Compute an estimate of the reciprocal condition number of in
     * 4- If the condition number is larger than a given threshold, rcnth, invert in using the LU factorization
     * Three work buffers have to be provided from outside:
     *   pivot1: Mx1
     *   pivot2: Mx1
     *   work:   P*M, with P>=3 (size should be fixed in terms of dsytrf
     * Returns true if the inversion is possible, false otherwise
     * BEWARE!!! in is modified in-place
     */
    IndexType checkAndInvertMxSymArray(
        BufferType in,
        const SizeType M,
        ElementType rcnth,
        IndexBuffer pivot1,
        IndexBuffer pivot2,
        const SizeType lwork,
        BufferType work
    )
    {
        IndexType info;
        ElementType anorm;
        ElementType rcond;
        ptrdiff_t MM = (ptrdiff_t)M;
        ptrdiff_t llwork = (ptrdiff_t)lwork;
        /** First, compute the norm of the input matrix */
        anorm = dlansy( "1", "U", &MM, in, &MM, work ); // in remains untouched
        /** Second, compute L-U factorization using dsytrf */
        dsytrf( "U", &MM, in, &MM, pivot1, 
                work, &llwork, &info ); // in is modified here!!!
        if(info!=0)
            return (IndexType)(-1);
        /** Third, estimate the reciprocal condition number */
        dsycon( "U", &MM, in, &MM, pivot1,
                &anorm, &rcond, work, pivot2, &info ); // in remains untouched
        if(info!=0)
            return (IndexType)(-2);
        if( rcond<rcnth )
            return (IndexType)(-3);
        /** Finally, invert the matrix */
        dsytri( "U", &MM, in, &MM, pivot1, work, &info );
        if(info!=0)
            return (IndexType)(-4);
        // This is a generic function, and we don't know if the
        // lower part of in will be further used. Hence, we must
        // fill it just in case (dsytri won't do it for us)
        for( IndexType c=0; c<(IndexType)M; ++c ){
            for( IndexType r=c+1; r<(IndexType)M; ++r )
                in[M*c+r] = in[M*r+c];
        }
        return 0;
    }
    
    /**
     * The following function uses Lapack to:
     * 1- Compute the norm of the symmetric, positive definite input matrix in, whose size is MxM
     * 2- Compute the Cholesky factorization of in
     * 3- Compute an estimate of the reciprocal condition number of in
     * 4- If the condition number is larger than a given threshold, rcnth, invert in using the factorization
     * Three work buffers have to be provided from outside:
     *   pivot1: Mx1
     *   pivot2: Mx1
     *   work:   3*M
     * Returns true if the inversion is possible, false otherwise
     * BEWARE!!! in is modified in-place
     */
    IndexType checkAndInvertMxPosArray(
        BufferType in,
        const SizeType M,
        ElementType rcnth,
        IndexBuffer pivot1,
        IndexBuffer pivot2,
        const SizeType lwork,
        BufferType work
    )
    {
        IndexType info;
        ElementType anorm;
        ElementType rcond;
        ptrdiff_t MM = (ptrdiff_t)M;
        ptrdiff_t llwork = (ptrdiff_t)lwork;
        /** First, compute the norm of the input matrix */
        anorm = dlansy( "1", "U", &MM, in, &MM, work ); // in remains untouched
        /** Second, compute Cholesky's factorization using dpotrf */
        dpotrf( "U", &MM, in, &MM, &info ); // in is modified here!!!
        if(info!=0) // Meaning: the matrix is not P.D.
            return (IndexType)(-1);
        /** Third, estimate the reciprocal condition number */
        dpocon( "U", &MM, in, &MM,
                &anorm, &rcond, work, pivot2, &info ); // in remains untouched
        if(info!=0)
            return (IndexType)(-2);
        if( rcond<rcnth )
            return (IndexType)(-3);
        /** Finally, invert the matrix */
        dpotri( "U", &MM, in, &MM, &info );
        if(info!=0)
            return (IndexType)(-4);
        // This is a generic function, and we don't know if the
        // lower part of in will be further used. Hence, we must
        // fill it just in case (dpotri won't do it for us)
        for( IndexType c=0; c<(IndexType)M; ++c ){
            for( IndexType r=c+1; r<(IndexType)M; ++r )
                in[M*c+r] = in[M*r+c];
        }
        return 0;
    }
    
    /** Matrix division A\B, i.e. A^(-1)B
     * A (in1) has size M x M
     * B (in2) has size M x N
     * C (out) has size M x N
     * (MxM)\(MxN) = (MxN)
     */
    // Taken from example matrixDivide.c
    IndexType divideMxArrays(
        const BufferType in1,
        const BufferType in2,
        BufferType out,
        const SizeType M,
        const SizeType N )
    {
        // -----------
        IndexType* iPivot = new IndexType[M];
        // -----------
        BufferType bwork = new ElementType[M*M];
        memcpy( bwork, in1, M*M*sizeof(ElementType) );
        memcpy( out,   in2, M*N*sizeof(ElementType) );
        // -----------
        IndexType info;
        // -----------
        ptrdiff_t MM = (ptrdiff_t)M;
        ptrdiff_t NN = (ptrdiff_t)N;
        dgesv( &MM, &NN, bwork, &MM, iPivot, out, &MM, &info );
        // -----------
        delete[] bwork;
        delete[] iPivot;
        // -----------
        return info;
    }
    
    /** Sum all elements of an N-D mxArray */
    ElementType sumMxNDArray(
        const BufferType in,
        const SizeType M )
    {
        double sum = 0.0f;
        for( SizeType m=0; m<M; ++m )
            sum += (double)(in[m]);
        return sum;
    }
    
    /** Compute all non-null elements of an ND array */
    SizeType nonnullsMxNDArray(
        const BufferType in,
        const SizeType M )
    {
        SizeType nonnulls = 0;
        for( SizeType m=0; m<M; ++m ){
            if( notNull(in[m]) )
                ++nonnulls;
        }
        return nonnulls;
    }
    
    bool notNull(const ElementType in)
    {
        return( (in>10*mxGetEps()) || (in<-10*mxGetEps()) );
    }
    
    /** Compute the trace of a square matrix */
    ElementType traceMxArray( const BufferType matrix, const SizeType M )
    {
        ElementType trace = 0.0;
        for( IndexType m=0; m<(IndexType)M; ++m )
            trace += matrix[M*m+m];
        return trace;
    }
    
    /** Add a constant to the diagonal of a square matrix
     * THIS IS DONE IN PLACE!!
     */
    void addIdentityMxArray( BufferType in, const ElementType damp, const SizeType M )
    {
        for( IndexType m=0; m<(IndexType)M; ++m )
            in[M*m+m] += damp;
    }
    
} // end namespace mataux

#endif // #ifndef _matrixCalculus_cxx
