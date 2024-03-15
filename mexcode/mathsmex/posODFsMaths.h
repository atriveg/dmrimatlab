/** 
 * 
 */

#ifndef _posODFsMaths_h_
#define _posODFsMaths_h_

#include "mexToMathsTypes.h"
#include "matrixCalculus.h"

//#define DEBUG_CODE

namespace posODFs
{
    typedef struct ProblemFeatures{
        SizeType N = 0;                  // Number of dMRI samples per voxel
        unsigned int L = 0;              // Maximum order of SH expansions for the (squared root of the) ODF
        double nu = 0.006;               // Laplacian penalty weighting
        BufferType Y = (BufferType)NULL; // The values of the SH for which the ODF is expanded, size Kp*N, for Kp=(2L+1)(2L+2)/2
    } ProblemFeatures;
    
    typedef struct WignerSymbols{ // To be retrieved with shmaths::computeNumberOfSquaredSHFactors in sphericalHarmonics.cxx
        SizeType      P  = 0;
        BufferType    xi = (BufferType)NULL;    // Size P, to be retrieved with shmaths::computeSquaredSHFactors in sphericalHarmonics.cxx
        IndexBuffer   k1 = (IndexBuffer)NULL;   // Size P, to be retrieved with shmaths::computeSquaredSHFactors in sphericalHarmonics.cxx
        IndexBuffer   k2 = (IndexBuffer)NULL;   // Size P, to be retrieved with shmaths::computeSquaredSHFactors in sphericalHarmonics.cxx
        IndexBuffer   k3 = (IndexBuffer)NULL;   // Size P, to be retrieved with shmaths::computeSquaredSHFactors in sphericalHarmonics.cxx
        unsigned int* l3 = (unsigned int*)NULL; // Size P, to be retrieved with shmaths::computeSquaredSHFactors in sphericalHarmonics.cxx
    } WignerSymbols;
    
    typedef struct VoxelFeatures{
        BufferType  E = (BufferType)NULL;       // The dMRI measurements, size N
        // The convolution factors are intended to be designed "shell-wise", i.e.
        // the N dMRI measurements are grouped in M subsets or "shells". Within 
        // each of these groups, the convolution factors to be applied to the SH
        // coefficients will be the same. Hence, together with the convolution
        // factors themselves, we need a N-sized pointer with values in the range
        // 0...M-1 that tells us which shell each measurement belongs to:
        BufferType  lambda = (BufferType)NULL;  // The convolution factors, size (L+1)*M
        IndexBuffer pshell = (IndexBuffer)NULL; // Shell pointer, size N, values in 0...M
        BufferType  psi = (BufferType)NULL;     // The SH coefficients of (squared root of) the ODF, size K=(L+1)(L+2)/2
        ElementType mu = 0.0;                   // The Lagrange multiplier
    } VoxelFeatures;
    
    typedef struct GradientInputs{
        BufferType delta = (BufferType)NULL; // The individual residuals, size N
        BufferType Delta = (BufferType)NULL; // The individual Laplacian contributions, size Kp=(2L+1)(2L+2)/2
    } GradientInputs;
    
    typedef struct HessianInputs{
        BufferType deltap = (BufferType)NULL; // Derivatives of delta, size N*K
        BufferType Deltap = (BufferType)NULL; // Derivatives of delta, size Kp*K
    } HessianInputs;
    
    typedef struct ODFAlgorithmParameters{
        int T = 200;                  // Maximum number of iterations
        unsigned int maxfails = 5;    // For Levenberg-Marquardt's, maximum numbe of consecutive failed iterations before exit
        ElementType thCost = 0.001;   // For Levenberg-Marquardt's, minimum allowed (successful) relative decrease in the cost function
        ElementType thGrad = 0.001;   // For Newton-Raphson's, minimum allowed relative change in the modulus of the gradient
        ElementType thCond = 0.001;   // For Newton-Raphson's, minimum allowed relative change in the unit-norm constraint
        ElementType rho0 = 1.0;       // The initial damping factor for Levenberg-Marquardt's
        ElementType minrcn = 1.0e-5;  // The minimum reciprocal condition number before a matrix is considered singular
        ElementType psi0eps = 0.0001; // The minimum squared value allowed for psi[0] (only for Levenberg-Marquardt's)
    } ODFAlgorithmParameters;
    
    typedef struct NRWorkBuffers{
        GradientInputs ginputs;
        HessianInputs  hinputs;
        BufferType gradient0 = (BufferType)NULL; // (K+1) x 1, the output of ComputeGradient()
        //BufferType gradient1 = (BufferType)NULL; // (K+1) x 1 , the working estimate of the gradient
        BufferType hessian0 = (BufferType)NULL;  // (K+1) x (K+1) symmetric, the output of ComputeHessian()
        BufferType hessian1 = (BufferType)NULL;  // (K+1) x (K+1) symmetric, temporary copy in case Cholesky's factorization fails
        BufferType x0 = (BufferType)NULL;        // (K+1) x 1, the previous estimate of the ODF (and the Lagrange multiplier)
        BufferType x1 = (BufferType)NULL;        // (K+1) x 1, the current estimate of the ODF (and the Lagrange multiplier)
        /** The following buffers are necessary for lapack's calls */
        SizeType    lwork;    // MUST be >= max( 4*(K+1), NB*(K+1) ), where NB is the optimum blocksize computed with Lapacks's ilaenv
        BufferType  work = (BufferType)NULL;     // lwork
        IndexBuffer pivot0 = (IndexBuffer)NULL;  // (K+1)*1
        IndexBuffer pivot1 = (IndexBuffer)NULL;  // (K+1)*1
    } NRWorkBuffers;
    
    typedef struct LMWorkBuffers{
        BufferType delta0 = (BufferType)NULL;    // (N+Kp) x 1, as provided by ComputeResidual()
        BufferType delta1 = (BufferType)NULL;    // (N+Kp) x 1, as provided by ComputeResidual()
        BufferType jacobian = (BufferType)NULL;  // (N+Kp) x (K-1), as provided by ComputeJacobian()
        BufferType hessian0 = (BufferType)NULL;  // (K-1) x (K-1) symmetric, the pseudo-Hessian computed as J^T*J
        BufferType hessian1 = (BufferType)NULL;  // (K-1) x (K-1) symmetric, the damped pseudo-Hessian
        BufferType x0 = (BufferType)NULL;        // (K-1) x 1, the previous estimate of the ODF (except the DC coefficient)
        BufferType x1 = (BufferType)NULL;        // (K-1) x 1, the current estimate of the ODF (except the DC coefficient)
        BufferType gradient = (BufferType)NULL;  // (K-1) x 1
        /** The following buffers are necessary for lapack's calls */
        SizeType    lwork;    // MUST be >= max( 4*(K-1), NB*(K-1) ), where NB is the optimum blocksize computed with Lapacks's ilaenv
        BufferType  work = (BufferType)NULL;     // lwork
        IndexBuffer pivot0 = (IndexBuffer)NULL;  // (K-1)*1
        IndexBuffer pivot1 = (IndexBuffer)NULL;  // (K-1)*1
    } LMWorkBuffers;
    
    void createWignerSymbols( WignerSymbols*, const unsigned int );
    
    void destroyWignerSymbols( WignerSymbols* );
    
    void createVoxelFeatures( VoxelFeatures*, const SizeType, const SizeType, const unsigned int );
    
    void destroyVoxelFeatures( VoxelFeatures* );
    
    void createGradientInputs( GradientInputs*, const SizeType, const unsigned int );
    
    void destroyGradientInputs( GradientInputs* );
    
    void createHessianInputs( HessianInputs*, const SizeType, const unsigned int );
    
    void destroyHessianInputs( HessianInputs* );
    
    void createNRWorkBuffers( NRWorkBuffers*, const SizeType, const unsigned int );
    
    void destroyNRWorkBuffers( NRWorkBuffers* );
    
    void createLMWorkBuffers( LMWorkBuffers*, const SizeType, const unsigned int );
    
    void destroyLMWorkBuffers( LMWorkBuffers* );
    
    /** For Newton-Raphson's approach */
    
    void ComputeLagrangian( const ProblemFeatures*, const WignerSymbols*, const VoxelFeatures*, GradientInputs*, ElementType* );
    
    void ComputeGradient( const ProblemFeatures*, const WignerSymbols*, const VoxelFeatures*, const GradientInputs*, HessianInputs*, BufferType );
    
    void ComputeHessian( const ProblemFeatures*, const WignerSymbols*, const VoxelFeatures*, const GradientInputs*, const HessianInputs*, BufferType );
    
    int NRFitPositiveODF( const ProblemFeatures*, const WignerSymbols*, VoxelFeatures*, const ODFAlgorithmParameters*, NRWorkBuffers*, ElementType*, ElementType*, ElementType* );
    
    /** For Levenberg-Marquardt's approach */
    
    void ComputeResidual( const ProblemFeatures*, const WignerSymbols*, const VoxelFeatures*, ElementType*, BufferType );
    
    void ComputeJacobian( const ProblemFeatures*, const WignerSymbols*, const VoxelFeatures*, BufferType, BufferType );
    
    int LMFitPositiveODF( const ProblemFeatures*, const WignerSymbols*, VoxelFeatures*, const ODFAlgorithmParameters*, LMWorkBuffers*, ElementType*, ElementType* );
    
} // End namespace posODFs

#endif // _posODFsMaths_h_
