#ifndef _posODFsMaths_cxx
#define _posODFsMaths_cxx

#include "posODFsMaths.h"
#include "math.h"
#include "sphericalHarmonics.h"

#ifdef DEBUG_CODE
#include "stdio.h"
#include "string.h"
#endif

namespace posODFs
{
    void createWignerSymbols( WignerSymbols* wigner, const unsigned int L )
    {
        wigner->P = shmaths::computeNumberOfSquaredSHFactors(L);
        wigner->xi = new ElementType[wigner->P];
        unsigned int* l1 = new unsigned int[wigner->P];
        unsigned int* l2 = new unsigned int[wigner->P];
        wigner->l3 = new unsigned int[wigner->P];
        int* m1 = new int[wigner->P];
        int* m2 = new int[wigner->P];
        int* m3 = new int[wigner->P];
        shmaths::computeSquaredSHFactors( L, wigner->xi, l1, l2, wigner->l3, m1, m2, m3 );
        wigner->k1 = new IndexType[wigner->P];
        wigner->k2 = new IndexType[wigner->P];
        wigner->k3 = new IndexType[wigner->P];
        shmaths::unrollEvenSHIndices( wigner->P, l1,         m1, wigner->k1 );
        shmaths::unrollEvenSHIndices( wigner->P, l2,         m2, wigner->k2 );
        shmaths::unrollEvenSHIndices( wigner->P, wigner->l3, m3, wigner->k3 );
        delete[] l1;
        delete[] l2;
        delete[] m1;
        delete[] m2;
        delete[] m3;
    }
    
    void destroyWignerSymbols( WignerSymbols* wigner )
    {
        delete[] wigner->xi;
        delete[] wigner->k1;
        delete[] wigner->k2;
        delete[] wigner->k3;
        delete[] wigner->l3;
    }
    
    void createVoxelFeatures( VoxelFeatures* voxel, const SizeType N, const SizeType M, const unsigned int L )
    {
        unsigned int K = (L+1)*(L+2)/2;
        voxel->E = new ElementType[N];
        voxel->lambda = new ElementType[(L+1)*M];
        voxel->pshell = new IndexType[N];
        voxel->psi = new ElementType[K];
        voxel->mu = 0.0;
    }
    
    void destroyVoxelFeatures( VoxelFeatures* voxel )
    {
        delete[] voxel->E;
        delete[] voxel->lambda;
        delete[] voxel->pshell;
        delete[] voxel->psi;
    }
    
    void createGradientInputs( GradientInputs* ginputs, const SizeType N, const unsigned int L )
    {
        unsigned int Kp = (2*L+1)*(2*L+2)/2;
        ginputs->delta = new ElementType[N];
        ginputs->Delta = new ElementType[Kp];
    }
    
    void destroyGradientInputs( GradientInputs* ginputs )
    {
        delete[] ginputs->delta;
        delete[] ginputs->Delta;
    }
    
    void createHessianInputs( HessianInputs* hinputs, const SizeType N, const unsigned int L )
    {
        unsigned int K = (L+1)*(L+2)/2;
        unsigned int Kp = (2*L+1)*(2*L+2)/2;
        hinputs->deltap = new ElementType[N*(SizeType)K];
        hinputs->Deltap = new ElementType[(SizeType)Kp*(SizeType)K];
    }
    
    void destroyHessianInputs( HessianInputs* hinputs )
    {
        delete[] hinputs->deltap;
        delete[] hinputs->Deltap;
    }
    
    void createNRWorkBuffers( NRWorkBuffers* nrwork, const SizeType N, const unsigned int L )
    {
        unsigned int K = (L+1)*(L+2)/2;
        unsigned int Kp = (2*L+1)*(2*L+2)/2;        
        createGradientInputs( &(nrwork->ginputs), N, L );
        createHessianInputs( &(nrwork->hinputs), N, L );
        nrwork->gradient0 = new ElementType[K+1];
        //nrwork->gradient1 = new ElementType[K+1];
        nrwork->hessian0 = new ElementType[ (SizeType)(K+1)*(SizeType)(K+1) ];
        nrwork->hessian1 = new ElementType[ (SizeType)(K+1)*(SizeType)(K+1) ];
        nrwork->x0 = new ElementType[K+1];
        nrwork->x1 = new ElementType[K+1];
        // Allocate memory for buffers related to checkAndInvertMxPosArray()
        // and checkAndInvertMxSymArray()
        // With Newton-raphson's the matrix to invert is not necessarily
        // positive definite, so we must contemplate both case.
        //  - With checkAndInvertMxPosArray()
        //    + dlansy needs a work buffer >= (K+1)
        //    + dpotrf does not use a work buffer
        //    + dpocon needs a work buffer >= 3(K+1)
        //    + dpotr does not use a work buffer
        //  - With checkAndInvertMxSymArray()
        //    + dlansy needs a work buffer >= (K+1)
        //    + dsytrf needs a size of buffer to be determined by query search
        //    + dsycon needs a work buffer >= 2(K+1)
        //    + dsytri needs a work buffer >= (K+1)
        // So we will use query search with dsytrf and use the maximum between
        // this value and 3.
        ptrdiff_t K_ = (ptrdiff_t)(K+1);
        ptrdiff_t info = 0;
        ptrdiff_t llwork = -1; // So that query search is performed
        ElementType osize = K;
        nrwork->pivot0 = new IndexType[K+1];
        nrwork->pivot1 = new IndexType[K+1];
        dsytrf( "L", &K_, nrwork->hessian0, &K_, nrwork->pivot0, 
                &osize, &llwork, &info );
        if(info==0){
            if((SizeType)osize>3)
                nrwork->lwork = ((SizeType)osize)*((SizeType)(K+1));
            else
                nrwork->lwork = 3*((SizeType)(K+1));
        }
        else // Might be non-optimal, but won't crash
            nrwork->lwork = 3*((SizeType)(K+1));
        nrwork->work = new ElementType[nrwork->lwork];
    }
    
    void destroyNRWorkBuffers( NRWorkBuffers* nrwork )
    {
        destroyGradientInputs( &(nrwork->ginputs) );
        destroyHessianInputs( &(nrwork->hinputs) );
        delete[] nrwork->gradient0;
        //delete[] nrwork->gradient1;
        delete[] nrwork->hessian0;
        delete[] nrwork->hessian1;
        delete[] nrwork->x0;
        delete[] nrwork->x1;
        delete[] nrwork->work;
        delete[] nrwork->pivot0;
        delete[] nrwork->pivot1;
    }
    
    void createLMWorkBuffers( LMWorkBuffers* lmwork, const SizeType N, const unsigned int L )
    {
        unsigned int K = (L+1)*(L+2)/2;
        unsigned int Kp = (2*L+1)*(2*L+2)/2;
        
        lmwork->delta0 = new ElementType[N+(SizeType)Kp];
        lmwork->delta1 = new ElementType[N+(SizeType)Kp];
        lmwork->jacobian = new ElementType[ (N+(SizeType)Kp)*(SizeType)((int)K-1) ];
        lmwork->hessian0 = new ElementType[ (SizeType)((int)K-1)*(SizeType)((int)K-1) ];
        lmwork->hessian1 = new ElementType[ (SizeType)((int)K-1)*(SizeType)((int)K-1) ];
        lmwork->x0 = new ElementType[(int)K-1];
        lmwork->x1 = new ElementType[(int)K-1];
        lmwork->gradient = new ElementType[(int)K-1];
        // Allocate memory for buffers related to checkAndInvertMxPosArray.
        // With levenberg-Marquardt's, the matrix is either Positive Definite
        // or singular, so this is the only function we will call:
        //   + dlansy needs a work buffer >= (K-1)
        //   + dpotrf does not use a work buffer
        //   + dpocon needs a work buffer >= 3(K-1)
        //   + dpotr does not use a work buffer
        // The solution is obvious: allocate 3(K-1)
        lmwork->lwork = 3*(SizeType)((int)K-1);
        lmwork->work = new ElementType[lmwork->lwork];
        lmwork->pivot0 = new IndexType[((int)K-1)];
        lmwork->pivot1 = new IndexType[((int)K-1)];
    }
    
    void destroyLMWorkBuffers( LMWorkBuffers* lmwork )
    {
        delete[] lmwork->delta0;
        delete[] lmwork->delta1;
        delete[] lmwork->jacobian;
        delete[] lmwork->hessian0;
        delete[] lmwork->hessian1;
        delete[] lmwork->x0;
        delete[] lmwork->x1;
        delete[] lmwork->gradient;
        delete[] lmwork->work;
        delete[] lmwork->pivot0;
        delete[] lmwork->pivot1;
    }
    
    /**
     * The next three functions are intended to solve a generic nonlinear minimization
     * problem with an equality constraint: i.e. the function to minimize is the
     * Lagrangian, defined as the original cost function (a sum of squared differences)
     * plus the Lagrange multiplier times the equality constraint.
     * To do so, we aim at nullifying the gradient of the Lagrangian, which is expected
     * to be tackled with Newton-Raphson's algorithm. This implies computing the
     * derivative of the function to null (the gradient), i.e. computing the Hessian
     * matrix of the Lagrangian.
     * Each of the next three functions implement the computation of each of these
     * quantities.
     */
    void ComputeLagrangian(
        const ProblemFeatures* features,
        const WignerSymbols* wigner,
        const VoxelFeatures* voxel,
        GradientInputs* ginputs,
        ElementType* lagrangian
    )
    {
        SizeType K  = ((features->L)+1)*((features->L)+2)/2;
        SizeType Kp = (2*(features->L)+1)*(2*(features->L)+2)/2;
        *lagrangian = -(voxel->mu);
        memcpy( ginputs->delta, voxel->E, (features->N)*sizeof(ElementType) );
        mataux::setValueMxArray( ginputs->Delta, 1, Kp, (ElementType)0 );
        for( IndexType p=0; p<(IndexType)(wigner->P); ++p ){
            for( IndexType n=0; n<(IndexType)(features->N); ++n ){
                ginputs->delta[n] -= (voxel->lambda[(voxel->pshell[n])*(IndexType)(features->L+1)+(IndexType)(wigner->l3[p])/2])
                    * (wigner->xi[p]) * (voxel->psi[wigner->k1[p]]) * (voxel->psi[wigner->k2[p]])
                    * (features->Y[(wigner->k3[p])*(features->N)+n]);
            }
            ginputs->Delta[wigner->k3[p]] -= (wigner->l3[p]) * (wigner->l3[p]+1) * (wigner->xi[p])
                * (voxel->psi[wigner->k1[p]]) * (voxel->psi[wigner->k2[p]]);
        }
        for( IndexType n=0; n<(IndexType)(features->N); ++n )
            *lagrangian += 0.5 * (ginputs->delta[n]) * (ginputs->delta[n]);
        for( IndexType k3=0; k3<(IndexType)Kp; ++k3 )
            *lagrangian += (features->nu) * (ginputs->Delta[k3]) * (ginputs->Delta[k3]) / 2;
        for( IndexType k1=0; k1<(IndexType)K; ++k1 )
            *lagrangian += (voxel->mu) * (voxel->psi[k1]) * (voxel->psi[k1]);
        return;
    }
    
    void ComputeGradient(
        const ProblemFeatures* features,
        const WignerSymbols* wigner,
        const VoxelFeatures* voxel,
        const GradientInputs* ginputs,
        HessianInputs* hinputs,
        BufferType gradient
    )
    {
        SizeType K  = ((features->L)+1)*((features->L)+2)/2;
        SizeType Kp = (2*(features->L)+1)*(2*(features->L)+2)/2;
        mataux::setValueMxArray( hinputs->deltap, features->N, K, (ElementType)0 );
        mataux::setValueMxArray( hinputs->Deltap, Kp, K, (ElementType)0 );
        mataux::setValueMxArray( gradient, 1, K+1, (ElementType)0 );
        for( IndexType p=0; p<(IndexType)(wigner->P); ++p ){
            for( IndexType n=0; n<(IndexType)(features->N); ++n ){
                ElementType tmp = -2 * (voxel->lambda[(voxel->pshell[n])*(IndexType)(features->L+1)+(IndexType)(wigner->l3[p])/2])
                    * (wigner->xi[p]) * (voxel->psi[wigner->k2[p]]) 
                    * (features->Y[(wigner->k3[p])*(features->N)+n]);
                hinputs->deltap[(wigner->k1[p])*(IndexType)(features->N)+n] += tmp;
                gradient[wigner->k1[p]] += (ginputs->delta[n])*tmp;
            }
            ElementType tmp2 = -2.0 * (wigner->l3[p]) * (wigner->l3[p]+1) * (wigner->xi[p]) * (voxel->psi[wigner->k2[p]]);
            hinputs->Deltap[(wigner->k1[p])*(IndexType)Kp+wigner->k3[p]] += tmp2;
            gradient[wigner->k1[p]] += (features->nu) * (ginputs->Delta[wigner->k3[p]]) * tmp2;
        }
        for( IndexType k1=0; k1<(IndexType)K; ++k1 ){
            gradient[k1] += 2 * (voxel->mu) * (voxel->psi[k1]);
            gradient[K]  += (voxel->psi[k1]) * (voxel->psi[k1]);
        }
        gradient[K] -= 1;
        return;
    }
    
    void ComputeHessian(
        const ProblemFeatures* features,
        const WignerSymbols* wigner,
        const VoxelFeatures* voxel,
        const GradientInputs* ginputs,
        const HessianInputs* hinputs,
        BufferType hessian
    )
    {
        SizeType K  = ((features->L)+1)*((features->L)+2)/2;
        SizeType Kp = (2*(features->L)+1)*(2*(features->L)+2)/2;
        mataux::setValueMxArray( hessian, K+1, K+1, (ElementType)0 );
        for( IndexType k1=0; k1<(IndexType)K; ++k1 ){
            for( IndexType k2=k1; k2<(IndexType)K; ++k2 ){
                for( IndexType n=0; n<(IndexType)(features->N); ++n )
                    hessian[k2*(IndexType)(K+1)+k1] += (hinputs->deltap[k1*(IndexType)(features->N)+n]) * (hinputs->deltap[k2*(IndexType)(features->N)+n]);
                for( IndexType k3=0; k3<(IndexType)Kp;++k3 )
                    hessian[k2*(IndexType)(K+1)+k1] += (features->nu) * (hinputs->Deltap[k1*(IndexType)Kp+k3]) * (hinputs->Deltap[k2*(IndexType)Kp+k3]);
                hessian[k1*(IndexType)(K+1)+k2] = hessian[k2*(IndexType)(K+1)+k1];
            }
        }
        for( IndexType p=0; p<(IndexType)(wigner->P); ++p ){
            hessian[(wigner->k2[p])*(IndexType)(K+1)+(wigner->k1[p])] += -2.0 * (features->nu) * (ginputs->Delta[wigner->k3[p]])
                * (wigner->l3[p]) * (wigner->l3[p]+1) * (wigner->xi[p]);
            for( IndexType n=0; n<(IndexType)(features->N); ++n ){
                hessian[(wigner->k2[p])*(IndexType)(K+1)+(wigner->k1[p])] += -2.0 * (ginputs->delta[n])
                    * (voxel->lambda[(voxel->pshell[n])*(IndexType)(features->L+1)+(IndexType)(wigner->l3[p])/2]) * (wigner->xi[p])
                    * (features->Y[(wigner->k3[p])*(features->N)+n]);
            }
        }
        for( IndexType k1=0; k1<(IndexType)K; ++k1 ){
            hessian[k1*(IndexType)(K+1)+k1] += 2 * (voxel->mu);
            hessian[K*(IndexType)(K+1)+k1] = 2 * (voxel->psi[k1]);
            hessian[k1*(IndexType)(K+1)+K] = 2 * (voxel->psi[k1]);
        }
        return;
    }
    
    /**
     * Fit a positive ODF using Newton-Raphson's approach. The initial guess is passed
     * in voxel->psi, and the final outcome should be retrieved from this same buffer
    */
    int NRFitPositiveODF(
        const ProblemFeatures* features,
        const WignerSymbols* wigner,
        VoxelFeatures* voxel,
        const ODFAlgorithmParameters* parameters,
        NRWorkBuffers* work,
        ElementType* lagrangian0,
        ElementType* q0,
        ElementType* q1
    )
    {
#ifdef DEBUG_CODE
        FILE* fid = fopen("nr_trace.dat","w");
#endif
        // Useful parameters to loop the buffers:
        SizeType    K   = ((features->L)+1)*((features->L)+2)/2;
        SizeType    Kp  = (2*(features->L)+1)*(2*(features->L)+2)/2;
        /** Compute a naive approximation to the Lagrange multiplier: */
        // When mu=0, the gradient of the Lagrangian is just the gradient of the cost to optimize:
        voxel->mu = 0.0;
        ComputeLagrangian( features, wigner, voxel, &(work->ginputs), lagrangian0 );
        ComputeGradient( features, wigner, voxel, &(work->ginputs), &(work->hinputs), work->gradient0 );
        // At the optimum, the gradient of the cost function plus mu times the gradient
        // of the constraint should equal 0:
        // grad{Q} + mu * grad{G} = 0 => grad{G}^T * grad{Q} + mu * grad{G}^T * grad{G} = 0
        // then we will guess mu = - grad{G}^T * grad{Q} / || grad{G} ||^2
        // But the gradient of G is just 2*psi, so that || grad{G} ||^2 should be just 4
        for( IndexType k=0; k<(IndexType)K; ++k )
            voxel->mu -= 0.5 * (voxel->psi[k]) * (work->gradient0[k]);
        // The initial iteration is taken from voxel->psi:
        memcpy( work->x0, voxel->psi, K*sizeof(ElementType) );
        work->x0[K] = voxel->mu;
        // Initialize the iterations:
        bool success     = false;
        int  t           = 0;
        while( (t<parameters->T) & (!success) ){
            ++t;
            // Compute the current values of the Lagrangian and the gradient:
            ComputeLagrangian( features, wigner, voxel, &(work->ginputs), lagrangian0 );
            ComputeGradient( features, wigner, voxel, &(work->ginputs), &(work->hinputs), work->gradient0 );
            // Compute the modules of the gradient to check the stop condition:
            *q0  = 0.0;
            for(IndexType k=0; k<(IndexType)K; ++k )
                *q0 += (work->gradient0[k])*(work->gradient0[k]);
            *q0 /= ((ElementType)K);
            *q1  = (work->gradient0[K])*(work->gradient0[K]);
#ifdef DEBUG_CODE
            fprintf(fid,"%i %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f %1.15f ",
                    t,voxel->mu,::sqrt(*q0),work->gradient0[K],*lagrangian0,0.0,parameters->thGrad,parameters->thCond);
#endif
            // Check the stop condition:
            if( (::sqrt(*q0)<=parameters->thGrad) && (::sqrt(*q1)<=parameters->thCond) ){
                success = true;
#ifdef DEBUG_CODE
                fprintf(fid,"1\n");
#endif
                break;
            }
            // If we didn't reach convergence yet, compute the Hessian:
            ComputeHessian( features, wigner, voxel, &(work->ginputs), &(work->hinputs), work->hessian0 );
            // Try to invert the Hessian. First, we will try to use
            // Cholesky's factorization, but this may fail if the
            // Hessian is not positive definite. In the latter case,
            // we will resort to LU factorization. NOTE: matrix
            // iversion matrices modify the input matrix in-place, so
            // that we need a "backup copy" of the Hessian just in
            // case Cholesky's factorization fails.
            memcpy( work->hessian1, work->hessian0, (K+1)*(K+1)*sizeof(ElementType) );
            IndexType iRes = mataux::checkAndInvertMxPosArray( work->hessian0, K+1, parameters->minrcn, work->pivot0,
                                                               work->pivot1, (ptrdiff_t)(work->lwork), work->work );
            if( iRes==-1 ){ // The Hessian is not P.D.
                // We have a second chance: use the function for inverting generic
                // symmetric matrixes:
                memcpy( work->hessian0, work->hessian1, (K+1)*(K+1)*sizeof(ElementType) ); // Restore the Hessian
                iRes = mataux::checkAndInvertMxSymArray( work->hessian0, K+1, parameters->minrcn, work->pivot0,
                                                         work->pivot1, (ptrdiff_t)(work->lwork), work->work );
            }
            if( iRes!=0 ){ // The matrix coouldn't be inverted by any means
                // If we have entered, work->hessian0 is nearly singular, hence
                // the algorithm has failed:
                success = false;
#ifdef DEBUG_CODE
                fprintf(fid,"2\n");
#endif
                break;
            }
            // If we reach here, the Hessian was successfully inverted and we can multiply to get the new step:
            mataux::multiplyMxArrays( work->hessian0, work->gradient0, work->x1, K+1, K+1, 1 );
            for( IndexType r=0; r<(IndexType)(K+1); ++r )
                work->x1[r] = work->x0[r] - work->x1[r];
            memcpy( voxel->psi, work->x1, K*sizeof(ElementType) );
            voxel->mu = work->x1[K];
            memcpy( work->x0, work->x1, (K+1)*sizeof(ElementType) );
#ifdef DEBUG_CODE
            fprintf(fid,"3\n");
#endif
        }
        // The algorithm succeeded it we reach convergence before performing
        // the maximum number of iterations:
#ifdef DEBUG_CODE
        fclose(fid);
#endif
        if(success)
            return t;
        else
            return -1;
    }
    
    /**
     * The next two functions are intended for a different approach: the size
     * of the problem is reduced in one variable by writing the first component
     * of the vector of unknowns (psi_0) as psi_0=sqrt(1 - sum_k psi_k^2).
     * This way, the problem turns out to be an unconstrained least
     * squares one, for which Levenberg-Marquardt's method should be 
     * efficient.
     * The first function computes the residual for this method.
     * The second one computes the Jacobian matrix used in this method.
     * Since we have N equations coming from the dMRI measurements plus
     * Kp=(2L+1)(2L+2)/2 equations coming from the Laplace-Beltrami
     * penalty, the Jacobian has N+Kp rows.
     * Since we have K-1 = (L+1)(L+2)/2-1 variables, the Jacobian has
     * K-1 columns.
     * Hence, the amount of memory to be externally allocated becomes:
     * 
     *    cost: 1 * sizeof(ElementType)
     *    jacobian: ( features->N + Kp )
     *              * ( K - 1 )
     *              * sizeof(ElementType)
     *    delta: ( features->N + Kp )
     *              * sizeof(ElementType)
     * 
     * (the second one is a either the vector of residuals, when computing
     * the cost, or a temporary buffer used to store the derivatives
     * of each term of the sum of squares w.r.t. psi_0, when computing the
     * Jacobian).
     */
    void ComputeResidual(
        const ProblemFeatures* features,
        const WignerSymbols* wigner,
        const VoxelFeatures* voxel,
        ElementType* cost,
        BufferType delta
    )
    {
        SizeType K  = ((features->L)+1)*((features->L)+2)/2;
        SizeType Kp = (2*(features->L)+1)*(2*(features->L)+2)/2;
        *cost = 0.0;
        mataux::setValueMxArray( delta, (features->N+Kp), 1, (ElementType)0 );
        memcpy( delta, voxel->E, (features->N)*sizeof(ElementType) );
        // -----
        ElementType psi0 = 1.0;
        for( IndexType k=1; k<(IndexType)K; ++k )
            psi0 -= (voxel->psi[k]) * (voxel->psi[k]);
        psi0 = ( psi0>=0.0 ? ::sqrt(psi0) : 0.0 );
        // -----
        for( IndexType p=0; p<(IndexType)(wigner->P); ++p ){
            ElementType psik1 = ( wigner->k1[p]==0 ? psi0 : voxel->psi[wigner->k1[p]] );
            ElementType psik2 = ( wigner->k2[p]==0 ? psi0 : voxel->psi[wigner->k2[p]] );
            for( IndexType n=0; n<(IndexType)(features->N); ++n ){
                delta[n] -= (voxel->lambda[(voxel->pshell[n])*(IndexType)(features->L+1)+(IndexType)(wigner->l3[p])/2])
                    * (wigner->xi[p]) * psik1 * psik2
                    * (features->Y[(wigner->k3[p])*(features->N)+n]); 
            }
            delta[ (IndexType)(features->N) + wigner->k3[p] ] -=
                ::sqrt(features->nu) * (wigner->l3[p]) * (wigner->l3[p]+1)
                * (wigner->xi[p]) * psik1 * psik2;
        }
        for( IndexType r=0; r<(IndexType)(features->N+Kp); ++r )
            *cost += delta[r]*delta[r];
        *cost /= 2;
        return;
    }
    
    void ComputeJacobian(
        const ProblemFeatures* features,
        const WignerSymbols* wigner,
        const VoxelFeatures* voxel,
        BufferType jacobian,
        BufferType delta
    )
    {
        SizeType K  = ((features->L)+1)*((features->L)+2)/2;
        SizeType Kp = (2*(features->L)+1)*(2*(features->L)+2)/2;
        mataux::setValueMxArray( jacobian, (features->N+Kp), K-1, (ElementType)0 );
        mataux::setValueMxArray( delta, (features->N+Kp), 1, (ElementType)0 );
        // -----
        ElementType psi0 = 1.0;
        for( IndexType k=1; k<(IndexType)K; ++k )
            psi0 -= (voxel->psi[k]) * (voxel->psi[k]);
        psi0 = ( psi0>=0.0 ? ::sqrt(psi0) : 0.0 );
        // -----
        for( IndexType p=0; p<(IndexType)(wigner->P); ++p ){
            ElementType psik2 = ( wigner->k2[p]==0 ? psi0 : voxel->psi[wigner->k2[p]] );
            for( IndexType n=0; n<(IndexType)(features->N); ++n ){
                if(wigner->k1[p]==0){
                    delta[n] +=
                        2.0 * (voxel->lambda[(voxel->pshell[n])*(IndexType)(features->L+1)+(IndexType)(wigner->l3[p])/2])
                        * (wigner->xi[p]) * psik2
                        * (features->Y[(wigner->k3[p])*(features->N)+n]);
                }
                else{
                    jacobian[(wigner->k1[p]-1)*(IndexType)(features->N+Kp) + n] +=
                        2.0 * (voxel->lambda[(voxel->pshell[n])*(IndexType)(features->L+1)+(IndexType)(wigner->l3[p])/2])
                        * (wigner->xi[p]) * psik2
                        * (features->Y[(wigner->k3[p])*(features->N)+n]);
                }
            }
            if(wigner->k1[p]==0){
                delta[ (IndexType)(features->N) + wigner->k3[p] ] -=
                    2.0 * ::sqrt(features->nu)
                    * (wigner->l3[p]) * (wigner->l3[p]+1) * (wigner->xi[p])
                    * psik2;
            }
            else{
                jacobian[(wigner->k1[p]-1)*(IndexType)(features->N+Kp) + (IndexType)(features->N) + wigner->k3[p] ] -=
                    2.0 * ::sqrt(features->nu)
                    * (wigner->l3[p]) * (wigner->l3[p]+1) * (wigner->xi[p])
                    * psik2;
            }
        }
        if( (psi0 > 10*mxGetEps()) || (psi0 < -10*mxGetEps()) ){
            for( IndexType i=0; i<(IndexType)(features->N+Kp); ++i ){
                for( IndexType k=0; k<(IndexType)(K-1); ++k )
                    jacobian[k*(IndexType)(features->N+Kp)+i] -= delta[i] * (voxel->psi[k+1]) / psi0;
            }
        }
        return;
    }
    
    /**
     * Fit a positive ODF using Levenberg-Marquardt's approach. The initial guess is passed
     * in voxel->psi, and the final outcome should be retrieved from this same buffer
    */
    int LMFitPositiveODF(
        const ProblemFeatures* features,
        const WignerSymbols* wigner,
        VoxelFeatures* voxel,
        const ODFAlgorithmParameters* parameters,
        LMWorkBuffers* work,
        ElementType* Q0,
        ElementType* q0
    )
    {
#ifdef DEBUG_CODE
        FILE* fid = fopen("lm_trace.dat","w");
#endif
        // Useful parameters to loop the buffers:
        SizeType    K   = ((features->L)+1)*((features->L)+2)/2;
        SizeType    Kp  = (2*(features->L)+1)*(2*(features->L)+2)/2;
        ElementType rho = parameters->rho0;
        // The initial iteration is taken from voxel->psi:
        memcpy( work->x0, &(voxel->psi[1]), (K-1)*sizeof(ElementType) );
        // Initialize the outputs passed by reference:
        ComputeResidual( features, wigner, voxel, Q0, work->delta0 );
        *q0 = 0.0;
        ElementType Q1;
        // Initialize the iterations: in the first one, it is always
        // necessary to compute the derivatives:
        bool success = false;
        bool computeDers = true;
        int  t = 0;
        unsigned int fails = 0; // Consecutive iterations failed
        while( (t<parameters->T) & (!success) ){
            ++t;
            if(computeDers){
                ComputeJacobian( features, wigner, voxel, work->jacobian, work->delta1 ); // It's important not to re-use delta0!!
                mataux::transposeMultiplyMxArray( work->jacobian, work->hessian0, (features->N+Kp), K-1 );
                // Compute delta^*jacobian instead of jacobian^T*delta to avoid explicitly transposing:
                mataux::multiplyMxArrays( work->delta0, work->jacobian, work->gradient, 1, (features->N+Kp), K-1 );
                *q0 = 0.0;
                for( IndexType k=0; k<(IndexType)K-1; ++k )
                    *q0 += (work->gradient[k]) * (work->gradient[k]);
            }
#ifdef DEBUG_CODE
            fprintf(fid,"%i %1.15f %1.15f %1.15f ",t,*Q0,*q0,rho);
#endif
            ElementType damp = mataux::traceMxArray( work->hessian0, K-1 )*rho/(K-1);
            memcpy( work->hessian1, work->hessian0, (K-1)*(K-1)*sizeof(ElementType) );
            mataux::addIdentityMxArray( work->hessian1, damp, K-1 );
            // In Levenberg-Marquardt's algorithm we invert a matrix
            // which is positive semidefintie by construction, so that
            // we can always use Cholesky's factorization. If that fails,
            // is because the matrix is indeed singular, so that we won't
            // do any furthe effort
            if( mataux::checkAndInvertMxPosArray(work->hessian1,K-1,parameters->minrcn,work->pivot0,work->pivot1,(ptrdiff_t)(work->lwork),work->work) !=0 ){
                // If we reach here, the matrix cannot be inverted
                computeDers = false;
                rho *= 2.0;
#ifdef DEBUG_CODE
                ElementType dnorm = 1.0;
                for( IndexType r=0; r<(IndexType)(K-1); ++r )
                    dnorm -= (work->x0[r]) * (work->x0[r]);
                dnorm = ::sqrt( (dnorm>0.0 ? dnorm : 0.0 ) );
                fprintf(fid,"%1.15f %1.15f %1.15f\n",-1.0,dnorm,parameters->thCost);
#endif
                continue;
            }
            // If we reach here, the pseudo-Hessian was successfully inverted and we can multiply to get the new step:
            mataux::multiplyMxArrays( work->hessian1, work->gradient, work->x1, K-1, K-1, 1 );
            // The next step is Levenberg-Marquardt's update. 
            ElementType norm = 0.0;
            for( IndexType r=0; r<(IndexType)(K-1); ++r ){
                work->x1[r] += work->x0[r];
                norm += (work->x1[r]) * (work->x1[r]);
            }
            // This block ensures that psi[0] remains real. In "regular voxels"
            // it shouldn't take any effect, since psi[0] (the DC coefficient)
            // should be much greater than any other psi[k]:
            if( norm>1.0-parameters->psi0eps ){
                for( IndexType r=0; r<(IndexType)(K-1); ++r )
                    work->x1[r] = (work->x1[r]) * ::sqrt(1-parameters->psi0eps) / norm;
            }
            // Place the new iteration in the psi buffer to re-compute the cost:
            // (it is unnecessary to fix voxel->psi[0], since this value will
            // never be used by ComputeResidual() ):
            memcpy( &(voxel->psi[1]), work->x1, (K-1)*sizeof(ElementType) );
            ComputeResidual( features, wigner, voxel, &Q1, work->delta1 );
            // Check if the iteration was successful, i.e. if we effectively
            // decreased the cost function:
            if( Q1>*Q0 ){
                // Nope. Increase the damping factor and repeat.
                rho *= 2.0;
                ++fails;
                computeDers = false;
                if( fails > parameters->maxfails )
                    success = true; //This will stop the iterations
#ifdef DEBUG_CODE
                ElementType dnorm = 1.0;
                for( IndexType r=0; r<(IndexType)(K-1); ++r )
                    dnorm -= (work->x0[r]) * (work->x0[r]);
                dnorm = ::sqrt( (dnorm>0.0 ? dnorm : 0.0 ) );
                fprintf(fid,"%1.15f %1.15f %1.15f\n",-2.0,dnorm,parameters->thCost);
#endif
                continue;
            }
            // We did! Check the stop condition:
            ElementType ccost = (*Q0-Q1)/(*Q0+1000*mxGetEps());
            if( ccost<=parameters->thCost )
                success = true; // This will stop the iterations
#ifdef DEBUG_CODE
            ElementType dnorm = 1.0;
            for( IndexType r=0; r<(IndexType)(K-1); ++r )
                dnorm -= (work->x1[r]) * (work->x1[r]);
            dnorm = ::sqrt( (dnorm>0.0 ? dnorm : 0.0 ) );
            fprintf(fid,"%1.15f %1.15f %1.15f\n",ccost,dnorm,parameters->thCost);
#endif
            // Update the iterations:
            *Q0 = Q1;
            memcpy( work->delta0, work->delta1, (features->N+Kp)*sizeof(ElementType) );
            memcpy( work->x0, work->x1, (K-1)*sizeof(ElementType) );
            // We are pointing in the right direction, so that the damping factor
            // may be decreased to take larger steps:
            rho /= 2.0;
            fails = 0;
            // We have updated the estimate, so it is necessary to re-compute 
            // the jacobian for a new guess 
            computeDers = true;
        }
        // The output of the algorithm must be taken from voxel->psi, which
        // could have been wrongly updated to work->x1 in an unsuccessful
        // iteration. Besides, we didn't pay attention to psi[0], so we have
        // to fix both two issues here:
        ElementType psi0 = 1.0;
        for( IndexType r=0; r<(IndexType)(K-1); ++r )
            psi0 -= (work->x0[r]) * (work->x0[r]);
        voxel->psi[0] = ( psi0>0 ? ::sqrt(psi0) : 0 );
        memcpy( &(voxel->psi[1]), work->x0, (K-1)*sizeof(ElementType) );
#ifdef DEBUG_CODE
        fclose(fid);
#endif
        // The algorithm succeeded it we reach convergence before performing
        // the maximum number of iterations:
        if(success)
            return t;
        else
            return -1;
    }    
} // end namespace posODFs

#endif // #ifndef _posODFsMaths_cxx
