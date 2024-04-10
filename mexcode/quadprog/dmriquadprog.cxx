#ifndef _dmriquadprog_cxx_
#define _dmriquadprog_cxx_

#include "dmriquadprog.h"

namespace dmriqpp
{
    /**
     * Utility function to allocate buffers of doubles
     */
    void allocateDouble( double*& buffer, const unsigned int M, const unsigned int N )
    {
        if(M*N>0)
            buffer = new double[M*N];
        else
            buffer = NULL;
        return;
    }
    
    /**
     * Utility function to free buffers of doubles
     */
    void freeDouble( double*& buffer  )
    {
        if( buffer!=NULL ){
            delete[] buffer;
            buffer = NULL;
        }
        return;
    }
    
    /**
     * Allocate all necessary memory to store the problem.
     * TAKE CARE:
     *   - ci and ce must equal n or be 0, otherwise
     *     the function will return an error
     *   - Aeq and A will be stored transposed (hence the
     *     names Aeqt and At) so that computing C'*y is 
     *     more effcient
     */
    int allocateQPProblem( QPProblem& problem,
                            const unsigned int n, const unsigned int ce,
                            const unsigned int ci, const unsigned int cl,
                            const unsigned int cu )
    {
        problem.Qi = problem.f = NULL;
        problem.Aeqt = problem.beq = NULL;
        problem.At = problem.b = NULL;
        problem.lb = problem.ub = NULL;
        problem.x = problem.leq = problem.l
           = problem.mu = problem.eta = NULL;
        problem.leq0 = problem.l0
           = problem.mu0 = problem.eta0 = NULL;
        
        int status = QPP_OK;
        
        problem.N  = n;
        problem.CE = ce;
        problem.CI = ci;
        problem.CL = cl;
        problem.CU = cu;
        
        if(n==0)
            status = QPP_BAD_ARGS;
        if( (cl>0) && (cl!=n) )
            status = QPP_BAD_ARGS;
        if( (cu>0) && (cu!=n) )
            status = QPP_BAD_ARGS;
        
        allocateDouble(problem.Qi,n,n);
        allocateDouble(problem.f,n,1);
        allocateDouble(problem.Aeqt,n,ce);
        allocateDouble(problem.beq,ce,1);
        allocateDouble(problem.At,n,ci);
        allocateDouble(problem.b,ci,1);
        allocateDouble(problem.lb,cl,1);
        allocateDouble(problem.ub,cu,1);
        allocateDouble(problem.x,n,1);
        allocateDouble(problem.leq,ce,1);
        allocateDouble(problem.l,ci,1);
        allocateDouble(problem.mu,cl,1);
        allocateDouble(problem.eta,cl,1);
        allocateDouble(problem.leq0,ce,1);
        allocateDouble(problem.l0,ci,1);
        allocateDouble(problem.mu0,cl,1);
        allocateDouble(problem.eta0,cl,1);
        
        problem.iters = 0;
        problem.cost = 0.0;
        problem.step = 1.0;
        
        return status;
    }
    
    /**
     * Free the previously allocated buffers
     */
    void freeQPProblem( QPProblem& problem )
    {
        freeDouble(problem.Qi);
        freeDouble(problem.f);
        freeDouble(problem.Aeqt);
        freeDouble(problem.beq);
        freeDouble(problem.At);
        freeDouble(problem.b);
        freeDouble(problem.lb);
        freeDouble(problem.ub);
        freeDouble(problem.x);
        freeDouble(problem.leq);
        freeDouble(problem.l);
        freeDouble(problem.mu);
        freeDouble(problem.eta);
        freeDouble(problem.leq0);
        freeDouble(problem.l0);
        freeDouble(problem.mu0);
        freeDouble(problem.eta0);
        problem.N = problem.CE = problem.CI
            = problem.CL = problem.CU = 0;
        return;
    }
    
    /**
     * Set the inverse of Q, size N x N
     */
    void setQi( QPProblem& problem, const double* Qi )
    {
        if(problem.N>0)
            memcpy( problem.Qi, Qi, (problem.N)*(problem.N)*sizeof(double) );
        return;
    }
    
    /**
     * Set the linear term in the cost function, N x 1
     */
    void setf( QPProblem& problem, const double* f )
    {
        if(problem.N>0)
            memcpy( problem.f, f, (problem.N)*sizeof(double) );
        return;
    }
    
    /**
     * Set the matrix of equality constraints, Aeq. NOTE this
     * matrix is stored transposed, so things are not as easy
     * as directly copying the memory. Input is CE x N, and it
     * has to be reshaped to N x CE 
     */
    void setAeq( QPProblem& problem, const double* Aeq )
    {
        for( unsigned int r=0; r<problem.CE; ++r ){
            for( unsigned int c=0; c<problem.N; ++c )
                problem.Aeqt[c+(problem.N)*r] = Aeq[r+(problem.CE)*c];
        }
        return;
    }
    
    /**
     * Directly set Aeq transposed, N x CE
     */
    void setAeqt( QPProblem& problem, const double* Aeqt )
    {
        if( problem.CE>0 )
            memcpy( problem.Aeqt, Aeqt, (problem.CE)*(problem.N)*sizeof(double) );
        return;
    }
    
    /**
     * Set the vector of equality constraints, CE x 1
     */
    void setbeq( QPProblem& problem, const double* beq )
    {
        if( problem.CE>0 )
            memcpy( problem.beq, beq, (problem.CE)*sizeof(double) );
        return;
    }
    
    /**
     * Set the matrix of inequality constraints, A. NOTE this
     * matrix is stored transposed, so things are not as easy
     * as directly copying the memory. Input is CI x N, and it
     * has to be reshaped to N x CI
     */
    void setA( QPProblem& problem, const double* A )
    {
        for( unsigned int r=0; r<problem.CI; ++r ){
            for( unsigned int c=0; c<problem.N; ++c )
                problem.At[c+(problem.N)*r] = A[r+(problem.CI)*c];
        }
        return;
    }
    
    /**
     * Directly set A transposed, N x CI
     */
    void setAt( QPProblem& problem, const double* At )
    {
        if( problem.CI>0 )
            memcpy( problem.At, At, (problem.CI)*(problem.N)*sizeof(double) );
        return;
    }
    
    /**
     * Set the vector of inequality constraints, CI x 1
     */
    void setb( QPProblem& problem, const double* b )
    {
        if( problem.CI>0 )
            memcpy( problem.b, b, (problem.CI)*sizeof(double) );
        return;
    }
    
    /**
     * Set the vector of lower bounds, lb, CL x 1
     */
    void setlb( QPProblem& problem, const double* lb )
    {
        if(problem.CL>0)
            memcpy( problem.lb, lb, (problem.CL)*sizeof(double) );
        return;
    }
    
    /**
     * Set the vector of upper bounds, ub, CU x 1
     */
    void setub( QPProblem& problem, const double* ub )
    {
        if ( problem.CU>0 )
            memcpy( problem.ub, ub, (problem.CU)*sizeof(double) );
        return;
    }
    
    /**
     * Allocate all necessary memory to store the intermediate.
     * results of the algorithm. NOTE: this function should be
     * called right after an appropriate initialization of the
     * problem (via allocateQPProblem) and an instantiation
     * and customization of a params structure QPParams */
    int allocateQPAuxiliar( const QPProblem& problem, 
                            const QPParams& params, QPAuxiliar& auxiliar )
    {
        const unsigned int n  = problem.N;
        const unsigned int ce = problem.CE;
        const unsigned int ci = problem.CI;
        const unsigned int cl = problem.CL;
        const unsigned int cu = problem.CU;
        const unsigned int c  = (ce+ci+cl+cu);
        
        auxiliar.fdualeq = auxiliar.fdual
           = auxiliar.fduallb = auxiliar.fdualub = NULL;
        auxiliar.dyeq = auxiliar.dy
           = auxiliar.dylb = auxiliar.dyub = NULL;
        auxiliar.CQieqt = auxiliar.CQit
           = auxiliar.CQilbt = auxiliar.CQiubt = NULL;
        auxiliar.resultNx1 = NULL;
        auxiliar.Ctyeq = auxiliar.Cty
           = auxiliar.Ctylb = auxiliar.Ctyub = NULL;
        auxiliar.ueq = auxiliar.u
           = auxiliar.ulb = auxiliar.uub = NULL;
           
        int status = QPP_OK;
        
        if(n==0)
            status = QPP_BAD_ARGS;
        if( (cl>0) && (cl!=n) )
            status = QPP_BAD_ARGS;
        if( (cu>0) && (cu!=n) )
            status = QPP_BAD_ARGS;
        
        allocateDouble(auxiliar.fdualeq,ce,1);
        allocateDouble(auxiliar.fdual,ci,1);
        allocateDouble(auxiliar.fduallb,cl,1);
        allocateDouble(auxiliar.fdualub,cu,1);
        
        allocateDouble(auxiliar.dyeq,ce,1);
        allocateDouble(auxiliar.dy,ci,1);
        allocateDouble(auxiliar.dylb,cl,1);
        allocateDouble(auxiliar.dyub,cu,1);
        
        allocateDouble(auxiliar.CQieqt,n,ce);
        allocateDouble(auxiliar.CQit,n,ci);
        allocateDouble(auxiliar.CQilbt,n,cl);
        allocateDouble(auxiliar.CQiubt,n,cu);
        
        allocateDouble(auxiliar.resultNx1,n,1);
        allocateDouble(auxiliar.Ctyeq,n,1);
        allocateDouble(auxiliar.Cty,n,1);
        allocateDouble(auxiliar.Ctylb,n,1);
        allocateDouble(auxiliar.Ctyub,n,1);
        
        allocateDouble(auxiliar.ueq,ce,1);
        allocateDouble(auxiliar.u,ci,1);
        allocateDouble(auxiliar.ulb,cl,1);
        allocateDouble(auxiliar.uub,cu,1);
        
        return status;
    }
    
    /**
     * Free the previously allocated buffers
     */
    void freeQPAuxiliar( QPAuxiliar& auxiliar )
    {
        freeDouble(auxiliar.fdualeq);
        freeDouble(auxiliar.fdual);
        freeDouble(auxiliar.fduallb);
        freeDouble(auxiliar.fdualub);
        
        freeDouble(auxiliar.dyeq);
        freeDouble(auxiliar.dy);
        freeDouble(auxiliar.dylb);
        freeDouble(auxiliar.dyub);
        
        freeDouble(auxiliar.CQieqt);
        freeDouble(auxiliar.CQit);
        freeDouble(auxiliar.CQilbt);
        freeDouble(auxiliar.CQiubt);
        
        freeDouble(auxiliar.resultNx1);
        freeDouble(auxiliar.Ctyeq);
        freeDouble(auxiliar.Cty);
        freeDouble(auxiliar.Ctylb);
        freeDouble(auxiliar.Ctyub);
        
        freeDouble(auxiliar.ueq);
        freeDouble(auxiliar.u);
        freeDouble(auxiliar.ulb);
        freeDouble(auxiliar.uub);
        
        return;
    }
    
    /**
     * A function to pre-compute the linear term of the
     * quadratic cost function in the dual problem: 
     * (f'*Q^(-1)*C'+u')'. We will store it in four 
     * separated buffers for each set of constraints.
     */
    void setDualf( const QPProblem& problem, QPAuxiliar& auxiliar )
    {
        // Compute first (f')*Q^(-1), which has size 1 x N:
        mataux::multiplyMxArrays( problem.f, problem.Qi, auxiliar.resultNx1,
                                 1, problem.N, problem.N);
        // Now, for each of the four terms comprising the Lagrange
        // multipliers, compute [ (f')*Q^(-1) ] * C' + u', which has size
        // (1 x N) x (N x Cr) = 1 x Cr:
        if( problem.CE>0 ){ // If we actually have equality constraints
            mataux::multiplyMxArrays( auxiliar.resultNx1, problem.Aeqt,
                                      auxiliar.fdualeq, 1, problem.N, problem.CE ); // 1 x CE
            mataux::addMxArrays( auxiliar.fdualeq, problem.beq, auxiliar.fdualeq, problem.CE, 1 ); // 1 x CE
        }
        if( problem.CI>0 ){ // If we actually have inequality constraints
            mataux::multiplyMxArrays( auxiliar.resultNx1, problem.At,
                                      auxiliar.fdual, 1, problem.N, problem.CI ); // 1 x CI
            mataux::addMxArrays( auxiliar.fdual, problem.b, auxiliar.fdual, problem.CI, 1 ); // 1 x CE
        }
        if( problem.CL>0 ){ // If we actually have lower bounds, the matrix is -I_N
            for( unsigned int c=0; c<problem.CL; ++c )
                auxiliar.fduallb[c] = -problem.lb[c]-auxiliar.resultNx1[c]; // 1 x CL = 1 x N
        }
        if( problem.CU>0 ){ // If we actually have upper bounds, the matrix is I_N
            memcpy( auxiliar.fdualub, auxiliar.resultNx1, (problem.CU)*sizeof(double) ); // 1 x CU = 1 x N
            mataux::addMxArrays( auxiliar.fdualub, problem.ub, auxiliar.fdualub, problem.CU, 1 ); // 1 x CU = 1 x N
        }
        return;
    }
    
    /**
     * A function to pre-compute the matrix (C*Q^(-1))', which is
     * necessary to compute the gradient of the dual cost
     * function. It has size N x (CE+CI+CL+CU)
     */
    void setCQi( const QPProblem& problem, QPAuxiliar& auxiliar )
    {
        // Proceed with each of the four types of constraints
        if( problem.CE>0 ){ // If we actually have equality constraints
            // Take advantage of the fact that Qi is symmetric to
            // compute Qi*Aeqt instead of Aeq*Qi
            mataux::multiplyMxArrays( problem.Qi, problem.Aeqt, auxiliar.CQieqt,
                                      problem.N, problem.N, problem.CE ); // N x CE
        }
        if( problem.CI>0 ){ // If we actually have inequality constraints
            // Take advantage of the fact that Qi is symmetric to
            // compute Qi*At instead of A*Qi
            mataux::multiplyMxArrays( problem.Qi, problem.At, auxiliar.CQit,
                                      problem.N, problem.N, problem.CI ); // N x CI
        }
        if( problem.CL>0 ){ // If we actually have lower bounds, C is just -I_N
            for( unsigned long p=0; p<(problem.N)*(problem.CL); ++p )
                auxiliar.CQilbt[p] = -problem.Qi[p]; // N x CL = N x N
        }
        if( problem.CU>0 ){ // If we actually have upper bounds, C is just is I_N
            memcpy( auxiliar.CQiubt, problem.Qi, (problem.N)*(problem.CU)*sizeof(double) ); // N x CU = N x N
        }
    }
    
    /**
     * A function to compute the gradient of the quadratic cost
     * function in the dual problem
     */
    void computeDualGradient( const QPProblem& problem, QPAuxiliar& auxiliar )
    {
        // The gradient reads:
        //   (C*Q^(-1)) * (C'*y) + (f'*Q^(-1)*C'+u')
        // where the term (C*Q^(-1))' is precomputed in auxiliar.CQi_xxx_t
        // and (f'*Q^(-1)*C'+u') is stored in auxiliar.fdual_xxx
        // Hence, it is a matter of checking the constraints we actually
        // have and store each part of the gradient.
        // 1- Compute C'*y, Nx1 (the result is stored in auxiliar.resultNx1
        computeCty( problem, auxiliar );
        // 2- From the four Nx1 terms, compute [ C*Q^(-1) ] * (C'*y) for each of
        //    the four parts. Note we have pre-computed [ C*Q^(-1) ]', hence it is
        //    better to compute (C'*y)' * [ C*Q^(-1) ]' and transpose
        if( problem.CE>0 ){
            mataux::multiplyMxArrays( auxiliar.resultNx1, auxiliar.CQieqt, auxiliar.ueq,
                                      1, problem.N, problem.CE ); // 1 x CE
        }
        if( problem.CI>0 ){
            mataux::multiplyMxArrays( auxiliar.resultNx1, auxiliar.CQit, auxiliar.u,
                                      1, problem.N, problem.CI ); // 1 x CI
        }
        if( problem.CL>0 ){
            mataux::multiplyMxArrays( auxiliar.resultNx1, auxiliar.CQilbt, auxiliar.ulb,
                                      1, problem.N, problem.CL ); // 1 x CL = 1 x N
        }
        if( problem.CU>0 ){
            mataux::multiplyMxArrays( auxiliar.resultNx1, auxiliar.CQiubt, auxiliar.uub,
                                      1, problem.N, problem.CU ); // 1 x CU = 1 x N
        }
        // Add the term [ C*Q^(-1) ] * (C'*y) to the linear term
        // (f'*Q^(-1)*C'+u') to obtain each part of the gradient:
        if(problem.CE>0)
            mataux::addMxArrays( auxiliar.ueq, auxiliar.fdualeq, auxiliar.dyeq, problem.CE, 1 );
        if(problem.CI>0)
            mataux::addMxArrays( auxiliar.u,   auxiliar.fdual,   auxiliar.dy,   problem.CI, 1 );
        if(problem.CL>0)
            mataux::addMxArrays( auxiliar.ulb, auxiliar.fduallb, auxiliar.dylb, problem.CL, 1 );
        if(problem.CU>0)
            mataux::addMxArrays( auxiliar.uub, auxiliar.fdualub, auxiliar.dyub, problem.CU, 1 );
        // We're all done
        return;
    }
    
    /**
     * A function to compute the current cost of the dual problem
     */
    double computeDualCost( const QPProblem& problem, QPAuxiliar& auxiliar )
    {
        double cost = 0.0;
        // 1 - Begin by computing the linear term:
        if(problem.CE>0){
            for( unsigned int c=0; c<problem.CE; ++c )
                cost += auxiliar.fdualeq[c]*problem.leq[c];
        }
        if(problem.CI>0){
            for( unsigned int c=0; c<problem.CI; ++c ){
                if( problem.l[c] > 0.0 )
                    cost += auxiliar.fdual[c]*problem.l[c];
            }            
        }
        if(problem.CL>0){
            for( unsigned int c=0; c<problem.CL; ++c ){
                if( problem.mu[c] > 0.0 )
                    cost += auxiliar.fduallb[c]*problem.mu[c];
            }            
        }
        if(problem.CU>0){
            for( unsigned int c=0; c<problem.CU; ++c ){
                if( problem.eta[c] > 0.0 )
                    cost += auxiliar.fdualub[c]*problem.eta[c];
            }            
        }
        // 2- Compute (C')*y taking advantage of the sparsity:
        computeCty( problem, auxiliar );
        // 3- Compute Q^(-1)*(C'*y), which has size N x 1
        mataux::multiplyMxArrays( problem.Qi, auxiliar.resultNx1, auxiliar.Cty,
                            problem.N, problem.N, 1 ); // N x 1
        // 4- Compute (y'*C)*Q^(-1)*(C'*y), which is 1 x 1
        cost *= 2.0;
        for( unsigned int n=0; n<problem.N; ++n )
            cost += (auxiliar.Cty[n])*(auxiliar.resultNx1[n]);
        cost *= 0.5;
        // We're all done
        return cost;
    }

    /**
     * This function actually proceeds withe the projected gradient descent: the
     * Lagrange multipliers are subtracted "step" times the gradient of the cost
     * function. All the inequality multipliers that become negative are clipped
     * to 0 so that the dual solution stays within the feasible region
     */
    void gradientDescent( QPProblem& problem, const QPAuxiliar& auxiliar, const double& step )
    {
        if( problem.CE>0 ){ // Equality multipliers are not constrained
            for( unsigned int c=0; c<problem.CE; ++c )
                problem.leq[c] = problem.leq0[c] - step*(auxiliar.dyeq[c]);
        }
        if( problem.CI>0 ){
            for( unsigned int c=0; c<problem.CI; ++c ){
                problem.l[c] = problem.l0[c] - step*(auxiliar.dy[c]);
                if(problem.l[c]<0.0)
                    problem.l[c] = 0.0;
            }
        }
        if( problem.CL>0 ){
            for( unsigned int c=0; c<problem.CL; ++c ){
                problem.mu[c] = problem.mu0[c] - step*(auxiliar.dylb[c]);
                if(problem.mu[c]<0.0)
                    problem.mu[c] = 0.0;
            }
        }
        if( problem.CU>0 ){
            for( unsigned int c=0; c<problem.CU; ++c ){
                problem.eta[c] = problem.eta0[c] - step*(auxiliar.dyub[c]);
                if(problem.eta[c]<0.0)
                    problem.eta[c] = 0.0;
            }
        }
        return;
    }
    
    /**
     * Revert the values of the Lagrange multipliers to their previous states
     */
    void revertMultipliers(  QPProblem& problem )
    {
        if( problem.CE>0 )
            memcpy( problem.leq, problem.leq0, (problem.CE)*sizeof(double) );
        if( problem.CI>0 )
            memcpy( problem.l, problem.l0, (problem.CI)*sizeof(double) );
        if( problem.CL>0 )
            memcpy( problem.mu, problem.mu0, (problem.CL)*sizeof(double) );
        if( problem.CU>0 )
            memcpy( problem.eta, problem.eta0, (problem.CU)*sizeof(double) );
        return;
    }

    /**
     * Forward the values of the Lagrange multipliers to their previous states
     */
    void forwardMultipliers(  QPProblem& problem )
    {
        if( problem.CE>0 )
            memcpy( problem.leq0, problem.leq, (problem.CE)*sizeof(double) );
        if( problem.CI>0 )
            memcpy( problem.l0, problem.l, (problem.CI)*sizeof(double) );
        if( problem.CL>0 )
            memcpy( problem.mu0, problem.mu, (problem.CL)*sizeof(double) );
        if( problem.CU>0 )
            memcpy( problem.eta0, problem.eta, (problem.CU)*sizeof(double) );
        return;
    }

    /**
     * A function to compute the solution of the primal problem
     * from the solution of the dual problem:
     *    x = -Q^(-1)*(f+C'*y)
     */
    void computePrimalSolution( const QPProblem& problem, QPAuxiliar& auxiliar )
    {
        // Compute (C')*y taking advantage of the sparsity:
        computeCty( problem, auxiliar );
        // Add to f:
        mataux::addMxArrays( auxiliar.resultNx1, problem.f, auxiliar.resultNx1, problem.N, 1 );
        // Pre-multiply by Qi^(-1):
        mataux::multiplyMxArrays( problem.Qi, auxiliar.resultNx1, problem.x,
                                      problem.N, problem.N, 1 ); // N x 1
        // Invert
        mataux::scalaropMxArray( problem.x, problem.N, 1, -1.0, mataux::MULTIPLY );
        return;
    }
    
    /**
     * A function to efficiently compute (C')*y taking advantage of
     * the potential sparsity of the Lagrange multipliers. The
     * output is stored in auxiliar.resultNx1
     */
    void computeCty( const QPProblem& problem, QPAuxiliar& auxiliar )
    {
        for( unsigned int n=0; n<problem.N; ++n ){ auxiliar.resultNx1[n] = 0.0; }
        if(problem.CE>0){
            mataux::multiplyMxArrays( problem.Aeqt, problem.leq, auxiliar.Ctyeq,
                                      problem.N, problem.CE, 1 ); // N x 1
            mataux::addMxArrays( auxiliar.resultNx1, auxiliar.Ctyeq, auxiliar.resultNx1, problem.N, 1 ); // N x 1
        }
        if(problem.CI>0){
            for( unsigned int n=0; n<problem.N; ++n ){auxiliar.Cty[n] = 0.0;}
            for( unsigned int c=0; c<problem.CI; ++c ){
                if( problem.l[c]>0.0 ){
                    for( unsigned int n=0; n<problem.N; ++n )
                        auxiliar.Cty[n] += (problem.l[c])*(problem.At[n+(problem.N)*c]);
                }
            }
            mataux::addMxArrays( auxiliar.resultNx1, auxiliar.Cty, auxiliar.resultNx1, problem.N, 1 ); // N x 1
        }
        if( problem.CL>0 ){ // If we actually have lower bounds
            // Compute (C')*y, which has size N x 1. We assume mu is typically
            // sparse, hence we add only the required terms. Besides, C' is
            // just -I_N
            for( unsigned int n=0; n<problem.N; ++n ){auxiliar.Ctylb[n] = 0.0;}
            for( unsigned int c=0; c<problem.CL; ++c ){ // CL=N
                if( problem.mu[c]>0.0 )
                    auxiliar.Ctylb[c] = -(problem.mu[c]);
            }
            mataux::addMxArrays( auxiliar.resultNx1, auxiliar.Ctylb, auxiliar.resultNx1, problem.N, 1 ); // N x 1
        }
        if( problem.CU>0 ){ // If we actually have upper bounds
            // Compute (C')*y, which has size N x 1. We assume eta is typically
            // sparse, hence we add only the required terms. Besides, C' is
            // just I_N
            for( unsigned int n=0; n<problem.N; ++n ){auxiliar.Ctyub[n] = 0.0;}
            for( unsigned int c=0; c<problem.CU; ++c ){ // CU=N
                if( problem.eta[c]>0.0 )
                    auxiliar.Ctyub[c] = problem.eta[c];
            }
            mataux::addMxArrays( auxiliar.resultNx1, auxiliar.Ctyub, auxiliar.resultNx1, problem.N, 1 ); // N x 1
        }
        return;
    }
    
    /**
     * Initalize all Lagrange multipliers to 0,
     * which corresponds to the unconstrained solution
     * of the primal problem
     */
    void initializeMultipliers( QPProblem& problem )
    {
        if(problem.CE>0){
            mataux::setValueMxArray(problem.leq,problem.CE,1,0.0);
            mataux::setValueMxArray(problem.leq0,problem.CE,1,0.0);
        }
        if(problem.CI>0){
            mataux::setValueMxArray(problem.l,problem.CI,1,0.0);
            mataux::setValueMxArray(problem.l0,problem.CI,1,0.0);
        }
        if(problem.CL>0){
            mataux::setValueMxArray(problem.mu,problem.CL,1,0.0);
            mataux::setValueMxArray(problem.mu0,problem.CL,1,0.0);
        }
        if(problem.CU>0){
            mataux::setValueMxArray(problem.eta,problem.CU,1,0.0);
            mataux::setValueMxArray(problem.eta0,problem.CU,1,0.0);
        }
        return;
    }
    
    /**
     * This function normalizes both the equality and the inequality
     * constraints so that all the rows of the constraints matrix have
     * RMS value 1. This may be crucial to obtain a proper scaling of
     * the gradient of the cost function.
     * NOTE: for the lower and upper bounds, the normalization is
     * already performed by construction
     */
    void normalizeConstraints( QPProblem& problem )
    {
        if( problem.CE>0 ){
            // Since Aeqt is stored transposed, we need
            // to normalize each of its CE columns:
            for( unsigned int c=0; c<problem.CE; ++c ){
                // Compute the RMS value of the c-th column:
                double rms = mataux::rmsMxNDArray( &(problem.Aeqt[c*(problem.N)]), problem.N );
                // Normalize the c-th column with this value:
                mataux::scalaropMxArray( &(problem.Aeqt[c*(problem.N)]), problem.N, 1, 1/rms, mataux::MULTIPLY );
                // Normalize the value of the constraint:
                problem.beq[c] /= rms;
            }
        }
        if( problem.CI>0 ){
            // Since At is stored transposed, we need
            // to normalize each of its CI columns:
            for( unsigned int c=0; c<problem.CI; ++c ){
                // Compute the RMS value of the c-th column:
                double rms = mataux::rmsMxNDArray( &(problem.At[c*(problem.N)]), problem.N );
                // Normalize the c-th column with this value:
                mataux::scalaropMxArray( &(problem.At[c*(problem.N)]), problem.N, 1, 1/rms, mataux::MULTIPLY );
                // Normalize the value of the constraint:
                problem.b[c] /= rms;
            }
        }
        return;
    }
    
    /**
     * This function computes an approximation to the initial step size
     * of the gradient descent based on the features of the quadratic 
     * problem (the size of the linear and quadratic terms of the cost
     * function)
     */
    void computeStep0( QPProblem& problem, const QPAuxiliar& auxiliar )
    {
        // It seems that:
        //     rms(g) / trace(Qi) / 1000,
        // with g the linear term of the dual problem, works fairly good
        problem.step = 0.0;
        if(problem.CE>0){
            double stepe = mataux::rmsMxNDArray( auxiliar.fdualeq, problem.CE );
            problem.step += stepe*stepe;
        }
        if(problem.CI>0){
            double stepi = mataux::rmsMxNDArray( auxiliar.fdual, problem.CI );
            problem.step += stepi*stepi;
        }
        if(problem.CL>0){
            double stepl = mataux::rmsMxNDArray( auxiliar.fduallb, problem.CL );
            problem.step += stepl*stepl;
        }
        if(problem.CU>0){
            double stepu = mataux::rmsMxNDArray( auxiliar.fdualub, problem.CU );
            problem.step += stepu*stepu;
        }
        problem.step  = 0.001*sqrt(problem.step);
        problem.step /= mataux::traceMxArray( problem.Qi, problem.N );
        return;
    }
    
    /**
     * The actual routine to solve the Quadratic Program
     */
    int solveQuadraticProgram( QPProblem& problem, QPAuxiliar& auxiliar, const QPParams& params )
    {
        int result = QPP_OK;
        double cost0; 
        double stepfct = 1.0;
        unsigned int ngoods = 0;
        // If required, normalize the constraints:
        if(params.normalize)
            normalizeConstraints( problem );
        // Initialize all Lagrange multipliers to 0. Note the equivalent
        // primal solution is the unconstrained one:
        initializeMultipliers( problem );
        // Compute the initial approach to the primal problem and the initial cost
        computePrimalSolution( problem, auxiliar ); // It is in problem.x
        cost0 = computeDualCost( problem, auxiliar );
        // In case there are no constraints, we are actually done:
        if(   (problem.CE==0) && (problem.CI==0) && (problem.CL==0) && (problem.CU==0)    )
            return QPP_UNCONSTRAINED;
        // Precompute all that can be precomputed:
        setCQi( problem, auxiliar );   // The four terms C*Qi
        setDualf( problem, auxiliar ); // The linear term of the dual cost
        // If necessary, estimate the initial step size to be used:
        if(params.computes0)
            computeStep0( problem, auxiliar );
        // Iterate:
        for( problem.iters=0; problem.iters<params.maxiters; ++problem.iters ){
            // Compute the gradient for the current value of the multipliers
            computeDualGradient( problem, auxiliar );
            // Step along the negative gradient direction (and project on
            // boundaries) until the cost is reduced:
            problem.cost = std::numeric_limits<double>::infinity();
            while(problem.cost>cost0){
                // Update the multipliers:
                gradientDescent( problem, auxiliar, (problem.step)*stepfct );
                // Re-compute the cost function:
                problem.cost = computeDualCost( problem, auxiliar );
                // Check if we succeeded to decrease the cost
                if( problem.cost<=cost0 ){ // Yes, we succeeded
                    if(ngoods<params.streak){++ngoods;}
                    if(ngoods==params.streak)
                        stepfct *= 2.0;
                    // Check if the relative drecrease in the cost is
                    // too small, then exit:
                    if( (cost0-problem.cost)/std::abs(cost0) < params.costtol ){ // The cost may be negative!
                        computePrimalSolution( problem, auxiliar );
                        return QPP_SMALLCOST;
                    }
                    else{ // Otherwise, update the multipliers and iterate
                        cost0 = problem.cost;
                        forwardMultipliers( problem );
                    }
                }
                else{ // We failed. Decrease the step and try again
                    stepfct *= 0.5;
                    ngoods = 0;
                    // We should check if the step is too small, then
                    // exit:
                    if( stepfct<params.steptol ){
                        revertMultipliers( problem ); // y = y0
                        computePrimalSolution( problem, auxiliar );
                        return QPP_SMALLSTEP;
                    }
                }
            }
        }
        // If we end up here, the maximum number of iterations was reached.
        computePrimalSolution( problem, auxiliar );
        return QPP_DIDNT_CONVERGE;
    }
    
} // namespace dmriqpp

#endif // #ifndef _dmriquadprog_cxx_
