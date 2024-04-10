/**
 * This header file implements a basic Quadratic Program solver based
 * on a projected gradient descent solution of the dual problem. It is
 * not intended as a general purpose method, but as an ad hoc sover for
 * problems related to Diffussion MRI fitting of signals.
 * 
 * In this kind of problems, the number of variables is typically 
 * moderated, one has only one equatity constraint, and either 
 * A*x >= bin or x >= lb inequality constraints, but not both at the
 * same time. Besides, it is assumed that most of the inequality
 * constraints will be NOT active
 * 
 * The most general form of the problem considered is:
 * 
 *      min_x    0.5 x'*Q*x + f'*x
 *      s.t.:    Aeq*x = beq
 *               A*x <= b
 *               lb <= x <= ub
 * 
 * with Q symmetric and positive definite (hence, Cholesky-
 * factorizable). The dual formulation of this problem is:
 * 
 *    min_{leq,l,mu,eta}
 *        0.5 y' * C*Q^(-1)*C' * y + (f'*Q^(-1)*C'+u') * y
 *    s.t.:
 *        l   >= 0
 *        mu  >= 0
 *        eta >= 0
 * 
 * where y = [leq',l',mu',eta']' is the vector of Lagrange multipliers
 * (unknowns of the dual problem); C = [Aeq',A',-I_N',I_N]';
 * u = [beq',b',-lb',ub']'.
 * 
 * NOTE: in dMRI problems: Aeq has typically one row, so that leq
 * reduces to a scalar; ub (hence eta) is typically not present;
 * l and mu are typically highly sparse, since many constraints will
 * not be active. All this means that the product C'*y can be effciently
 * computed by working only with non-null components of y.
 */

#ifndef _dmriquadprog_h_
#define _dmriquadprog_h_

#include <cstring>
#include <cmath>
#include <limits>
#include "../mathsmex/matrixCalculus.h"
#include "../mathsmex/mexToMathsTypes.h"

namespace dmriqpp
{

#define QPP_SMALLCOST 3
#define QPP_SMALLSTEP 2
#define QPP_UNCONSTRAINED 1
#define QPP_OK 0
#define QPP_BAD_ARGS -1
#define QPP_DIDNT_CONVERGE -2
    
    /**
     * A structure to store the input data to the problem:
     */
    typedef struct QPProblem{
        // Problem dimensions
        unsigned int N;  // Problem size; Q is N x N
        unsigned int CE; // Number of equality constraints, Aeq is CE x N
        unsigned int CI; // Number of inequality constraints, A is CI x N
        unsigned int CL; // Number of lower bounds, either N or 0
        unsigned int CU; // Number of upper bounds, either N or 0
        // Problem input data:
        double* Qi;   // N x N, the inverse of Q
        double* f;    // N x 1
        double* Aeqt; // N x CE
        double* beq;  // CE x 1
        double* At;   // N x CI
        double* b;    // CI x 1
        double* lb;   // CL x 1
        double* ub;   // CU x 1
        // Problem solutions:
        double* x;   // N x 1, solution of the primal problem
        double* leq; // CE x 1, Lagrange multiplier of Aeq*x = beq
        double* l;   // CI x 1, Lagrange multiplier of A*x = b
        double* mu;  // CL x 1, Lagrange multiplier of lb <= x
        double* eta; // CU x 1, Lagrange multiplier of x <= ub
        // Problem solutions in the previous iteration:
        double* leq0; // CE x 1, Lagrange multiplier of Aeq*x = beq
        double* l0;   // CI x 1, Lagrange multiplier of A*x = b
        double* mu0;  // CL x 1, Lagrange multiplier of lb <= x
        double* eta0; // CU x 1, Lagrange multiplier of x <= ub
        // Scalar variables:
        unsigned int iters;
        double cost;
        double step;
    } QPProblem;
    
    /**
     * A structure to store intermediate results of the 
     * algorithm:
     */
    typedef struct QPAuxiliar{
        // The linear term of the quadratic cost function
        // in the dual problem is: (f'*Q^(-1)*C'+u')'. We
        // will store it in four separated buffers for 
        // each set of constraints:
        double* fdualeq; // CE x 1
        double* fdual;   // CI x 1
        double* fduallb; // CL x 1
        double* fdualub; // CU x 1
        // The vector to store the gradient of the dual
        // cost function. It is splitted in four parts:
        double* dyeq; // CE x 1
        double* dy;   // CI x 1
        double* dylb; // CL x 1
        double* dyub; // CU x 1
        // These are four buffers to store the pre-
        // computed result of (C*Q^(-1))', which is
        // necessary to compute the gradient of the
        // dual cost function. It is stored in four
        // blcoks, each of size Cr x N
        double* CQieqt; // N x CE
        double* CQit;   // N x CI
        double* CQilbt; // N x CL
        double* CQiubt; // N x CU
        // These are auxiliar buffers to compute matrix
        // products with size N x 1:
        double* resultNx1; // N x 1
        double* Ctyeq;     // N x 1
        double* Cty;       // N x 1
        double* Ctylb;     // N x 1
        double* Ctyub;     // N x 1
        // These are auxiliar buffers to compute matrix
        // products with size Cr x 1:
        double* ueq;  // CE x 1
        double* u;    // CI x 1
        double* ulb;  // CL x 1
        double* uub;  // CU x 1
    } QPAuxiliar;
    
    typedef struct QPParams{
        // Maximum number of iterations:
        unsigned long maxiters = 1000;
        // Every time the number of successful iterations
        // in a row is greater or equal than "streak", the
        // step of the gradient descent is doubled (it is
        // halved upon every unsuccessful iteration):
        unsigned int streak = 2;
        // We assume the algorithm has converged if either
        // the step size is smaller than steptol times the
        // initial one or if the relative change in the cost
        // function after a successful iteration is smaller
        // than costtol:
        double steptol = 1.0e-12;
        double costtol = 1.0e-12;
        // If normalize is true, them all the constraints
        // are normalized to get a RMS value of 1. This may
        // become crucial to avoid issues with badly scaled
        // gradients:
        bool normalize = true;
        // If this is set true, an approximation to the initial
        // step size of the gradient descent is computed based
        // on the features of the quadratic problem. Otherwise,
        // it is assumed that problem.step has been externally
        // fixed:
        bool computes0 = true;
    } QPParams;
    
    void allocateDouble( double*& buffer, const unsigned int M, const unsigned int N );
    
    void freeDouble( double*& buffer  );
    
    int allocateQPProblem( QPProblem& problem,
                            const unsigned int n, const unsigned int ce,
                            const unsigned int ci, const unsigned int cl,
                            const unsigned int cu );
    
    void freeQPProblem( QPProblem& problem );
    
    void setQi( QPProblem& problem, const double* Qi );
    
    void setf( QPProblem& problem, const double* f );
    
    void setAeq( QPProblem& problem, const double* Aeq );
    
    void setAeqt( QPProblem& problem, const double* Aeqt );
    
    void setbeq( QPProblem& problem, const double* beq );
    
    void setA( QPProblem& problem, const double* A );
    
    void setAt( QPProblem& problem, const double* At );
    
    void setb( QPProblem& problem, const double* b );
    
    void setlb( QPProblem& problem, const double* lb );
    
    void setub( QPProblem& problem, const double* ub );
    
    int allocateQPAuxiliar( const QPProblem& problem, 
                            const QPParams& params, QPAuxiliar& auxiliar );
    
    void freeQPAuxiliar( QPAuxiliar& auxiliar );
    
    void setDualf( const QPProblem& problem, QPAuxiliar& auxiliar );
    
    void setCQi( const QPProblem& problem, QPAuxiliar& auxiliar );
    
    void computeDualGradient( const QPProblem& problem, QPAuxiliar& auxiliar );
    
    double computeDualCost( const QPProblem& problem, QPAuxiliar& auxiliar );

    void gradientDescent( QPProblem& problem, const QPAuxiliar& auxiliar, const double& step );

    void revertMultipliers(  QPProblem& problem );

    void forwardMultipliers(  QPProblem& problem );
    
    void computePrimalSolution( const QPProblem& problem, QPAuxiliar& auxiliar );
    
    void computeCty( const QPProblem& problem, QPAuxiliar& auxiliar );
    
    void initializeMultipliers( QPProblem& problem );
    
    void normalizeConstraints( QPProblem& problem );
    
    void computeStep0( QPProblem& problem, const QPAuxiliar& auxiliar );
    
    int solveQuadraticProgram( QPProblem& problem, QPAuxiliar& auxiliar, const QPParams& params );
    
} // namespace dmriqpp

#endif // #ifndef _dmriquadprog_h_
