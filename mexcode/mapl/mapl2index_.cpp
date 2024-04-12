/*==========================================================
 * mapl2index_.c
 *
 * This is a core function to mapl2index, and should only
 * be called from therein
 *
 * This is a MEX-file for MATLAB.
 * Copyright 2024- Antonio Tristán Vega
 *
 *========================================================*/

#include "mex.h"
#include "matrix.h"
#include "math.h"
#include <cmath>
#include <algorithm>
#include "../mathsmex/matrixCalculus.h"
#include "../threads/threadHelper.h"
#include "hermitePols.h"
#include "maplMaths.h"
#include "../mathsmex/sanityCheckDTI.h"

typedef struct MAPLIOData{
    // Input:
    BufferType mapl;
    BufferType dti;
    // Output:
    BufferType index;
    
} MAPLIOData;

typedef struct MAPLParameters{
    unsigned char type; // The particular index we ask for
    unsigned int Nmax;  // The maximum degree of Hermite polynomials
    double tau;         // The effective diffusion time
    ElementType ADC0;   // Free-water diffusivity
} MAPLParameters;

class ThArgs : public DMRIThreader
{
public:
    MAPLIOData* io;
    MAPLParameters* params;
};

ElementType computeRTOP( const BufferType mapl, const BufferType Bnnn, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax );

ElementType computeRTAP( const BufferType mapl, const BufferType Bnnn, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax );

ElementType computeRTPP( const BufferType mapl, const BufferType Bnnn, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax );

ElementType computePADTI( const ElementType& ux,
                          const ElementType& uy, const ElementType& uz );

ElementType computeNG( const BufferType mapl, const unsigned int Nmax );

ElementType computeE0( const BufferType mapl, const BufferType Bnnn, const unsigned int Nmax );

ElementType computeMSD( const BufferType mapl, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax );

ElementType computeQIV( const BufferType mapl, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax );

ElementType computeLEnergy( const BufferType mapl, const BufferType U, 
                            BufferType Uc, const unsigned int Nmax );

ElementType computeu0( const ElementType& ux, const ElementType& uy, const ElementType& uz );

int polyRoots( const double& a0, const double& a1, const double& a2, double* WR, double* WI );

THFCNRET mapl2index_process_fcn( void* );

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    /** INPUTS (No error checking is performed):
     * 
     * prhs[0]:  mapl, the MAPL cofficients, M x N
     * prhs[1]:  dti, the tensor fit of that signal, 6 x N
     * prhs[2]:  ADC0, the free-water diffusivity at body temperature, 1 x 1
     * prhs[3]:  tau, the effective diffusion time in seconds, 1x1
     * prhs[4]:  type, the index to be computed (char)
     * prhs[5],  maxthreads, the maximum number of threads, 1x1
     *
     *  OUTPUTS:
     *
     * plhs[0]: atti, G x N, the reconstructed E(q)
     * 
     */
    //=======================================================================================
    /** Make sure this function is called in a "controlable" way*/
    mxArray* callstack[1];
    mexCallMATLAB(1, callstack, 0, NULL, "dbstack");
    char *callerFunc = mxArrayToString( mxGetField(callstack[0], 0, "name") );
    if( callerFunc==(char*)NULL )
        mexErrMsgIdAndTxt("MyToolbox:mapl2index_:callstack","This function should only be called from mapl2index");
    else if( strcmp(callerFunc,"mapl2index") )
        mexErrMsgIdAndTxt("MyToolbox:mapl2index_:callstack","This function should only be called from mapl2index");
    //=======================================================================================
    if(nrhs!=6)
        mexErrMsgIdAndTxt("MyToolbox:mapl2index_:nrhs","Exactly 5 input arguments are required");
    //=======================================================================================
    SizeType N = mxGetN(prhs[0]); // The number of voxels to process
    SizeType M = mxGetM(prhs[0]); // The number of basis functions
    unsigned int Nmax = 0;
    while( mapl::numBasisFunctions(Nmax) < M )
        Nmax += 2;
    if( mapl::numBasisFunctions(Nmax) != M )
        mexErrMsgIdAndTxt("MyToolbox:mapl2index_:dims","Input 0 doesn't seem a MAPL volume");
    //=======================================================================================
    MAPLParameters params;
    params.Nmax = Nmax;
    params.ADC0 = mxGetScalar(prhs[2]);
    params.tau = mxGetScalar(prhs[3]);
    params.type = mxGetChars(prhs[4])[0];    
    //=======================================================================================
    MAPLIOData io;
    // ------ Inputs
    io.mapl = mxGetDoubles(prhs[0]);
    io.dti = mxGetDoubles(prhs[1]);
    // ------ Outputs
    if(nlhs>0){
        plhs[0] = mxCreateDoubleMatrix( 1, N, mxREAL );
        io.index = mxGetDoubles(plhs[0]);
    }
    else
        return;
    //=======================================================================================
    unsigned int maxthreads = get_number_of_threads( (unsigned int)mxGetScalar(prhs[5]) );
    unsigned long chunksz = N/maxthreads/20;
    if(chunksz<1)
        chunksz = 1;
    //=======================================================================================
    // Use the helper class to pass arguments. Inherited values:
    ThArgs threader;
    threader.setProcessSize( N, chunksz );
    // Own values:
    threader.io         = &io;
    threader.params     = &params;
    //=======================================================================================
    // Do the threaded job:
    threader.threadedProcess( maxthreads, mapl2index_process_fcn );
    //=======================================================================================
    return;
}

THFCNRET mapl2index_process_fcn( void* inargs )
{
    // Retrieve the structure with all the parameters.
    ThArgs* args = (ThArgs*)inargs;
    MAPLIOData* io = args->io;
    const MAPLParameters* params = args->params;
    const unsigned int Nmax = params->Nmax;
    const unsigned char type = params->type;
    const double ADC0 = params->ADC0;
    const double tau = params->tau;
    const SizeType M = mapl::numBasisFunctions( Nmax );
        
    // Note: this call is crucial so that subsequent calls to
    // Lapack/BLAS won't create their own threads that blow up
    // the total amount of threads putting down the overall
    // peformance.
    unsigned int blas_threads = blas_num_threads(1);
    // -----------------------------------------------------------------------------
    // Auxiliary buffers to compute DTI spectrum:
    ElementType dti[6];
    ElementType eigval[3];
    ElementType eigvec[9];
    // -----------------------------------------------------------------------------
    // The values of the phi functions at (0,0,0)
    BufferType Bnnn = new ElementType[M];
    mapl::computePhi0Values(  Bnnn, Nmax );
    // -----------------------------------------------------------------------------
    // Prepare buffers for the regularization matrix,
    // U:             M x M
    // nx, ny, nz:    M x 1
    // Snm, Tnm, Unm: (Nmax+1) x (Nmax+1)
    BufferType U, Snm, Tnm, Unm, Uc;
    unsigned int *nx, *ny, *nz;
    if( type=='l' ){
        mapl::allocateRegularizationMatrix( U, Snm, Tnm, Unm,
                                            nx, ny, nz, Nmax );
        Uc = new ElementType[M];
        // nx, ny, and nz are data-independent, hence they can be pre-computed here:
        mapl::computenxyz( nx, ny, nz, Nmax );
        // The same goes for Snm, Tnmm, and Unm
        mapl::computeSTUnm( Snm, Tnm, Unm, Nmax );
    }
    // -----------------------------------------------------------------------------
    // That's all we can pre-compute. From now on, we will perform 
    // voxel-especific operations:
    IndexType start = 0; // First index of data blocks to process
    IndexType end   = 0; // Post-last index of data blocks to process
    do{
        // Claim a new block of data to process within this thread:
        args->claimNewBlock( &start, &end );
        // While there is still data not processed by any thread...
        for(  IndexType i=start; i<end; ++i ){ // Process the current chunk
            // ---------------------------------------------------------------------
            // 1- Check the tensor model and compute eigenvalues and eigenvectors
            //    to shape the transformed (anatomical) space (will use Lapack's
            //    dsyev). This is necessary only if the computation of this index
            //    depends on ux, uy, uz
            double ux, uy, uz;
            switch(type)
            {
                case 'o': // RTOP
                case 'a': // RTAP
                case 'p': // RTPP
                case 'd': // PA_DTI
                case 'm': // MSD
                case 'q': // QIV
                case 'l': // Laplacian energy
                case 'u': // u0, most similar isotropic propagator
                    memcpy( dti, &(io->dti[6*i]), 6*sizeof(ElementType) );
                    dtisc::sanityCheckDTI( (BufferType)dti, (BufferType)eigval, (BufferType)eigvec, args->params->ADC0, 'D' );
                    ux = sqrt(2.0*eigval[0]*tau);
                    uy = sqrt(2.0*eigval[1]*tau);
                    uz = sqrt(2.0*eigval[2]*tau);
                    break;
                case 'g': // NG
                case 'e': // E(0)
                    break;
                default:
                    break;
            }
            // ---------------------------------------------------------------------
            // 2- Compute the required index:
            double index;
            switch(type)
            {
                case 'o': // RTOP
                    index = computeRTOP( &(io->mapl[i*M]), Bnnn, ux, uy, uz, Nmax );
                    break;
                case 'a': // RTAP
                    index = computeRTAP( &(io->mapl[i*M]), Bnnn, ux, uy, uz, Nmax );
                    break;
                case 'p': // RTPP
                    index = computeRTPP( &(io->mapl[i*M]), Bnnn, ux, uy, uz, Nmax );
                    break;
                case 'd': // PA_DTI
                    index = computePADTI( ux, uy, uz );
                    break;
                case 'g': // NG
                    index = computeNG( &(io->mapl[i*M]), Nmax );
                    break;
                case 'e': // E(0)
                    index = computeE0( &(io->mapl[i*M]), Bnnn, Nmax );
                    break;
                case 'm': // MSD
                    index = computeMSD( &(io->mapl[i*M]), ux, uy, uz, Nmax );
                    break;
                case 'q': // QIV
                    index = computeQIV( &(io->mapl[i*M]), ux, uy, uz, Nmax );
                    break;
                case 'l': // Laplacian energy
                    mapl::computeRegularizationMatrix( U, Snm, Tnm, Unm,
                                       nx, ny, nz, ux, uy, uz, Nmax );
                    index = computeLEnergy( &(io->mapl[i*M]), U, Uc, Nmax );
                    break;
                case 'u':
                    index = computeu0( ux, uy, uz );
                    break;
                default:
                    break;
            }
            // ---------------------------------------------------------------------
            // 3- Assign the output
            io->index[i] = index;
            // ---------------------------------------------------------------------
        }
    }
    while( start < args->getN() );

    // Revert BLAS threads usage to its default:
    blas_num_threads(blas_threads);

    // Free memory previously allocated
    delete[] Bnnn;
    if( type=='l' ){
        mapl::destroyRegularizationMatrix( U, Snm, Tnm, Unm, nx, ny, nz );
        delete[] Uc;
    }
    
    return (THFCNRET)NULL;
}


ElementType computeRTOP( const BufferType mapl, const BufferType Bnnn, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax )
{
    int sign = 1;
    ElementType rtop = 0.0;
    unsigned long pos = 0;
    for( unsigned int L=0; L<=Nmax; L+=2 ){ // Loop through even integers
        for( unsigned int n1=0; n1<=L; n1+=2 ){ // The Bnnn are non-null only if all n1, n2, n3 are even
            for( unsigned int n2=0; n2<=L-n1; n2+=2 ){ // The Bnnn are non-null only if all n1, n2, n3 are even
                unsigned int n3 = L-n1-n2;
                rtop += mapl[pos]*Bnnn[pos] * sign;
                ++pos;
            }
        }
        sign = -sign;
    }
    rtop /= (sqrt(8.0*PI*PI*PI)*ux*uy*uz);
    return rtop;
}

ElementType computeRTAP( const BufferType mapl, const BufferType Bnnn, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax )
{
    int sign;
    unsigned long pos = 0;
    ElementType rtap = 0.0;
    for( unsigned int L=0; L<=Nmax; L+=2 ){ // Loop through even integers
        for( unsigned int n1=0; n1<=L; n1+=2 ){ // The Bnnn are non-null only if all n1, n2, n3 are even
            for( unsigned int n2=0; n2<=L-n1; n2+=2 ){ // The Bnnn are non-null only if all n1, n2, n3 are even
                unsigned int n3 = L-n1-n2;
                unsigned int test = (n2+n3)/2;
                if(2*(test/2)==test)
                    sign = 1;
                else
                    sign = -1;
                rtap += mapl[pos]*Bnnn[pos] * sign;
                ++pos;
            }
        }
    }
    rtap /= (2.0*PI*uy*uz);
    return rtap;
}

ElementType computeRTPP( const BufferType mapl, const BufferType Bnnn, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax )
{
    int sign;
    unsigned long pos = 0;
    ElementType rtpp = 0.0;
    for( unsigned int L=0; L<=Nmax; L+=2 ){ // Loop through even integers
        for( unsigned int n1=0; n1<=L; n1+=2 ){ // The Bnnn are non-null only if all n1, n2, n3 are even
            for( unsigned int n2=0; n2<=L-n1; n2+=2 ){ // The Bnnn are non-null only if all n1, n2, n3 are even
                unsigned int test = n1/2;
                if(2*(test/2)==test)
                    sign = 1;
                else
                    sign = -1;
                rtpp += mapl[pos]*Bnnn[pos] * sign;
                ++pos;
            }
        }
    }
    rtpp /= ((2.0*PI)*ux);
    return rtpp;
}

ElementType computePADTI( const ElementType& ux,
                          const ElementType& uy, const ElementType& uz )
{
    ElementType padti;
    ElementType u0 = computeu0( ux, uy, uz );
    if(u0<0.0)
        return NAN;
    else{
        padti = 8*(u0*u0*u0)*ux*uy*uz;
        padti /= (ux*ux+u0*u0);
        padti /= (uy*uy+u0*u0);
        padti /= (uz*uz+u0*u0);
    }
    if(padti>1.0)
        padti = 1.0;
    padti = sqrt(1-padti);
    ElementType t = padti;
    padti  = pow(t,3*0.4);
    padti /= ( 1.0 - 3.0*pow(t,0.4) + 3.0*pow(t,2*0.4) );
    return padti;
}

ElementType computeNG( const BufferType mapl, const unsigned int Nmax )
{
    const SizeType M = mapl::numBasisFunctions( Nmax );
    ElementType ng = 0.0;
    for( IndexType m=0; m<M; ++m )
        ng += mapl[m]*mapl[m];
    ng = mapl[0]*mapl[0] / ng;
    ng = sqrt(1.0-ng);
    return ng;
}

ElementType computeE0( const BufferType mapl, const BufferType Bnnn, const unsigned int Nmax )
{
    const SizeType M = mapl::numBasisFunctions( Nmax );
    ElementType e0 = 0.0;
    for( IndexType m=0; m<M; ++m )
        e0 += mapl[m]*Bnnn[m];
    return e0;
}

ElementType computeMSD( const BufferType mapl, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax )
{
    ElementType msd = 0.0;
    unsigned long pos = 0;
    for( unsigned int L=0; L<=Nmax; L+=2 ){ // Loop through even integers
        for( unsigned int n1=0; n1<=L; ++n1 ){
            for( unsigned int n2=0; n2<=L-n1; ++n2 ){
                unsigned int n3 = L-n1-n2;
                if( (n1==2*(n1/2)) && (n2==2*(n2/2)) && (n3==2*(n3/2)) ){
                    unsigned int ns = L/2;
                    ElementType num = (1+2*n1)*ux*ux + (1+2*n2)*uy*uy + (1+2*n3)*uz*uz;
                    if( ns != 2*(ns/2) )
                        num = -num;
                    ElementType den = -(ElementType)L*log(2.0) + lgamma((ElementType)n1+1) 
                       + lgamma((ElementType)n2+1) + lgamma((ElementType)n3+1);
                    den *= 0.5;
                    den += lgamma(0.5*(1-(ElementType)n1)) 
                       + lgamma(0.5*(1-(ElementType)n2)) + lgamma(0.5*(1-(ElementType)n3));
                    den = exp(-den);
                    msd += num*den*mapl[pos];
                }
                ++pos;
            }
        }
    }
    return sqrt(PI*PI*PI)*msd;
}

ElementType computeQIV( const BufferType mapl, const ElementType& ux, 
                    const ElementType& uy, const ElementType& uz, const unsigned int Nmax )
{
    ElementType qiv = 0.0;
    unsigned long pos = 0;
    for( unsigned int L=0; L<=Nmax; L+=2 ){ // Loop through even integers
        for( unsigned int n1=0; n1<=L; ++n1 ){
            for( unsigned int n2=0; n2<=L-n1; ++n2 ){
                unsigned int n3 = L-n1-n2;
                if( (n1==2*(n1/2)) && (n2==2*(n2/2)) && (n3==2*(n3/2)) ){
                    ElementType num = lgamma((ElementType)n1+1) + lgamma((ElementType)n2+1) + lgamma((ElementType)n3+1);
                    num *= 0.5;
                    num += lgamma(0.5*(1-(ElementType)n1)) 
                       + lgamma(0.5*(1-(ElementType)n2)) + lgamma(0.5*(1-(ElementType)n3));
                    ElementType den = 0.5*((ElementType)L-1)*log(2.0);
                    num = exp(num-den);
                    den = (1+2*n1)*uy*uy*uz*uz + (1+2*n2)*ux*ux*uz*uz + (1+2*n3)*ux*ux*uy*uy;
                    qiv += (num/den)*mapl[pos];
                }
                ++pos;
            }
        }
    }
    return (8.0*PI*PI)*(ux*ux*ux)*(uy*uy*uy)*(uz*uz*uz)*qiv;
}

ElementType computeLEnergy( const BufferType mapl, const BufferType U, 
                            BufferType Uc, const unsigned int Nmax )
{
    const SizeType M = mapl::numBasisFunctions( Nmax );
    ElementType lenerg = 0.0;
    mataux::multiplyMxArrays( U, mapl, Uc, M, M, 1 ); // M x 1
    mataux::multiplyMxArrays( mapl, Uc, &lenerg, 1, M, 1 ); // 1 x 1
    return lenerg;
}

ElementType computeu0( const ElementType& ux, const ElementType& uy, const ElementType& uz )
{
    ElementType u0;
    double X = ux*ux;
    double Y = uy*uy;
    double Z = uz*uz;
    // Normalize the polynomial:
    double a3 = -3.0;
    double a0 = 3.0*X*Y*Z / a3;
    double a1 = (X*Y+X*Z+Y*Z) / a3;
    double a2 = -(X+Y+Z) / a3;
    double WR[3];
    double WI[3];
    int result = polyRoots( a0, a1, a2, (double*)WR, (double*)WI );
    if( result!=0 )
        return (ElementType)result;
    // If we reach this point, the roots of the polynomial
    // are stored in WR (real part) and WI (imaginary part)
    // Like in the dipy package, we keep the solution with
    // the largest real part:
    u0 = std::max( WR[0], std::max(WR[1],WR[2]) );
    if(u0<=0.0)
        return -3.0;
    else
        u0 = sqrt(u0);    
    return u0;
}

int polyRoots( const double& a0, const double& a1, const double& a2, double* WR, double* WI )
{
    // Compute the companion matrix for the polynomial, which is 3x3:
    double companion[9];
    companion[0] = companion[1] = companion[4] = companion[6] = 0.0;
    companion[3] = companion[7] = 1.0;
    companion[2] = -a0;
    companion[5] = -a1;
    companion[8] = -a2;
    // The eigenvalues of the companion matrix are the roots of the
    // polynomial. We first reduce the matrix to its Hessenberg form:
    double tau[2];
    ptrdiff_t N = 3;
    ptrdiff_t O = 1;
    ptrdiff_t L = 9;
    ptrdiff_t info = 0;
    double work[9];
    dgehrd( &N, &O, &N, (double*)companion, &N, 
            (double*)tau, (double*)work, &L, &info );
    if(info==0){ // Now we can actually compute the eigenvalues
        dhseqr( "E", "N", &N, &O, &N, (double*)companion, &N,
                (double*)WR, (double*)WI, (double*)NULL, &N,
                (double*)work, &L, &info );
        if(info!=0)
            return -2;
    }
    else
        return -1;
    return 0;
}
