#ifndef _sanityCheckDTI_cxx_
#define _sanityCheckDTI_cxx_

#include "sanityCheckDTI.h"

namespace dtisc{
    
/**
 * Given a Diffusion tensor represented by its unique 6 components,
 * compute its eigenvalues and eigenvectors. Then, check the eigenvalues
 * to ensure they are within a valid range. Finally, sort them in
 * ascending or descenfing order
*/
void sanityCheckDTI( double* dti, double* eigval, double* eigvec, const double& ADC0, const char order )
{
    // Use Lapack's dspev:
    const ptrdiff_t dim = 3;
    ptrdiff_t info = 0;
    double work[9];
    dspev( "V", "L", &dim, dti, eigval, eigvec, &dim, work, &info );
    if(info!=0){
        // Computation failed for some reason. Fill eigenvalues with
        // free-water value:
        eigval[0] = eigval[1] = eigval[2] = ADC0;
        eigvec[0] = eigvec[4] = eigvec[8] = 1.0f;
        eigvec[1] = eigvec[2] = eigvec[3] = eigvec[5] = eigvec[6] =  eigvec[7] = 0.0f;
    }
    else{
        // Sanity checks...
        for( unsigned int p=0; p<3; ++p ){
            eigval[p] = std::abs(eigval[p]);
            if( eigval[p]<ADC0/100 )
                eigval[p] = ADC0/100;
            if( eigval[p]>ADC0 )
                eigval[p] = ADC0;
        }
        // Ordering (descending order)
        double eval;
        double evec[3];
        for( unsigned int p=0; p<2; ++p ){
            for( unsigned int q=p+1; q<3; ++q ){
                bool swap;
                if(order=='D')
                    swap = (eigval[p]<eigval[q]);
                else
                    swap = (eigval[p]>eigval[q]);
                if(swap){ // Must swap
                    eval = eigval[p];
                    evec[0] = eigvec[3*p+0];
                    evec[1] = eigvec[3*p+1];
                    evec[2] = eigvec[3*p+2];
                    eigval[p] = eigval[q];
                    eigval[q] = eval;
                    eigvec[3*p+0] = eigvec[3*q+0];
                    eigvec[3*p+1] = eigvec[3*q+1];
                    eigvec[3*p+2] = eigvec[3*q+2];
                    eigvec[3*q+0] = evec[0];
                    eigvec[3*q+1] = evec[1];
                    eigvec[3*q+2] = evec[2];
                }
            }
        }
        // Make sure eigvec is a rotation matrix
        // by making e1 = e2 x e3
        eigvec[0] = eigvec[4]*eigvec[8]-eigvec[5]*eigvec[7];
        eigvec[1] = eigvec[5]*eigvec[6]-eigvec[3]*eigvec[8];
        eigvec[2] = eigvec[3]*eigvec[7]-eigvec[4]*eigvec[6];
    }
    return;
}   

} // End namespace dtisc

#endif
