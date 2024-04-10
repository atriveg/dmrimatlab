#ifndef _sanityCheckDTI_h_
#define _sanityCheckDTI_h_

#include "../mathsmex/matrixCalculus.h"
#include <cmath>

namespace dtisc{
    
    void sanityCheckDTI( double* dti, double* eigval, double* eigvec, const double& ADC0, const char order );

}

#endif
