/** 
 * This is an implementation of the Gaussian Hypergeometric function 2F1
 * taken from:
 *
 *    Andrew Horchler (2020). 
 *    hypergeomq (https://www.github.com/horchler/hypergeomq), 
 *    GitHub. Retrieved March 18, 2020.
 *    Revision: 1.0, 5-19-14
 *
 * Since the algorithm converges only for abs(z)\<1, for other (real) values
 * of z we have to rely on the linear transformation formulas described in:
 *    
 *    N. Michel and M.V. Stoitsov
 *    "Fast computation of the Gauss hypergeometric function with all its
 *     parameters complex with application to the Pöschl-Teller-Ginocchio 
 *     potential wave functions"
 *    arXiv:0708.0116v2 (2007),
 *
 * eqs. (8)-(13). Note these expressions are not well defined for certain
 * combinations of parameters leading to evaluations of the Gamma function
 * at negative, integer values. However, for the particular purpose of this
 * toolbox, we can work around these situations.
 */

#ifndef _hypergeom2F1_h_
#define _hypergeom2F1_h_

#include <cmath>
#include "mexToMathsTypes.h"

namespace hypegeo
{

/**
 * Return values on different error conditions:
 */
#define H2F1_SUCCESS 0
#define H2F1_BADARGS -1
#define H2F1_NONCONVERGENT -2
#define H2F1_GAMMAUNDEF -3
#define H2F1_UNDEF -4
    
#ifndef PI
#define PI 3.14159265358979323846
#endif

    const double HYPEF21AZTH = 1.0 - 10*__eps__;
    
    /** The original implementation as provided (in Matlab code) by 
        Andrew Horchler for abs(z)<1 */
    int hyperGeom2F1_azl1( const double&, const double&, const double&, const double&, double& );
    
    /** The general function using the linear transformation formulas
        described by Michel and Stoitsov, valid for "almost all" z */
    int hyperGeom2F1( const double&, const double&, const double&, const double&, double& );
    
    /** A especialized function for the case we are interested in,
        when a=1/2, b=g/2, c=3/2 */
    int hyperGeom2F1( const double&, const double&, double& );
    
    bool isnegint(const double& );
    
    bool badCombination( const double&, const double&, const double& );
    
} // End namespace shmaths

#endif //_sphericalHarmonics_h_
