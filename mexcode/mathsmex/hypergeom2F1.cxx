#ifndef _hypergeom2F1_cxx
#define _hypergeom2F1_cxx

#include "hypergeom2F1.h"
#include "math.h"
#include <limits>

namespace hypegeo
{
    /** The original implementation as provided (in Matlab code) by 
        Andrew Horchler for abs(z)\<1 . */
    int hyperGeom2F1_azl1( const double& a, const double& b, const double& c, const double& z, double& f )
    {
        if( abs(z)>=1.0 ){ // This approximation is not defined in this case
            f = NAN;
            return H2F1_BADARGS;
        }
        
        // Tolerance to assess convergence:
        const double tol = __eps__;         // DOUBLE_EPS = 2.2204e-16
        const unsigned int itermax = 32768; // 2^15
        double y, yi;
        unsigned int i;
        
        if( a==1.0 ){
            if( b==c ){
                yi = z;
                y  = 1+yi;
                for( i=1; i<itermax; ++i ){
                    yi *= z;
                    y  += yi;
                    if( abs(yi) < tol )
                        break;
                }
            }
            else{ // General case
                if( abs(c) < 10*__eps__ ){
                    f = std::numeric_limits<double>::infinity();
                    return H2F1_BADARGS;
                }
                yi = b*z/c;
                y  = 1+yi;
                for( i=1; i<itermax; ++i ){
                    yi *= (b+i)*z/(c+i);
                    y  += yi;
                    if( abs(yi) < tol )
                        break;
                }
            }
        }
        else if( b==1.0 ){
            if( a==c ){
                yi = z;
                y  = 1+yi;
                for( i=1; i<itermax; ++i ){
                    yi *= z;
                    y  += yi;
                    if( abs(yi) < tol )
                        break;
                }
            }
            else{
                if( abs(c) < 10*__eps__ ){
                    f = std::numeric_limits<double>::infinity();
                    return H2F1_BADARGS;
                }
                yi = a*z/c;
                y  = 1+yi;
                for( i=1; i<itermax; ++i ){
                    yi *= (a+i)*z/(c+i);
                    y  += yi;
                    if( abs(yi) < tol )
                        break;
                }
            }
        }
        else if( a==c ){
            yi = b*z;
            y  = 1+yi;
            for( i=1; i<itermax; ++i ){
                yi *= (b+i)*z/(i+1);
                y  += yi;
                if( abs(yi) < tol )
                    break;
            }
        }
        else if( b==c ){
            yi = a*z;
            y  = 1+yi;
            for( i=1; i<itermax; ++i ){
                yi *= (a+i)*z/(i+1);
                y  += yi;
                if( abs(yi) < tol )
                    break;
            }
        }
        else{
            if( abs(c) < 10*__eps__ ){
                f = std::numeric_limits<double>::infinity();
                return H2F1_BADARGS;
            }
            yi = a*b*z/c;
            y  = 1+yi;
            for( i=1; i<itermax; ++i ){
                yi *= (a+i)*(b+i)*z/((i+1)*(c+i));
                y  += yi;
                if( abs(yi) < tol )
                    break;
            }            
        }
        f = y;
        if(i==itermax)
            return H2F1_NONCONVERGENT;
        else
            return H2F1_SUCCESS;
    }
    
    /** The general function using the linear transformation formulas
        described by Michel and Stoitsov, valid for "almost all" z
        and "almost all" a, b, and c*/
    int hyperGeom2F1( const double& a, const double& b, const double& c, const double& z, double& f )
    {
        int result;
        double v1, v2;
        const double max_negative_threshold = -100000.0f;
        const double z_1_threshold = 1000*sqrt(__eps__);
        if( z == 1.0 ){
            /**
             * This is a very special case where z is exactly 1. In this
             * case we:
             *   - Use Michel & Stoitsov's eq. (10) with the transformation
             *     s = z/(z-1), which now tends to -infinity
             *   - Use a power series expansion obtained with Mathematica
             *   - Substitute s = z/(z-1) in the power series.
             *   - particularize for z=1
             * Bad news is that the power series expansion involves Gamma
             * functions, hence we might end up with bad combinations of
             * parameters of a, b, and c that will produce evaluations of
             * Gamma functions at negative integers
             */
            if( isnegint(a) || isnegint(b) || isnegint(c) ||
                isnegint(c-b-a) || isnegint(a+b-c) ||
                isnegint(c-a) || isnegint(c-b) ){
                f = NAN;
                return H2F1_GAMMAUNDEF;
            }
            // Everything allright, we can compute:
            double gc   = tgamma(c);
            double gcb  = tgamma(c-b);
            double gca  = tgamma(c-a);
            double gcba = tgamma(c-b-a);
            f = (gcba*gc)/(gcb*gca);
            return H2F1_SUCCESS;
        }
        else if(   ( z > 1.0-z_1_threshold ) && (z<1.0)   ){
            /**
             * In case z is a number smaller than, but very close
             * to, 1.0, Horchler's approach might not converge, hence
             * we will use a different strategy
             */
            // Try to use the regular series:
            result = hyperGeom2F1_azl1( a, b, c, z, f );
            if( result == H2F1_SUCCESS ) // It worked! nothing more to do:
                return H2F1_SUCCESS;
            else{ // Didn't work :-(
                /**
                 * The strategy to use here is:
                 *   - Use Michel & Stoitsov's eq. (10) with the transformation
                 *     s = z/(z-1), which now tends to -infinity
                 *   - Use a power series expansion obtained with Mathematica
                 * Bad news is that the power series expansion involves Gamma
                 * functions, hence we might end up with bad combinations of
                 * parameters of a, b, and c that will produce evaluations of
                 * Gamma functions at negative integers
                 */
                if( isnegint(a) || isnegint(b) || isnegint(c) ||
                    isnegint(c-b-a) || isnegint(a+b-c) ||
                    isnegint(c-a) || isnegint(c-b) ){
                    f = NAN;
                    return H2F1_GAMMAUNDEF;
                }
                else{
                    /**
                     * Note if we reach here means a+b-c+1.0 != 0  and c-b-a+1.0 != 0,
                     * since otherwise a+b-c and/or c-b-a would be negative integers
                     * and the test over the computations of Gamma function would have
                     * failed
                     */
                    /** NOTE:
                     * This block is problematic in case c-b-a < 0, since
                     * pow(var,c-b-a) may become arbitrarily large...
                     */
                    // Everything allright, we can compute:
                    double ga   = tgamma(a);
                    double gb   = tgamma(b);
                    double gc   = tgamma(c);
                    double gcb  = tgamma(c-b);
                    double gca  = tgamma(c-a);
                    double gcba = tgamma(c-b-a);
                    double gabc = tgamma(a+b-c);
                    double var  = (1.0-z)/z;
                    v1 = gcba/(gcb*gca) * (   1.0 -     a * (a-c+1.0) / (a+b-c+1.0) * var   ); // Can divide by (a+b-c+1.0), see above
                    v2 =   gabc/(ga*gb) * (   1.0 - (c-b) *   (1.0-b) / (c-b-a+1.0) * var   ); // Can divide by (c-b-a+1.0), see above
                    f  = gc * pow(z,-a) * (v1+v2*pow(var,c-b-a));
                    return H2F1_SUCCESS;
                }
            }
        }
        else if(   (z>0) && ( z <= 1.0-z_1_threshold )   ){
            /**
             * Here Horchler's approach should converge without issues
             */
            return hyperGeom2F1_azl1( a, b, c, z, f );
        }
        else if( z<=0.0 ){
            /**
             * If z is negative, then the module of z/(z-1) is always smaller
             * than 1, so that Horchler's approach should work fine after the
             * transformation defined in eq. (10) of Michel & Stoitsov's paper,
             * which is always numerically stable. We will recursively call
             * this function, and thence we already contemplate the case where
             * z -> -infty and z/(z-1) -> 1:
             */
            result = hyperGeom2F1( a, c-b, c, z/(z-1.0), f );
            if(result!=H2F1_SUCCESS){ return result; }
            f *= pow(1.0-z,-a);
            return H2F1_SUCCESS;
        }
        else if( z>1.0 ){
            /**
             * If we end up here, we have a positive number strictly greater
             * than 1. The only thing we can do is working with eq. (13) of
             * Michel & Stoitsov's paper, which involves the computation of
             * Gamma functions (and, therefore, numerical unstabilities...)
             */
            // Check first if we have evaluations of integer, negative
            // values for Gamma functions:
            if( isnegint(a) || isnegint(b) || isnegint(c) ||
                    isnegint(a-b) || isnegint(b-a) ||
                    isnegint(c-a) || isnegint(c-b) ){
                f = NAN;
                return H2F1_GAMMAUNDEF;
            }
            else{
                // Otherwise, directly implement eq. (13):
                result = hyperGeom2F1( a, 1.0-c+a, 1.0-b+a, 1.0/z, v1 );
                if( result!=H2F1_SUCCESS ){ f=v1; return result;}
                result = hyperGeom2F1( b, 1.0-c+b, 1.0-a+b, 1.0/z, v2 );
                if( result!=H2F1_SUCCESS ){ f=v2; return result;}
                double c1 = ( tgamma(c)*tgamma(b-a) )/( tgamma(b)*tgamma(c-a) );
                double c2 = ( tgamma(c)*tgamma(a-b) )/( tgamma(a)*tgamma(c-b) );
                if(z<0.0)
                    f = c1 * pow(-z,-a) * v1 + c2 * pow(-z,-b) * v2;
                else // Return the real part
                    f = c1 * pow(z,-a) * cos(a*PI) * v1 + c2 * pow(z,-b) * cos(b*PI) * v2;
                return H2F1_SUCCESS;
            }
        }
        else{ // Should never get here
            return H2F1_UNDEF;
        }
    }
    
    /**
     * A especialized function for the case we are interested in,
     * when a=1/2, b=g/2, c=3/2
     *
     * NOTE: The function becomes unstable if g>=2  and z->1.
     *       This shouldn't be an important issue within the
     *       Matlab toolbox: first, because x->1 is a degenerate
     *       case; second, because integer values of g (such as
     *       g=2) are managed in a different way, with recursive
     *       rules.
     */
    int hyperGeom2F1( const double& g, const double& z, double& f )
    {
        if( z > -100.0 ){
            /**
             * Simply call the general purpose function, which
             * might be problematic in the case z>1 because of
             * the use of Gamma functions. Note, howvever, that
             * the Matlab toolbox should never invoke this case:
             *   - because the argument z should never be
             *     greater than 1.
             *   - because the cases for "bad integer combinations"
             *     of a, b, and c are dealt with independently
             */
            return hyperGeom2F1( 0.5, g/2, 1.5, z, f );
        }
        else{ // z <= -100
            /**
             * This is a very large, negative number. We will use an
             * assymptotic expansion. Bad news is we need to distinguish
             * the cases g=1.0 and g~=1.0. Therefore, we will have to
             * check values of g close 1.0, not exactly 1.0, which we deal
             * with by linear interpolation
             */
            double mz = -z;
            double arg, ser;
            if( g==1.0 ){
                // Series expansion obtained with Mathematica. Since
                // abs(z)>=100, each new term gains at least one digit
                // of precission (because of the sqrt). The number of
                // terms guarantees precission up to __eps__
                arg  = sqrt(mz);
                ser  = log(4*mz) / 2.0 / arg;          arg *= mz;
                ser += 1.0/4.0/arg;                    arg *= mz;
                ser -= 3.0/32.0/arg;                   arg *= mz;
                ser += 5.0/96.0/arg;                   arg *= mz;
                ser -= 35.0/1024.0/arg;                arg *= mz;
                ser += 63.0/2560.0/arg;                arg *= mz;
                ser -= 77.0/4096.0/arg;                arg *= mz;
                ser += 429.0/28672.0/arg;              arg *= mz;
                ser -= 6435.0/524288.0/arg;            arg *= mz;
                ser += 12155.0/1179648.0/arg;          arg *= mz;
                ser -= 46189.0/5242880.0/arg;          arg *= mz;
                ser += 88179.0/11534336.0/arg;         arg *= mz;
                ser -= 676039.0/100663296.0/arg;       arg *= mz;
                ser += 1300075.0/218103808.0/arg;      arg *= mz;                
                ser -= 5014575.0/939524096.0/arg;      arg *= mz;
                ser += 646323.0/134217728.0/arg;       arg *= mz;
                ser -= 300540195.0/68719476736.0/arg;  arg *= mz;
                ser += 583401555.0/146028888064.0/arg;
                f = ser;
                return H2F1_SUCCESS;
            }
            else if( abs(g-1.0) < sqrt(__eps__) ){
                // Find by interpolation:
                double ref, val0, val1;
                int res;
                res = hyperGeom2F1( 1.0, z, val0 );
                if(res!=H2F1_SUCCESS){ f=val0; return res; }
                if(g>1.0)
                    ref = 1.0+sqrt(__eps__);
                else
                    ref = 1.0-sqrt(__eps__);
                res = hyperGeom2F1( ref, z, val1 );
                if(res!=H2F1_SUCCESS){ f=val1; return res; }
                // If both two evaluations were successfull, we can interpolate:
                double w1 = abs(g-1.0);
                double w0 = abs(g-ref);
                f = (w0*val0+w1*val1)/(w0+w1);
                return H2F1_SUCCESS;
            }
            else{
                const double tol = __eps__;         // DOUBLE_EPS = 2.2204e-16
                const unsigned int itermax = 32768; // 2^15
                unsigned int i;
                double corr, cumprod, fact;
                int sign = 1;
                // Series expansion obtained with Mathematica:
                arg = 1.0;
                ser = 1.0/(1.0-g);
                cumprod = g;
                fact    = 2.0;
                for( i=0; i<itermax; ++i ){
                    arg     *= mz;
                    corr     = sign * (cumprod/fact) * (g+2*i+1) / arg;
                    ser     += corr;
                    if(abs(corr)<tol)
                        break;
                    fact    *= (2*i+4);
                    cumprod *= (g+2*i+2);
                    sign     = -sign;
                }
                ser *= pow( mz, -g/2 );
                ser += sqrt(PI)/2.0 * tgamma((g-1)/2.0) / tgamma(g/2.0) / sqrt(mz);
                f = ser;
                if( i==itermax )
                    return H2F1_NONCONVERGENT;
                else
                    return H2F1_SUCCESS;
            }
        } // end z <= -100
    }
    
    bool isnegint(const double& z){
        if( z >= 10*__eps__ )
            return false;
        else if( abs(z) < 10*__eps__ )
            return true;
        else{
            unsigned int zi = (unsigned int)(-z);
            return ( abs((double)zi+z) < 10*__eps__ );
        }
    }
    
} // End namespace hypegeo

#endif // #ifndef _hypergeom2F1_cxx
