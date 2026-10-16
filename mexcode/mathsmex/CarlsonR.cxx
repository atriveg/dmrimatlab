#ifndef _CarlsonR_cxx_
#define _CarlsonR_cxx_

#include "CarlsonR.h"
#include "math.h"
#include <limits>

namespace carlson
{

    int CarlsonRd( const ElementType& x0, const ElementType& y0, const ElementType& z0, ElementType& result )
    {
        // Avoid modifying x0, y0, z0 in-place
        ElementType x = x0;
        ElementType y = y0;
        ElementType z = z0;

        // Bad inputs:
        if( isnan(x) || isnan(y) || isnan(z) ){
            result = nan("");
            return CarlsonR_BADARGS;
        }
        if( (x<0) || (y<0) || (z<=0) ){
            result = nan("");
            return CarlsonR_BADARGS;
        }
        if( x+y == 0 ){
            result = nan("");
            return CarlsonR_BADARGS;
        }
        // Special cases:
        if( isinf(x) || isinf(y) || isinf(z) ){
            result = 0.0; // ??
            return CarlsonR_SUCCESS;
        }
        if( (x==0) && (y==0) ){
            result = INFINITY;
            return CarlsonR_SUCCESS;
        }
        if( y==z ){
            if( x==0 ){
                result = 3*PI/(4*y*sqrt(y));
                return CarlsonR_SUCCESS;
            }
            if( x==y ){
                result = 1/(x*sqrt(x));
                return CarlsonR_SUCCESS;
            }
        }
        // General case:
        int success = CarlsonR_NONCONVERGENT;
        const unsigned int maxiters = 10000000;
        ElementType sm = 0;
        ElementType fc = 1;
        ElementType rx, ry, rz, lm, m, dx, dy, dz;
        for( unsigned int it=0; it<maxiters; ++it ){
            rx = sqrt(x);
            ry = sqrt(y);
            rz = sqrt(z);
            lm = rx*(ry + rz) + ry*rz;
            sm = sm + fc/(rz*(z + lm));
            fc = fc/4;
            x  = (lm + x)/4;
            y  = (lm + y)/4;
            z  = (lm + z)/4;
            m  = (x + y + 3*z)/5;
            if(isinf(m)){
                result = nan("");
                return CarlsonR_NONCONVERGENT;
            }
            dx = (m - x)/m;
            dy = (m - y)/m;
            dz = (m - z)/m;
            if( (abs(dx)<Carlson_ERRTOLd) && (abs(dy)<Carlson_ERRTOLd) && (abs(dz)<Carlson_ERRTOLd) ){
                success = CarlsonR_SUCCESS;
                break;
            }
        }
        if( success==CarlsonR_SUCCESS ){
            ElementType ea = dx*dy;
            ElementType eb = dz*dz;
            ElementType ec = ea - eb;
            ElementType ed = ea - 6*eb;
            ElementType ee = ed + ec + ec;

            ElementType s  = ed*(C1d + C5d*ed - C6d*dz*ee);
            s   = s + dz*(C2d*ee + dz*(dz*C4d*ea - C3d*ec));

            result  = 3*sm + fc*(1 + s)/(m*sqrt(m));

        }
        return success;
    }

    int CarlsonRf( const ElementType& x0, const ElementType& y0, const ElementType& z0, ElementType& result )
    {
        // Avoid modifying x0, y0, z0 in-place
        ElementType x = x0;
        ElementType y = y0;
        ElementType z = z0;

        // Bad inputs:
        if( isnan(x) || isnan(y) || isnan(z) ){
            result = nan("");
            return CarlsonR_BADARGS;
        }
        if( (x<0) || (y<0) || (z<0) ){
            result = nan("");
            return CarlsonR_BADARGS;
        }
        if( (x+y==0) || (x+z==0) || (y+z==0) ){
            result = nan("");
            return CarlsonR_BADARGS;
        }
        // Special cases:
        if( isinf(x) || isinf(y) || isinf(z) ){
            result = 0.0;
            return CarlsonR_SUCCESS;
        }
        ElementType t;
        if(z==0){
            if(x > 0){
                // swap x and z
                t = x; x = z; z = t;
            }
        }
        if(y==0){
            if(x>0){
                // swap x and y
                t = x; x = y; y = t;
            }
        }
        if( (x==0) && (y==0) ){
            result = INFINITY;
            return CarlsonR_SUCCESS;
        }
        if( (x==y) && (y==z) ){
            result = 1/sqrt(x);
            return CarlsonR_SUCCESS;
        }
        if( y==z ){
            if( x==0 ){
                result = PI/(2*sqrt(y));
                return CarlsonR_SUCCESS;
            }
            else{
                return CarlsonRc(x,y,result);
            }
        }
        // General case:
        int success = CarlsonR_NONCONVERGENT;
        const unsigned int maxiters = 10000000;
        ElementType rx, ry, rz, lm, av, dx, dy, dz;
        for( unsigned int it=0; it<maxiters; ++it ){
            rx = sqrt(x);
            ry = sqrt(y);
            rz = sqrt(z);
            lm = rx*(ry + rz) + ry*rz;
            x  = (lm + x)/4;
            y  = (lm + y)/4;
            z  = (lm + z)/4;
            av = (x + y + z)/3;
            if(isinf(av)){
                result = nan("");
                return CarlsonR_NONCONVERGENT;
            }
            dx = (av - x)/av;
            dy = (av - y)/av;
            dz = (av - z)/av;
            if( (abs(dx)<Carlson_ERRTOLf) && (abs(dy)<Carlson_ERRTOLf) && (abs(dz)<Carlson_ERRTOLf) ){
                success = CarlsonR_SUCCESS;
                break;
            }
        }
        if( success==CarlsonR_SUCCESS ){
            ElementType e2 = dx*dy - dz*dz;
            ElementType e3 = dx*dy*dz;
            ElementType  s = e2*(C1f*e2 - C2f - C3f*e3) + C4f*e3;
            result = (1.0 + s)/sqrt(av);
        }
        return success;
    }

    int CarlsonRc( const ElementType& x0, const ElementType& y0, ElementType& result )
    {
        // Avoid modifying x0, y0, z0 in-place
        ElementType x = x0;
        ElementType y = y0;

        // Bad inputs:
        if( isnan(x) || isnan(y) ){
            result = nan("");
            return CarlsonR_BADARGS;
        }
        if( (x<0) || (y==0) ){
            result = nan("");
            return CarlsonR_BADARGS;
        }

        // Special cases:
        if( isinf(x) || isinf(y) ){
            result = 0.0;
            return CarlsonR_SUCCESS;
        }
        if(y==0){
            // NOTE: we have already checked if y==0 and detected it
            // as a bad input condition. However, this is done also
            // in the original code:
            result = INFINITY;
            return CarlsonR_SUCCESS;
        }
        if( (x==0) && (y>0) ){
            result = PI/(2*sqrt(y));
            return CarlsonR_SUCCESS;
        }
        if( x==y ){
            result = 1.0/sqrt(x);
            return CarlsonR_SUCCESS;
        }

        // General case:
        ElementType w;
        if(y>0)
            w = 1;
        else{
            w = sqrt(x/(x - y));
            x = x - y;
            y = -y;
        }
        int success = CarlsonR_NONCONVERGENT;
        const unsigned int maxiters = 10000000;
        ElementType p,s,m;
        for( unsigned int it=0; it<maxiters; ++it ){
            p = 2*sqrt(x*y) + y;
            x = (x + p)/4;
            y = (y + p)/4;
            s = x + y + y;
            m = s/3;
            s = (y - x)/s;
            if(isnan(s)){
                result = nan("");
                return CarlsonR_NONCONVERGENT;
            }
            if(abs(s)<Carlson_ERRTOLc){
                success = CarlsonR_SUCCESS;
                break;
            }
        }
        if( success==CarlsonR_SUCCESS ){
            s  = (C1c + s*(C2c + s*(C3c + s*C4c)))*s*s;
            result = w*((1 + s)/sqrt(m));
        }
        return success;
    }
    
} // End namespace carlson

#endif // #ifndef _CarlsonR_cxx_
