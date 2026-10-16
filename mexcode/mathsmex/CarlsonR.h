/** 
 * This is an implementation of Carlson's elliptic integrals R_F and R_D,
 * used to compute mean kurtosis values. It has been directly coded from
 * Milan Batista's elfun18 toolbox:
 *
 * https://es.mathworks.com/matlabcentral/fileexchange/65915-elfun18/files/elFun18v1_5/elfun18/rd.m
 *
 * and the theoretical mathematical foundations of the method can be
 * checked in:
 *
 *  [1] B.C.Carlson,"Numerical computation of real or complex elliptic
 *      integrals", Numerical Algorithms 10, 1995, pp 13-26.
 *  [2] B.C.Carlson, "Elliptic Integrals" In: F.W.J.Olver et al., NIST
 *      Handbook of Mathematical Functions, Chap.19, Cambridge, 2010.
 *      Handbook of Mathematical Functions, Cambridge, 2010.
 *  [3] W.H.Press et al., "Numerical Recipes in Fortran 90", 2nd ed.,
 *      Cambridge, 1996
 */

#ifndef _CarlsonR_h_
#define _CarlsonR_h_

#include <cmath>
#include "mexToMathsTypes.h"

namespace carlson
{

/**
 * Return values on different error conditions:
 */
#define CarlsonR_SUCCESS 0
#define CarlsonR_BADARGS -1
#define CarlsonR_NONCONVERGENT -2
    
#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef Carlson_ERRTOLd
#define Carlson_ERRTOLd 0.002
#endif

#ifndef Carlson_ERRTOLf
#define Carlson_ERRTOLf 0.001
#endif

#ifndef Carlson_ERRTOLc
#define Carlson_ERRTOLc 3.0e-4
#endif

const ElementType C1d = -3.0/14.0;
const ElementType C2d = 1.0/6.0;
const ElementType C3d = 9.0/22.0;
const ElementType C4d = 3.0/26.0;
const ElementType C5d = 9.0/88.0;
const ElementType C6d = 9.0/52.0;

const ElementType C1f = 1.0/24.0;
const ElementType C2f = 1.0/10.0;
const ElementType C3f = 3.0/44.0;
const ElementType C4f = 1.0/14.0;

const ElementType C1c = 3.0/10.0;
const ElementType C2c = 1.0/7.0;
const ElementType C3c = 3.0/8.0;
const ElementType C4c = 9.0/22.0;

int CarlsonRd( const ElementType& x, const ElementType& y, const ElementType& z, ElementType& result );

int CarlsonRf( const ElementType& x, const ElementType& y, const ElementType& z, ElementType& result );

int CarlsonRc( const ElementType& x, const ElementType& y, ElementType& result );

    
} // End namespace carlson

#endif //_CarlsonR_h_
