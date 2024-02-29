function I2 = clean_moments_outliers(I,strength)
% function I2 = clean_moments_outliers(I,strength)
%
%   This function is used to avoid out-of-range values in computational
%   dMRI images assumed to eval in the range [0,infinity), mainly arbitrary
%   order moments: from a minimum percentile to a maximum percentile, the
%   algorithm checks if the jump in the values is excessively large, so
%   that all the correspondign values above such percentile are considered
%   outliers. The agressiveness of the procedure is tuned with the strength
%   parameter, which ranges [0,100]
%
%     I:        N-dimensional array with arbitrary dimensionality
%     strength: 1x1 double in [0,100]
%     I2:       "clean" image the same size as I

if(strength<100*eps)
    I2 = I;
    return;
end

I(I<0)      = 0;
I(isinf(I)) = 0;
I(isnan(I)) = 0;

NS       = 20;
RATIO    = interp1( (0:10:100), log10([1000,100,10,10-8.5*(1:8)/8]), strength );
RATIO    = 10^RATIO;
PMIN     = interp1( (0:10:100), [99.9,99.5,99,99-9*(1:8)/8], strength );
PMAX     = interp1( (0:10:100), log10([1.0e-4,1.0e-3,0.01,0.01+0.09*(1:8)/8]), strength );
PMAX     = 100 - 10^PMAX;
step_per = (PMAX-PMIN)*(1:NS)/NS + PMIN;
Val0     = prctile( I(:), PMIN );

for ii=1:NS
    Val1   = prctile( I(:), step_per(ii) );
    RatioV = Val1/Val0;
    if(RatioV<=RATIO)
        Val0 = Val1;
    else
        break;
    end
end

I2 = I;
I2(I>Val0) = Val0;
