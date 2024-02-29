function [atti,gi,bi] = dwi2unbiasedatti(dwi,gi,bi,sigma)
% function [atti,gi,bi] = dwi2unbiasedatti(dwi,gi,bi,sigma)
%
%   INPUTS:
%
%    Computes the attenuation signal from the dwi signal and corrects the
%    former (assuming Rician noise) to avoid the bias.
%
%       dwi: MxNxPxG double array, the original signal acquired at the
%          scanner.
%       gi: Gx3, the gradients table (including baselines).
%       bi: Gx1, the b-values (including baselines).
%       sigma: either a scalar or a MxNxP double array with the estandard
%          deviation of noise in the original complex domain.
%
%   OUTPUTS:
%
%       atti: MxNxPxG2 double array, the attenuation signal (gradients over
%          baselines) with bias correction.
%       gi: G2x3, the gradients table without baselines.
%       bi: G2x1, the b-values without baselines.

[M,N,P,G] = size(dwi);
assert(ismatrix(gi),'gi must be 2-dimensional');
assert(ismatrix(bi),'bi must be 2-dimensional');
assert(size(gi,1)==G,'The number of entries (rows) in gi must match the 4-th dimension of dwi');
assert(size(gi,2)==3,'gi must be Gx3');
assert(size(bi,1)==G,'The number of entries in bi must match the 4-th dimension of dwi');
assert(size(bi,2)==1,'bi must be Gx1');
if(isscalar(sigma))
    sigma = ones(M,N,P)*sigma;
else
    assert( isequal(size(sigma),[M,N,P]), 'sigma must be either scalar or math the first 3 dimensions of atti' );
end

pbsl = (bi<1);
pgrd = (bi>=1);

bsl  = mean(dwi(:,:,:,pbsl),4);
atti = dwi(:,:,:,pgrd);

bsl  = correct_Rician_bias(bsl,sigma);
atti = correct_Rician_bias(atti,sigma);

if(is_broadcast_available_test)
    atti = atti./bsl;
else
    atti = bsxfun( @(x,y)(x./y), atti, bsl );
end

bi = bi(pgrd,:);
gi = gi(pgrd,:);

function S0 = correct_Rician_bias(S1,sigma)
% The mean of a Rician distribution is:
%
%   E{X} = sigma*sqrt(pi/2)*Laguerre(1/2,-nu^2/2/sigma^2)
is_broadcast_available = is_broadcast_available_test;
if(is_broadcast_available)
    S0 = S1./sigma;
else
    S0 = bsxfun( @(x,y)(x./y), S1, sigma );
end
S0 = dmri_inverse_Laguerre_fcn(S0);
if(is_broadcast_available)
    S0 = S0.*sigma;
else
    S0 = bsxfun( @(x,y)(x.*y), S0, sigma );
end 
