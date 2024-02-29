function v = rice_nu2var(nu0)
% function v = rice_nu2var(nu0)
%
%    Computes the normalized variance of a Rician variable as a function of
%    its normalized central parameter. I.e:
%
%       Let nu be a real number and sigma a real positive number. Let nc
%       and ns be two Gaussian random variables with zero mean and standard
%       deviation 1, so that the Rician variable X is computed as:
%
%          X = sqrt( (nu+sigma*nc)^2+ns^2 )
%
%    Then the normalized variance of X, Var{X}/sigma^2, will depend just on
%    nu0=nu/sigma. This function computes this functional dependence. Note
%    that, for nu0->infinity, the Rician variable becomes closely Gaussian
%    and Var{X}/sigma^2 assymptotically approaches 1. On the opposite, for
%    nu0->0 the Rician variable becomes closely Rayleigh and Var{X}/sigma^2
%    assymptotically approaches (4-pi)/2. The transition occurs at nu/sigma
%    = 1.4191, for wich Var{X}/sigma^2 = ( (4-pi)/2 + 1 )/2.

nu0 = double(nu0);
% Compute the mean of the Rician variable:
mu1 = laguerreL_1_2(nu0);
% Compute the mean squared value:
mu2 = nu0.*nu0 + 2;
% Compute the variance:
v   = mu2 - mu1.*mu1;
