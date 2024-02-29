function [S0,nit] = dmri_inverse_Laguerre_fcn(S1,sigma)
% function A = dmri_inverse_Laguerre_fcn(X[,sigma])
%
%   This is a quite dumb correction for Rician noise corrupted signals.
%   Since the expectation of a Rician distributed variable for a noise-free
%   value A and a noise power sigma in the complex domain is:
%
%     X = sigma*sqrt(pi/2)*L_{1/2}(-A^2/2/sigma^2),
%
%   we aim at inverting the Laguerre polynomial to trace back the
%   noise-free signal A from the noisy measurement X.
%
%   You can call this function as either:
%
%      [1] A = dmri_inverse_Laguerre_fcn(X), and it is assumed that you
%          have properly normalized your measurements dividing them by the
%          noise power sigma, so that the recovered noise-free signal will
%          be returned as well.
%      [2] A = dmri_inverse_Laguerre_fcn(X,sigma), with a scalar value of
%          sigma, so that both the measurements and the output signal will
%          be unnormalized.
%      [3] A = dmri_inverse_Laguerre_fcn(x,sigma), with sigma the same size
%          as A, so that a different sigma will be applied to each
%          component of A.
%
%   This way, these two use cases aere exactly the same, wether sigma is a
%   scalar or an array the same size as A:
%
%      >> X = X./sigma;
%      >> A = dmri_inverse_Laguerre_fcn(X);
%      >> A = A.*sigma;
%
%      >> A = dmri_inverse_Laguerre_fcn(X,sigma);
%
%   NOTE: Input S1 is assumed to be always positive, so that the function
%         will internally set S1(S1<0) = 0;

S1       = double(S1);
S1(S1<0) = 0;
switch(nargin)
    case 1
        sz = size(S1);
        [S0,nit] = dmri_inverse_Laguerre_fcn_core(S1(:));
        S0  = reshape(S0,sz);
        nit = reshape(nit,sz);
    case 2
        sz = size(S1);
        assert( isscalar(sigma) || isequal(size(sigma),sz), ...
            'sigma must be either an scalar or match the size of S1');
        S1 = S1./sigma;
        [S0,nit] = dmri_inverse_Laguerre_fcn_core(S1(:));
        S0  = reshape(S0,sz).*sigma;
        nit = reshape(nit,sz);
    otherwise
        error('Only 1 or 2 input arguments are allowed');
end

function [S0,nit] = dmri_inverse_Laguerre_fcn_core(S1)
% ----------
p1      = (S1<=sqrt(pi/2));
S1(p1)  = 0;
p2      = (S1>50); % Beyond this value, matlab's built-in Laguerre returns NaN
S1(p2)  = ( S1(p2) + sqrt(S1(p2).*S1(p2)-2) )/2;
% ----------
p1      = find(~(p1|p2));
N       = length(p1);
chunksz = 10000;
S0      = S1;
nit     = zeros(size(S0,1),1);
for ck=1:ceil(N/chunksz)
    % -----------
    idi       = (ck-1)*chunksz+1;
    idf       = min(ck*chunksz,N);
    % -----------
    [S0_,nit_] = dmri_inverse_Laguerre_fcn_core_chunk(S1(p1(idi:idf),1));
    S0(p1(idi:idf),1)  = S0_;
    nit(p1(idi:idf),1) = nit_;
    % -----------
end

function [S0,nit] = dmri_inverse_Laguerre_fcn_core_chunk(S1)
% Inverts: y = sqrt(pi/2)*Laguerre(1/2,-x^2/2)
% -------------------------------------------------------------------------
% S1 is in the range (sqrt(pi/2),50], since this function is called from
% dmri_inverse_Laguerre_fcn_core, which ensures that condition.
S0  = S1;
% The derivative of the Laguerre function becomes null at zero, so we
% cannot use Newton-Raphson's iterations everywhere:
dth = 0.1;
S0(S0<laguerreL_1_2(dth)) = 0; % To begin iterations at the origin, where the Taylor expansion is very precise
% Initialize iterations:
S0n = S0;
nit = zeros(size(S0));
p2  = true(size(S1));
% -------------------------------------------------------------------------
for n=1:10000
    if(~any(p2)), break; end
    % ---------------------------------------------------------------------
    nit(p2) = nit(p2)+1;
    p21 = p2 & (S0<laguerreL_1_2(dth)); % Polynomial approximation of higher order
    p22 = p2 & (~p21);   % Newton's
    % ---------------------------------------------------------------------
    % Quadratic approximation:
    if(any(p21))
        [y,y1,y2] = laguerreL_1_2(S0(p21));
        S0n(p21)  = ( -y1 + sqrt(abs(y1.*y1-2*y2.*(y-S1(p21)))) )./y2 + S0(p21);
    end
    % ---------------------------------------------------------------------
    % Newton's:
    if(any(p22))
        [y,y1]   = laguerreL_1_2(S0(p22));
        S0n(p22) = S0(p22) - (y-S1(p22))./y1;
    end
    % ---------------------------------------------------------------------
    p2 = p2 & ( abs(laguerreL_1_2(S0n)-S1)./abs(S1) > 1.0e-12 );
    S0 = S0n;
    % ---------------------------------------------------------------------
end
% -------------------------------------------------------------------------
