function sh = sh2squaredsh( shsqrt, varargin )
% function sh = sh2squaredsh( shsqrt, 
%                                 'option1', value1, 'option2', value2... )
%
%   Computes the SH coefficients of a non-negative ODF, shodf, from the SH
%   coefficients of its squared root, shsqrt, using the analytical
%   relationships between them provided by Wigner's symbols.
%
%   NOTE: if the squared root of the ODF spans up to a maximum order L
%   within the SH basis, the ODF itself spans up to an order 2L.
%
%   NOTE: if the ODF needs to be unit-mass, then the SH coefficients of its
%   squared root must be a unit-norm vector.
%
% INPUTS:
%
%   shsqrt: N1 x N2 x N3 x ... x ND x K, with K=(L+1)(L+2)/2 for some
%      non-negative, even integer L. The SH coefficients of the squared
%      root of the ODF
%
% OUTPUTS:
%
%   shodf: N1 x N2 x N3 x ... x ND x Kp, with Kp=(2L+1)2(L+2)/2. The SH
%      coefficients of the ODF.
%
% OPTIONAL ARGUMENTS PASSED AS OPTION/VALUE PAIRS:
%
%   mask: N1 x N2 x N3 x ... x ND, a mask to avoid unnecesary computations.
%      Zero valued mask positions are filled with zeros in the output
%      array (default: all ones)
%   maxthreads: the algorithm is run with multiple threads. This is the 
%      maximum allowed number of threads, which can indeed be reduced if
%      it exceeds the number of logical cores (default: the number of 
%      logical cores in the machine).
%

% Sanity checks on the size of shsqrt:
sz  = size(shsqrt);
K   = sz(end);
fov = sz(1:end-1);
if(isscalar(fov))
    fovm = [fov,1];
else
    fovm = fov;
end
L   = (sqrt(8*(K-1)+9)-3)/2; % SH order
if(abs(L-round(L))>1.0e-9)
    error('Weird size of the SH volume. Its fourth dimension should have size 1, 6, 15, 28, 45, ..., (L+1)(L+2)/2, with even L');
end
Kp  = (2*L+1)*(2*L+2)/2;

% Parse the optional input arguments:
opt.mask = true(fovm);  optchk.mask = [true,true];
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true];
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Reshape the inputs to work comfortably:
shsqrt = reshape( shsqrt,   [prod(fov),K] );
mask   = reshape( opt.mask, [prod(fov),1] );

% Mask the input SH volume:
shsqrt = shsqrt( mask, : );

% Call the mex implementation:
sh_ = sh2squaredsh_(shsqrt,opt.maxthreads);

% Create a buffer to store the output:
sh  = zeros( prod(fov), Kp );
sh(mask,:) = sh_;

% Reshape the output:
sh = reshape( sh, [fov,Kp] );
end
