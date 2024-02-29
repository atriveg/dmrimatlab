function SH = atti2shadc( atti, gi, bi, varargin )
% function SH = atti2shadc( atti, gi, bi, 'opt1', value1, 'opt2', value2, ... )
%
%   Fits the normalized logarithm of a given attenuation signal (DWI signal 
%   over the baseline) in the basis of spherical harmonics by using
%   (Tikhonov regularized) least squares over a mono-exponential model.
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradient over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that multi-shell acquisitions are allowed.
%
%   The signal S_i is assumed to follow a mono-exponential model:
%
%      S_i(gi,bi) = S_0*exp(-bi*ADC(gi)),
%
%   so that:
%
%         [-log(S_i(gi,bi)) / S_0] / bi = ADC(gi)
%      -> dwi/bi = sum_{m=1}^{M} SH(m)*Y_m(gi)     [Linear Least Squares]
%
%   where ADC(gi) is the apparent diffusion coefficient, a function defined
%   over the unit sphere to be fit at each image voxel in the basis of
%   even spherical harmonics up to order L, so that the output of the 
%   function becomes:
%
%      SH: a MxNxPx((L+1)(L+2)/2) array with the coefficients computed at 
%         each voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      L: an even integer with the maximum order of the SH to be used
%         (default: 6).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem (default 0.006).
%      chunksz: the LLS problem reduces to the product of the dwi signal
%         by an inverse matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the LLS
%         inverse at each step, hence taking advantage of matlab's
%         capabilities (default: 1000).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-5, 1-1.0e-5).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% Check the mandatory input argments:
if(nargin<3)
    error('At least the atti volume, the gradient table, and the b-values must be supplied');
end
[M,N,P,G] = size(atti);
NV = M*N*P; % Total number of voxels to be processed
if(~ismatrix(gi)||~ismatrix(bi))
    error('gi and bi must be 2-d matlab matrixes');
end
if(size(gi,1)~=G)
    error('The number of rows in gi must match the 4-th dimension of dwi');
end
if(size(gi,2)~=3)
    error('The gradients table gi must have size Gx3');
end
if(size(bi,1)~=G)
    error('The number of b-values bi must match the 4-th dimension of dwi');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end

% Parse the optional input arguments:
opt.L = 6;              optchk.L = [true,true];       % always 1x1 double
opt.lambda = 0.006;     optchk.lambda = [true,true];  % always 1x1 double
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
opt.tl = 1.0e-5;        optchk.tl = [true,true];      % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];      % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});

% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------

% Compute the LS matix for SH fitting:
B   = GenerateSHMatrix( opt.L, gi );    % GxK, where K=(L+1)(L+2)/2
LR  = GenerateSHEigMatrix( opt.L );     % KxK
WLS = (B'*B+(opt.lambda).*LR.*LR)\(B');   % (KxK)^(-1) * (KxG) -> KxG
WLS = WLS'; % GxK, for convenience, see loop below

% Make sure the dwi lay within the proper range:
atti(atti>opt.tu) = opt.tu;
atti(atti<opt.tl) = opt.tl;
atti = -log(atti); % This is bi*ADC
% Normalize the channels to compute ADC instead of bi*ADC (for the
% mono-exponential model, the ADC does not depend on bi):
bi  = reshape(1./bi,[1,1,1,G]);
if(is_broadcast_available) % We can broadcast
    atti = atti.*bi;
else
    atti = bsxfun( @(x,y)(x.*y), atti, bi );
end

% Now, fit the data chunk-by-chunk via least squares where the mask is true
atti = reshape(atti,[NV,G]);  % NVxG
mask = opt.mask(:);           % NVx1
% Mask...
atti   = atti(mask,:); % PVxG
PV     = size(atti,1);
SHmask = zeros(PV,size(WLS,2)); % NVxK
for ck=1:ceil(PV/opt.chunksz)
    idi = (ck-1)*opt.chunksz+1;
    idf = min(ck*opt.chunksz,PV);
    SHmask(idi:idf,:) = atti(idi:idf,:)*WLS; % (chunksz x G) * (G x K) -> (chunksz x K)
end
% Cast the result to the proper size:
SH = zeros(NV,size(WLS,2)); % NVx(L+1)(L+2)/2
SH(mask,:) = SHmask;
SH = reshape(SH,[M,N,P,size(WLS,2)]);
