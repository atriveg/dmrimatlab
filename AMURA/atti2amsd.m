function msd = atti2amsd( atti, gi, bi, varargin )
% function msd = atti2amsd( atti, gi, bi, 'option1', value1, ... )
%
%   Computes the (apparent) Non-Gaussianity NG) according to a 
%   mono-exponential model for single-shell acquisitions as described in
%   AMURA:
%
%      atti(u,b) = exp(-b*D_0(u)),
%
%   where D_0(u) is the Apparent Diffusion Coefficient measured at b=b0 for
%   each direction 'u' within the unit sphere. The NG is defined in terms
%   of the 'angle' between the true attenuation signal and its Gaussian
%   counterpart, which in AMURA is modeled as that with a 'cropped' version
%   of the Diffusion Coefficent (with SH coefficients up to order L=2):
%   D_Gauss(u) = sum for l<=2 c_l^m Y_l^m(u)
%
%   Inputs:
%
%      atti: a MxNxPxG double array with the attenuation signal measured at
%         each voxel within the MxNxP field of view.
%      gi: a Gx3 double array with the gradients table, i.e. each row of gi
%         should be a unit vector describing the corresponding gradient
%         direction.
%      bi: a Gx1 double array with the corresponding b-values of each entry
%         in gi. NOTE: for the AMURA model to make sense, all bi should be
%         similar, otherwise the sampling scheme is a multi-shell and a
%         differente model (such as MiSFIT) should be used. Alternatively,
%         bi can be a single scalar describing the acquired shell.
%
%   Outputs:
%
%      msd: MxNxP double array with the computed NG.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).
%
%   Parameters related to SH computations:
%
%      L: the maximim order for the SH expansions used to represent
%         spherical functions, wich must be greater than 2 (default: 8).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem of fitting SH coefficients (note the order
%         L used for SH is internally computed) (default: 0.001).
%      tau: the effective diffusion time in seconds (default: 70.0-3).
%
%   Sanity checks on the attenuation signal:
%
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%
%   Other optional parameters:
%
%      chunksz: SH expansions are computed as matrix products, which are 
%         performed chunk-by-chunk to avoid memory issues. This parameter
%         sets the size of such chunks (default: 256).

% Check the mandatory input arguments:
if(nargin<3)
    error('At least atti, gi, and bi must be supplied');
end
[M,N,P,G] = size(atti);
assert(ismatrix(gi),'gi must be a Gx3 matrix');
assert(size(gi,2)==3,'The second dimension of gi must have size 3');
assert(size(gi,1)==G,'The first dimension of gi must match the fourth dimension of atti');
if(isscalar(bi))
    bi = repmat(bi,[G,1]);
else
    assert(ismatrix(bi),'bi must be either Gx1 or 1x1');
    assert(size(bi,2)==1,'bi should be a column vector');
    assert(size(bi,1)==G,'The first dimension of bi should be either singleton or match the fourth dimension of atti');
end
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
% -------------------------------------------------------------------------
opt.L = 8;              optchk.L =[true,true];         % always 1x1 double
opt.lambda = 0.001;     optchk.lambda = [true,true];   % always 1x1 double
opt.tau = 70.0e-3;      optchk.tau = [true,true];      % always 1x1 double
% -------------------------------------------------------------------------
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
% -------------------------------------------------------------------------
opt.chunksz = 256;      optchk.chunksz = [true,true];  % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
assert(abs(opt.L-floor(opt.L))<100*eps,'Optional argument ''L'' must be integer');
assert(abs(opt.L/2-floor(opt.L/2))<100*eps,'Optional argument ''L'' must be even');
% -------------------------------------------------------------------------
sh = atti2shadc( atti, gi, bi, ...
    'L', opt.L, 'lambda', opt.lambda, 'tl', opt.tl, 'tu', opt.tu, ...
    'mask', opt.mask, 'chunksz', opt.chunksz );
msd = 6*(opt.tau)*sh(:,:,:,1)/sqrt(4*pi);
end
