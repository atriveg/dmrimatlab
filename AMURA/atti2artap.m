function rtap = atti2artap( atti, gi, bi, varargin )
% function rtap = atti2artap( atti, gi, bi, 'option1', value1, ... )
%
%   Computes the (apparent) Return To Axis Probability according to a 
%   mono-exponential model for single-shell acquisitions as described in
%   AMURA:
%
%      atti(u,b) = exp(-b*D_0(u)),
%
%   where D_0(u) is the Apparent Diffusion Coefficient measured at b=b0 for
%   each direction 'u' within the unit sphere.
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
%      rtap: MxNxP double array with the computed RTAP.
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
%      tensor: wether (true) or not (false) using a tensor model instead of
%         a model-free ODF estimation to compute the maximum diffusion
%         direction. Using a tensor estimation is usually more robust
%         (default: true).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem of fitting SH coefficients (note the order
%         L used for SH is internally computed) (default 0.001).
%
%   Sanity checks on the attenuation signal:
%
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%
%   Advanced parameters:
%
%      tau: 1x1, the effective diffusion time of the dMRI sequence in
%         miliseconds (default: 70.0e-3).
%      chunksz: the evaluation of SH at desired directions is done by 
%         repeatedly calling GenerateSHMatrix. This is done chunk-by-chunk
%         for efficiency (default: 256).
%      clean: 1x1 double in the range [0,100]. This is a simple outlier 
%         rejection parameter to avoid out-of-range values: 0 means no
%         outlier rejection is applied; >0 means outlier rejection is
%         applied, the closer to 100 the more agressive (default: 50).

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
opt.tensor = true;      optchk.tensor = [true,true];   % always 1x1 boolean
opt.lambda = 0.001;     optchk.lambda = [true,true];   % always 1x1 double
% -------------------------------------------------------------------------
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
% -------------------------------------------------------------------------
opt.tau = 70.0e-3;      optchk.tau = [true,true];      % always 1x1 double
opt.chunksz = 256;      optchk.chunksz = [true,true];  % always 1x1 double
opt.clean = 50;         optchk.clean = [true,true];    % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Take the maximum available SH order:
L = (sqrt(8*(G-1)+9)-3)/2;
L = 2*floor(L/2);
% -------------------------------------------------------------------------
% Compute the maximum diffusivity direction
if(opt.tensor) % Use tensor model
    sh = atti2shadc( atti, gi, bi, 'L', 2, 'lambda', 0.001, ...
        'tl', 1.0e-5, 'tu', 1-1.0e-5, 'mask', opt.mask, ...
        'chunksz', opt.chunksz );                  % MxNxPx6
    tens   = shadc2dti(sh,'mask',opt.mask,...
        'chunksz',opt.chunksz,'unroll',false );    % MxNxPx6
    u0 = dti2xyz( tens, 'mask', opt.mask );        % MxNxPx3
else           % Use model-free method
    sh = atti2shadc( atti, gi, bi, 'L', min(L,8), 'lambda', 0.006, ...
        'tl', 1.0e-5, 'tu', 1-1.0e-5, 'mask', opt.mask, ...
        'chunksz', opt.chunksz );                  % MxNxPxK
    Di = sh2signal( sh, gi, ...
        'chunksz', opt.chunksz, 'mask', opt.mask); % MxNxPxG
    [~,id] = max(Di,[],4);
    % id has size M x N x P; its entries are the arg. max. in the
    % gradients table:
    gx = gi(:,1); % G x 1
    gy = gi(:,2); % G x 1
    gz = gi(:,3); % G x 1
    u0 = zeros(M,N,P,3);  % M x N x P x 3
    u0(:,:,:,1) = gx(id); % M x N x P
    u0(:,:,:,2) = gy(id); % M x N x P
    u0(:,:,:,3) = gz(id); % M x N x P
end
% -------------------------------------------------------------------------
% Compute the measurement using the generic fucntion:
rtap = atti2amura( atti, gi, bi, ...
    'type', 'rtap', 'u0', u0, 'mask', opt.mask, ...
    'L', L, 'lambda', opt.lambda, ...
    'tl', opt.tl, 'tu', opt.tu, ...
    'tau', opt.tau, 'chunksz', opt.chunksz, 'clean', opt.clean );
