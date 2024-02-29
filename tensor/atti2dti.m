function [tensor,S0] = atti2dti( atti, gi, bi, varargin )
% function [tensor,S0n] = atti2dti( atti, gi, bi, 'opt1', value1, 'opt2', value2, ... )
%
%   Given the attenuation signal atti, computed as S_i/S_0, the gradients
%   table gi, and the b-values bi it corresponds to, computes the diffusion
%   tensor model and the correction to S_0 that best explains it by
%   following this procedure:
%     1- Use ordinary least squares (OLS) to fit the logarithm of atti to a
%        diffusion tensor model linearly.
%     2- If asked, refine the previous solution using weighted least
%        squares (WLS) with weight proportional to A_i^2, where A_i are
%        iteratily estimated using the previous value of the estimated
%        tensor.
%     3- If asked, use the previous solution as the initial iteration of an
%        iterative, non-linear estimation procedure in the natural (as
%        opposed to the logarithmic) domain. Here, it is the squared root
%        of the diffusion tensor, instead of the diffusion tensor itself,
%        which is actually estimated, ensuring the final solution is
%        positive semi-definite.
%   In either 1-) or 2-), the function allows to project the solution in
%   the space of positive semi-definite matrixes by either clipping
%   negative eigenvalues to 0 or computing their absolute value. 
%
%      atti: a MxNxPxG double array containing the signal sampled at G
%         directions at each voxel within the MxNxP image frame. In case
%         you want DT-MRI, signal should be log(atti).
%      gi: a Gx3 matrix with the directions sampled table, each row 
%         corresponding to a unit vector (signal(gi) is assumed to be equal
%         to signal(-gi));
%      bi: a Gx1 vector with the b-values corresponding to each gradient
%         direction. Note this value will be used only if the 'wls' flag is
%         marked true, so you can otherwise pass an empty array [].
%
%      tensor: Can be either:
%         MxNxPx6 if the 'unroll' option is switchwd off
%         MxNxPx3x3 if the 'unroll' option is switched on.
%         In the former case, no duplicates of tensor entries are returned,
%         so that squeeze(tensor(x,y,x,:)) = [D11,D12,D13,D22,D23,D33]'.
%      S0n: MxNxP, the normalized baseline (a value near 1), so that you
%         can recover your ESTIMATED baseline as S0.*S0n, assuming S0 was
%         the value you used to normalize atti (atti=dwi./S0).
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      wls: either (true) or not (false) using weighted least squares
%         instead of ordinary least squares.
%         (default: true).
%      nonlinear: either (true) or not (false) using the nonlinear fit
%         (default: false).
%      wlsit: in case wls is switched on, the maximum number of iterations
%         for the problem (default: 5).
%      wsc: the minimum weight to be applied in the WLS problem with
%         respect to the maximum one (default: 0.01).
%      tol: Levenberg-Marquardt's tolerance for the nonlinear optimization
%         procedure, meaning: iterations will stop if the relative change
%         in the solution (for a seccessful iteration) is below this 
%         threshold (default: 1.0e-6).
%      maxiters: maximum number of iterations in Levenberg-Marquardt's
%         algorithm (default: 100).
%      rcondth: minimum allowed reciprocal condition number for matrix
%         inversions (default: 1.0e-6).
%      fixmode: either 'n', 'z', or 'a', respectively meaning: no
%         correction of negative eigenvalues, zero-clipping or absolute
%         value (default: 'n').
%      unroll: wether (true) or not (false) output the result as a 3x3
%         matrix at each voxel instead of a 6x1 vector (in the latter case
%         duplicates of the entries of the diffusion tensor are removed so
%         that only [D11,D12,D13,D22,D23,D33] are returned) (default:
%         false).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      maxthreads: only for POSIX systems, the maximum number of threads
%         allowed (default: automatically determine).

%%% -----------------------------------------------------------------------
% Check the mandatory input argments:
if(nargin<3)
    error('At lest the atti, the gradients table, and the b-values must be supplied');
end
[M,N,P,G] = size(atti);
NV = M*N*P; % Total number of voxels to be processed
assert(ismatrix(gi),'gi must be 2-d matlab matrix');
assert(size(gi,1)==G,'The number of rows in gi must match the 4-th dimension of dwi');
assert(size(gi,2)==3,'The gradients table gi must have size Gx3');
%%% -----------------------------------------------------------------------
% Parse the optional input arguments:
opt.wls = true;         optchk.wls = [true,true];        % always 1x1 boolean
opt.nonlinear = false;  optchk.nonlinear = [true,true];  % always 1x1 boolean
opt.wlsit = 5;          optchk.wlsit = [true,true];      % always 1x1 double
opt.wsc = 0.01;         optchk.wsc = [true,true];        % always 1x1 double
opt.tol = 1.0e-6;       optchk.tol = [true,true];        % always 1x1 double
opt.maxiters = 100;     optchk.maxiters = [true,true];   % always 1x1 double
opt.rcondth = 1.0e-6;   optchk.rcondth = [true,true];    % always 1x1 double
opt.fixmode = 'n';      optchk.rcondth = [true,true];    % always 1x1 char
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true]; % always 1x1 char
opt.unroll = false;     optchk.unroll = [true,true];     % always 1x1 boolean
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%% -----------------------------------------------------------------------
% Now, check bi if necessary:
if(opt.wls)
    assert(size(bi,2)==1,'bi must be a column vector if wls is true');
    assert(size(bi,1)==G,'bi must be Gx1 if wls is true');
end
%%% -----------------------------------------------------------------------
% Now, fit the data chunk-by-chunk via least squares where the mask is true
atti  = reshape(atti,[NV,G]); % NVxG
mask  = opt.mask(:);            % NVx1
% Mask...
atti = atti(mask,:); % PVxG
%%% -----------------------------------------------------------------------
% Call the mex function
options.wlsit = opt.wlsit;
options.maxiters = opt.maxiters;
options.wsc = opt.wsc;
options.tol = opt.tol;
options.rcondth = opt.rcondth;
switch(opt.fixmode)
    case 'n'
    case 'z'
    case 'a'
    otherwise
        error(['Unknown value for ''fixmode'' option: ',opt.fixmode]);
end
options.fixmode = opt.fixmode;
if(opt.nonlinear)
    options.mode = 'p';
elseif(opt.wls)
    options.mode = 'w';
else
    options.mode = 'o';
end
[tmask,smask] = atti2dti_( double(atti'), double(gi), double(bi), options, opt.maxthreads );
%%% -----------------------------------------------------------------------
% Cast the result to the proper size:
S0             = zeros(NV,1); % NVx1
tensor         = zeros(NV,6); % NVx6
tensor(mask,:) = tmask';
tensor         = reshape(tensor,[M,N,P,6]);
if(opt.unroll)
    tensor = tensor(:,:,:,[1,2,3,2,4,5,3,5,6]);
    tensor = reshape(tensor,[M,N,P,3,3]);
end
S0(mask)       = smask;
S0             = reshape(S0,[M,N,P]);
%%% -----------------------------------------------------------------------
