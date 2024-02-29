function tensor = signal2dti( signal, gi, bi, varargin )
% function tensor = signal2dti( signal, gi, bi, 'opt1', value1, 'opt2', value2, ... )
%
%   Takes a generic symmetric signal defined over the unit sphere and 
%   represents it as a simple tensor-driven attenuation by means of Least
%   Squares.
%
%      signal: a MxNxPxG double array containing the signal sampled at G
%         directions at each voxel within the MxNxP image frame. In case
%         you want DT-MRI, signal should be -log(atti)/bi.
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
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      wls: either (true) or not (false) using weighted least squares
%         instead of simple least squares. Notice the former can be
%         dramatically slower than the latter, especially under Matlab
%         (default: false).
%      wlsit: in case wls is switched on, the maximum number of iterations
%         for the problem (default: 1).
%      wsc: the minimum weight to be applied in the WLS problem with
%         respect to the maximum one (default: 1.0e-6).
%      ADC0: free water diffusivity at body temperature, used to supress
%         badly behaved voxels with WLS (default: 3.0e-3).
%      unroll: wether (true) or not (false) output the result as a 3x3
%         matrix at each voxel instead of a 6x1 vector (in the latter case
%         duplicates of the entries of the diffusion tensor are removed so
%         that only [D11,D12,D13,D22,D23,D33] are returned) (default:
%         false).
%      chunksz: the LS problem reduces to the product of the dwi signal
%         by an inverse matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the LS
%         inverse at each step, hence taking advantage of matlab's
%         capabilities (default: 100). Note this is not possible when WLS
%         are used.
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

%%% -----------------------------------------------------------------------
% Check the mandatory input argments:
if(nargin<2)
    error('At lest the signal and the directions table must be supplied');
end
[M,N,P,G] = size(signal);
NV = M*N*P; % Total number of voxels to be processed
assert(ismatrix(gi),'gi must be 2-d matlab matrix');
assert(size(gi,1)==G,'The number of rows in gi must match the 4-th dimension of dwi');
assert(size(gi,2)==3,'The gradients table gi must have size Gx3');
%%% -----------------------------------------------------------------------
% Parse the optional input arguments:
opt.wls = false;        optchk.wls = [true,true];     % always 1x1 boolean
opt.wlsit = 1;          optchk.wlsit = [true,true];   % always 1x1 double
opt.wsc = 1.0e6;        optchk.wsc = [true,true];     % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
opt.unroll = false;     optchk.unroll = [true,true];  % always 1x1 boolean
opt.lambda = 0.006;     optchk.lambda = [true,true];  % always 1x1 double
opt.chunksz = 100;      optchk.chunksz = [true,true]; % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%% -----------------------------------------------------------------------
% Now, check bi if necessary:
if(opt.wls)
    assert(size(bi,2)==1,'bi must be a column vector if wls is true');
    assert(size(bi,1)==G,'bi must be Gx1 if wls is true');
end
%%% -----------------------------------------------------------------------
% Compute the LS matrix for tensor fitting:
GR = [ gi(:,1).*gi(:,1), 2*gi(:,1).*gi(:,2), 2*gi(:,1).*gi(:,3), ...
    gi(:,2).*gi(:,2), 2*gi(:,2).*gi(:,3), gi(:,3).*gi(:,3) ];
LS = ((GR')*GR)\(GR'); % 6xG
LS = LS'; %Gx6
%%% -----------------------------------------------------------------------
% Now, fit the data chunk-by-chunk via least squares where the mask is true
signal  = reshape(signal,[NV,G]); % NVxG
mask    = opt.mask(:);            % NVx1
% Mask...
signal = signal(mask,:); % PVxG
PV     = size(signal,1);
tmask  = zeros(PV,6); % NVx6
for ck=1:ceil(PV/opt.chunksz)
    idi = (ck-1)*opt.chunksz+1;
    idf = min(ck*opt.chunksz,PV);
    tmask(idi:idf,:) = signal(idi:idf,:)*LS; % (chunksz x G) * (G x 6) -> (chunksz x 6)
end
%%% -----------------------------------------------------------------------
% In case iterative WLS is required, this is just a first iteration. We
% need to iterate to get the solution, which requires arranging and
% inverting a different WLS matrix at each voxel. Using parallel pools
% sounds a good idea here:
use_parallel = use_parallel_test;
wsc  = opt.wsc;
ADC0 = opt.ADC0;
if(opt.wls)
    is_broadcast_available = is_broadcast_available_test;
    for it=1:opt.wlsit % Each iteration
        if(use_parallel) % Parallel pool available
            parfor pv=1:PV
                tmask(pv,:) = signal2dti_wls( tmask(pv,:), signal(pv,:), GR, bi, wsc, ADC0, is_broadcast_available );
            end
        else             % Parallel pool not-available
            for pv=1:PV
                tmask(pv,:) = signal2dti_wls( tmask(pv,:), signal(pv,:), GR, bi, wsc, ADC0, is_broadcast_available );
            end
        end
    end
end
%%% -----------------------------------------------------------------------
% Cast the result to the proper size:
tensor         = zeros(NV,6); % NVx6
tensor(mask,:) = tmask;
tensor         = reshape(tensor,[M,N,P,6]);
if(opt.unroll)
    tensor = tensor(:,:,:,[1,2,3,2,4,5,3,5,6]);
    tensor = reshape(tensor,[M,N,P,3,3]);
end
%%% -----------------------------------------------------------------------

function tensor = signal2dti_wls( tensor, signal, GR, bi, wsc, ADC0, is_broadcast_available )
%%% -------------
% tensor: 1x6
% signal: 1xG
% GR:     Gx6
% bi:     Gx1
%%% -------------
% Compute the estimated noise-free signal for the current estimate of the
% diffusion tensor:
atti = GR*(tensor');   % G x 1
atti = max(atti,0);    % G x 1
atti = min(atti,ADC0); % G x 1
atti = exp(-atti.*bi); % G x 1
%%% -------------
% Compute the WLS weights:
atti = atti.*atti;     % G x 1
MSC  = wsc*max(atti);
atti(atti<MSC) = MSC;
%%% -------------
% Compute the WLS matrix:
if(is_broadcast_available)
    atti = atti.*GR;                     % G x 6
else
    atti = bsxfun(@(x,y)(x.*y),atti,GR); % G x 6
end
GR = ((GR')*atti)\(atti');               % 6 x G
%%% -------------
% Output:
tensor = (GR*(signal'))';                % 1 x 6
%%% -------------














