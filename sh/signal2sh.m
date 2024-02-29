function SH = signal2sh( signal, gi, varargin )
% function SH = signal2sh( signal, gi, 'opt1', value1, 'opt2', value2, ... )
%
%   Takes a generic symmetric signal defined over the unit sphere and 
%   represents it in the basis of spherical harmonics by using (Tikhonov 
%   regularized) least squares:
%
%      signal: a MxNxPxG double array containing the signal sampled at G
%         directions at each voxel within the MxNxP image frame.
%      gi: a Gx3 matrix with the directions sampled table, each row 
%         corresponding to a unit vector (signal(gi) is assumed to be equal
%         to signal(-gi));
%
%      SH: a MxNxPx((L+1)(L+2)/2) array with the coefficients computed at 
%         each voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      L: an even integer with the maximum order of the SH to be used
%         (default: 6).
%      lambda: the Tikhonov regularization parameter for the linear least
%         squares problem (default 0.006).
%      chunksz: the LLS problem reduces to the product of the dwi signal
%         by an inverse matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the LLS
%         inverse at each step, hence taking advantage of matlab's
%         capabilities (default: 1000).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% Check the mandatory input argments:
if(nargin<2)
    error('At lest the signal and the directions table must be supplied');
end
[M,N,P,G] = size(signal);
NV = M*N*P; % Total number of voxels to be processed
if(~ismatrix(gi))
    error('gi must be 2-d matlab matrix');
end
if(size(gi,1)~=G)
    error('The number of rows in gi must match the 4-th dimension of dwi');
end
if(size(gi,2)~=3)
    error('The gradients table gi must have size Gx3');
end

% Parse the optional input arguments:
opt.L = 6;              optchk.L = [true,true];       % always 1x1 double
opt.lambda = 0.006;     optchk.lambda = [true,true];  % always 1x1 double
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Compute the LS matix for SH fitting:
B   = GenerateSHMatrix( opt.L, gi );    % GxK, where K=(L+1)(L+2)/2
LR  = GenerateSHEigMatrix( opt.L );     % KxK
WLS = (B'*B+(opt.lambda).*LR^2)\(B');   % (KxK)^(-1) * (KxG) -> KxG
WLS = WLS'; % GxK, for convenience, see loop below

% Now, fit the data chunk-by-chunk via least squares where the mask is true
signal  = reshape(signal,[NV,G]); % NVxG
mask    = opt.mask(:);            % NVx1
% Mask...
signal = signal(mask,:); % PVxG
PV     = size(signal,1);
SHmask = zeros(PV,size(WLS,2)); % NVxK
for ck=1:ceil(PV/opt.chunksz)
    idi = (ck-1)*opt.chunksz+1;
    idf = min(ck*opt.chunksz,PV);
    SHmask(idi:idf,:) = signal(idi:idf,:)*WLS; % (chunksz x G) * (G x K) -> (chunksz x K)
end
% Cast the result to the proper size:
SH   = zeros(NV,size(WLS,2)); % NVx(L+1)(L+2)/2
SH(mask,:) = SHmask;
SH = reshape(SH,[M,N,P,size(WLS,2)]);

