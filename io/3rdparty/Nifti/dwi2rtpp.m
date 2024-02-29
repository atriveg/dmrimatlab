function rtpp = dwi2rtpp(dwi,gi0,bi0,tau,usetensor)
% function rtpp = dwi2rtpp(dwi,gi0,bi0,tau,usetensor)
%
%    dwi: MxNxPxG, the DWI volume
%    gi0: Gx3, the gradients table (place zero-norm for baseline)
%    bi0: 1x1 or Gx1, the b-value(s)
%    tau: 1x1, the effective diffusion time (optional; default: 20.0e-4)
%    usetensor: 1x1 logical, compute the maximum diffusion direction using
%       a tensor model? (optional; default: true, since it seems more
%       robust).

if(nargin<4)
    tau = 20.0e-3;
end
if(nargin<5)
    usetensor = true;
end

[M,N,P,G] = size(dwi);
if(numel(bi0)==1)
    bi0 = bi0*ones(G,1);
end
if( size(gi0,1)~=G || size(bi0,1)~=G )
    error('The size of the gradients table and b-values must match the 4-th dimension of the DWI');
end
if( size(gi0,2)~=3 )
    error('gi0 should have size (G+B)x3');
end
if( size(bi0,2)~=1 )
    error('bi0 should have size (G+B)x1');
end

modules = sqrt(sum(gi0.*gi0,2)); % Gx1
bsidx   = (modules<0.01);        % baselines
gridx   = ~bsidx;                % gradients

baseline = mean( dwi(:,:,:,bsidx), 4 ); % MxNxP
mask     = (baseline>10*eps);
baseline(~mask) = nan;

atti     = dwi(:,:,:,gridx);  % MxNxPxG
modules  = modules(gridx,:);  % Gx1
gi       = gi0(gridx,:);      % Gx3

atti = atti./baseline;      % MxNxPxG
gi   = gi./modules;         % Gx3
bi   = bi0(gridx).*modules; % Gx1

atti(atti<10*eps)   = 10*eps;
atti(atti>1-10*eps) = (1-10*eps);

G  = size(bi,1);
Di = -log(atti)./reshape(bi,[1,1,1,G]); % MxNxPxG

lambda = 0.001;
L2     = floor((sqrt(1+8*G)-3)/2);
SH2    = signal2sh( Di, gi, 'L', L2, 'lambda', lambda );
sig2   = sh2signal( SH2, gi );
if(usetensor)
    SH1  = signal2sh( Di, gi, 'L', 2, 'lambda', lambda );
    sig1 = sh2signal( SH1, gi );
else
    sig1 = sig2;
end

[~,id] = max(sig1,[],4);
sig2   = reshape(sig2,[M*N*P,G]);
id     = reshape(id,[M*N*P,1]);
id2    = ((1:M*N*P)') + (M*N*P)*(id-1);
sig2   = 4*sqrt(pi*tau)*sqrt(max(sig2(id2),0));

rtpp   = 1./sig2;
rtpp(sig2<10*eps) = nan;
rtpp   = reshape(rtpp,[M,N,P]);

% -------------------------------------------------------------------------
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
%         capabilities (default: 100).
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
opt.chunksz = 100;      optchk.chunksz = [true,true]; % always 1x1 double
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

% -------------------------------------------------------------------------
function signal = sh2signal( SH, gi, varargin )
% function signal = sh2signal( SH, gi, 'opt1', value1, 'opt2', value2, ... )
%
%   Computes a symmetric signal defined over the unit sphere (signal(gi) = 
%   signal(-gi)) from its SH coefficients at given directions gi:
%
%      SH: a MxNxPx((L+1)(L+2)/2) double array containing the coefficients
%         of the spherical harmonics expansionn at each voxel. Note L is
%         even, and the fourth dimension of this array must match a proper
%         size.
%      gi: a Gx3 matrix with the directions table, each row corresponding 
%         to a unit vector.
%
%      signal: a MxNxPxG double array with the signal at each voxel and for 
%         each direction.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      chunksz: the computation reduces to the product of the SH coeffs
%         by a matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the 
%         corresponding matrix, hence taking advantage of matlab's
%         capabilities (default: 100).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% Check the mandatory input arguments:
if(nargin<2)
    error('At lest the coefficients volume and the gradient table must be supplied');
end
[M,N,P,K] = size(SH);
NV = M*N*P;                % Total number of voxels to be processed
L = (sqrt(8*(K-1)+9)-3)/2; % SH order
if(abs(L-round(L))>1.0e-9)
    error('Weird size of the SH volume. Its fourth dimension should have size 1, 6, 15, 28, 45, ..., (L+1)(L+2)/2, with even L');
end
if(~ismatrix(gi))
    error('gi must be a 2-d matlab matrix');
end
if(size(gi,2)~=3)
    error('The directions table gi must have size Gx3');
end
G = size(gi,1);

% Parse the optional input arguments:
opt.chunksz = 100;      optchk.chunksz = [true,true]; % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Compute the LS matix:
B = GenerateSHMatrix( L, gi ); % GxK, where K = (L+1)*(L+2)/2
B = B'; % KxG, for convenience, see loop below

% Now, process the data chunk-by chunk where the mask is true:
SH   = reshape(SH,[NV,K]);   % NVxK
mask = opt.mask(:);          % NVx1
% Mask...
SH      = SH(mask,:); % PVxK
PV      = size(SH,1);
sigMask = zeros(PV,G); % PVxG
for ck=1:ceil(PV/opt.chunksz)
    idi = (ck-1)*opt.chunksz+1;
    idf = min(ck*opt.chunksz,PV);
    sigMask(idi:idf,:) = SH(idi:idf,:)*B; % (chunksz x K) * (K x G) -> (chunksz x G)
end
% Cast the result to the proper size:
signal = zeros(NV,G); % NVxG
signal(mask,:) = sigMask;
signal = reshape(signal,[M,N,P,G]);

