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
%         capabilities (default: 1000).
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
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
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
