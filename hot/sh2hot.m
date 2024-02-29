function [HOT,A] = sh2hot( SH, varargin )
% function [HOT,A] = sh2hot( SH, 'opt1', value1, 'opt2', value2, ... )
%
%   Computes the (unique) components of a hyper-symmetric, higher order
%   tensor of rank L from the SH expansion of a signal up to order L.
%
%      SH: a MxNxPx((L+1)(L+2)/2) double array containing the coefficients
%         of the spherical harmonics expansion at each voxel. Note L is
%         even, and the fourth dimension of this array must match a proper
%         size.
%      HOT: a MxNxPx((L+1)(L+2)/2) double array containing the unique 
%         coefficients of the equivalent HOT representation.
%      A: if a second output is requested, A is the conversion matrix
%         leading from SH to HOT (Size KxK, for K=(L+1)(L+2)/2).
%
%   For an order L, the rank-L HOT represetation is the sum of
%   K=(L+1)(L+2)/2 terms:
%
%
%      S(x,y,z) = sum_k T(k)*mu(k)*x^nx(k)*y^ny(k)*z^nz(k)
%               = T(1)*mu(1)*x^L*y^0*z^0 + ...
%               + T(2)*mu(2)*x^{L-1}*y^1*z^0 + ...
%               + T(3)*mu(3)*x^{L-2}*y^2*z^0 + ...
%               ...
%               + T(L)*mu(L)*x^{0}*y^L*z^0 + ...
%               + T(L+1)*mu(L+1)*x^{L-1}*y^0*z^1 + ...
%               ...
%               + T(K)*mu(K)*x^0*y^0*z^L + ...
%
%   where mu(k) are the multiplicities of each unique component of the HOT
%   and nx+ny+nz=L.
%
%   See:
%
%     M. Descoteaux, E. Angelino, S. Fitzgibbons, and R. Deriche. 
%     "Apparent Diffusion Profile estimation from High Angular 
%     Resolution Diffusion Images: estimation and applications."
%     Magnetic Resonance in Medicine, 56(2):395–410, 2006.
%
%   for further details.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      maxthreads: ONLY IN POSIX SYSTEMS the algorithm is run with
%         multiple threads. This is the maximum allowed number of threads,
%         which can indeed be reduced if it exceeds the number of logical
%         cores (default: the number of logical cores in the machine).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the mandatory input arguments:
if(nargin<1)
    error('At lest the SH coefficients volume must be supplied');
end
SH = double(SH);
[M,N,P,K] = size(SH);
NV = M*N*P;                % Total number of voxels to be processed
L = (sqrt(8*(K-1)+9)-3)/2; % SH order
if(abs(L-round(L))>1.0e-9)
    error('Weird size of the SH volume. Its fourth dimension should have size 1, 6, 15, 28, 45, ..., (L+1)(L+2)/2, with even L');
end
% Parse the optional input arguments:
opt.mask       = true(M,N,P); optchk.mask       = [true,true];    % boolean with the size of the image field
opt.maxthreads = 1.0e6;       optchk.maxthreads = [true,true];    % always 1x1 pos number
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH   = reshape(SH,[NV,K]);   % NVxK
mask = opt.mask(:);          % NVx1
SH   = SH(mask,:); % PVxK
[HOT2,A] = sh2hot_(SH,opt.maxthreads);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HOT = zeros(NV,K);
HOT(mask,:) = HOT2;
HOT = reshape(HOT,[M,N,P,K]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
