function [S,mu,pwrs] = hot2signal( HOT, gi, varargin )
% function [S,mu,pwrs] = hot2signal( HOT, gi, 'opt1', value1, 'opt2', value2, ... )
%
%   Computes the values of a volume of Higher Order Tensors at given space
%   orientations gi.
%
%      HOT: a MxNxPx((L+1)(L+2)/2) double array containing the unique 
%         coefficients of a hyper-symmetric HOT representation. Note L is
%         even, and the fourth dimension of this array must match a 
%         proper size.
%      gi: a Gx3 array, each row representing a space orientation (they
%         should be unit norm) for which the HOT will be evaluated.
%      S: a MxNxPxG, the values of the HOT at each voxel.
%      mu: ((L+1)(L+2)/2)x1, the multiplicities of each unique term of the
%         HOT (integer values).
%      pwrs: ((L+1)(L+2)/2)x3, each column being the powers of each of the
%         x, y, and z components of the space orientations associated to
%         the coefficients of the HOT, fulfilling sum(j,:)=L for all j.
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
%   and nx+ny+nz=L. Then:
%
%       S(:,:,:,g) = ...
%            sum( HOT.*reshape(mu,[1,1,1,K]).* ...
%                gi(g,1)^reshape(pwrs(:,1),[1,1,1,K]).* ...
%                gi(g,2)^reshape(pwrs(:,2),[1,1,1,K]).* ...
%                gi(g,3)^reshape(pwrs(:,3),[1,1,1,K]).*,     4    );
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
if(nargin<2)
    error('At lest the coefficients volume and the gradient table must be supplied');
end
HOT = double(HOT);
[M,N,P,K] = size(HOT);
NV = M*N*P;                % Total number of voxels to be processed
L = (sqrt(8*(K-1)+9)-3)/2; % SH order
if(abs(L-round(L))>1.0e-9)
    error('Weird size of the HOT volume. Its fourth dimension should have size 1, 6, 15, 28, 45, ..., (L+1)(L+2)/2, with even L');
end
if( ~ismatrix(gi) || (size(gi,2)~=3) )
    error('Second argument, gi, must be Gx3');
end
G = size(gi,1);
% Parse the optional input arguments:
opt.mask       = true(M,N,P); optchk.mask       = [true,true];    % boolean with the size of the image field
opt.maxthreads = 1.0e6;       optchk.maxthreads = [true,true];    % always 1x1 pos number
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HOT  = reshape(HOT,[NV,K]);   % NVxK
mask = opt.mask(:);          % NVx1
HOT  = HOT(mask,:); % PVxK
[S2,mu,pwrs] = hot2signal_(HOT,gi,opt.maxthreads);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = zeros(NV,G);
S(mask,:) = S2;
S = reshape(S,[M,N,P,G]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
