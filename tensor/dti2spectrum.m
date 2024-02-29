function [u1,u2,u3,l1,l2,l3] = dti2spectrum( tensor, varargin )
% function [u1,u2,u3,l1,l2,l3] = dti2spectrum( tensor, 'opt1', value1, 'opt2', value2, ... )
%
%   Computes the eigenvectors and eigenvalues of a diffusion tensor field
%   in descending order:
%
%      tensor: a MxNxPx6 double array containing the unique coefficients
%         of the diffusion tensor at each voxel: D11, D12, D13, D22, D23,
%         D33.
%
%      u1: a MxNxPx3 double array with the eigenvector associated to the
%         largest eigenvalue, l1.
%      u2: a MxNxPx3 double array with the eigenvector associated to the
%         second largest eigenvalue, l2.
%      u3: a MxNxPx3 double array with the eigenvector associated to the
%         smallest eigenvalue, l3.
%      l1: a MxNxP double array with the largest eigenvalue.
%      l2: a MxNxP double array with the second largest eigenvalue.
%      l3: a MxNxP double array with the smallest eigenvalue.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      maxthreads: the maximum number of threads to use in POSIX systems
%         (default: automatically determined).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with [0,0,1].

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<1)
    error('At lest the tensor volume must be supplied');
end
[M,N,P,K] = size(tensor);
NV = M*N*P; % Total number of voxels to be processed
if(K~=6)
    error('Weird size of the tensor volume. Its fourth dimension should have size 6');
end
% Parse the optional input arguments:
opt.maxthreads = 1e6;   optchk.maxthreads = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
u1 = zeros(NV,3); u1(:,1) = 1; % NV x 3
u2 = zeros(NV,3); u2(:,2) = 1; % NV x 3
u3 = zeros(NV,3); u3(:,3) = 1; % NV x 3
l1 = zeros(NV,1);              % NV x 1
l2 = zeros(NV,1);              % NV x 1
l3 = zeros(NV,1);              % NV x 1
% -------------------------------------------------------------------------
mask   = opt.mask(:);            % NV x 1
tensor = reshape(tensor,[NV,6]); % NV x 6
tensor = tensor(mask,:);         % NVM x 6
% -------------------------------------------------------------------------
[l1m,l2m,l3m,u1m,u2m,u3m] = dti2spectrum_(tensor',opt.maxthreads);
% -------------------------------------------------------------------------
u1(mask,:) = u1m';
u2(mask,:) = u2m';
u3(mask,:) = u3m';
l1(mask)   = l1m';
l2(mask)   = l2m';
l3(mask)   = l3m';
% -------------------------------------------------------------------------
u1 = reshape(u1,[M,N,P,3]);
u2 = reshape(u2,[M,N,P,3]);
u3 = reshape(u3,[M,N,P,3]);
l1 = reshape(l1,[M,N,P]);
l2 = reshape(l2,[M,N,P]);
l3 = reshape(l3,[M,N,P]);
% -------------------------------------------------------------------------
