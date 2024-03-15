function xyz = dti2xyz( tensor, varargin )
% function xyz = dti2xyz( tensor, 'opt1', value1, 'opt2', value2, ... )
%
%   Computes the direction of maximum diffusion within a tensor field by
%   just computing the eigenvector associated to the maximum eigenvalue. It
%   returns a vector field of unit-norm vectors as returned by Matlab's eig
%   function.
%
%      tensor: a MxNxPx6 double array containing the unique coefficients
%         of the diffusion tensor at each voxel: D11, D12, D13, D22, D23,
%         D33.
%
%      xyz: a MxNxPx3 double array with the unit vector providing the
%         maximum diffusion direction at each vooxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with [0,0,1].

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
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});

xyz    = zeros(NV,3);            % NV x 3
mask   = opt.mask(:);            % NV x 1
tensor = reshape(tensor,[NV,6]); % NV x 6
xyz(:,3)    = 1;                              % NV  x 3
xyz(mask,:) = dti2xyz_unroll(tensor(mask,:)); % NVM x 3
xyz         = reshape(xyz,[M,N,P,3]);

function xyzn = dti2xyz_unroll(tensorn)
use_parallel = use_parallel_test;
N    = size(tensorn,1);
xyzn = zeros(N,3);

pbad = any( isinf(tensorn) | isnan(tensorn) , 2 );
tensorn(pbad,:) = 0;

if(use_parallel)
    parfor n=1:N
        D = tensorn(n,:);
        D = [ ...
            D(1),D(2),D(3);
            D(2),D(4),D(5);
            D(3),D(5),D(6) ];
        [U,~] = eig(D);
        xyzn(n,:) = U(:,3)';
    end
else
    for n=1:N
        [U,~] = eig( [ ...
            tensorn(n,1), tensorn(n,2), tensorn(n,3);
            tensorn(n,2), tensorn(n,4), tensorn(n,5);
            tensorn(n,3), tensorn(n,5), tensorn(n,6)   ] );
        xyzn(n,:) = U(:,3)';
    end
end
