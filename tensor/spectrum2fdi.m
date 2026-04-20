function fdi = spectrum2fdi(u1,varargin)
%   fdi = spectrum2fdi(u1,varargin)
%
%   Computes the so-called fiber dispersion index (FDI) of a diffusion 
%   tensor image based on its main eigenvector, u1, as computed by the 
%   dti2spectrum function. The FDI is the average disagreement of the fiber
%   configuration inside each voxel as compared to its neighbors.
%
%   The mandatory input u1 is returned by dti2spectrum.
%
%   INPUTS:
%
%      u1: a MxNxPx3 double array with the unit-norm eigenvector associated
%         to the largest eigenvalue, l1.
%
%   OUPUTS:
%
%      fdi: a MxNxP volume with the fdi index at each voxel.
%
%
%   For example, you can plot the 17-th slice of a "tensor" volume as:
%
%      >> [u1,~,~,~,~,~] = dti2spectrum( tensor, 'mask', mask );
%      >> fdi = spectrum2fdi( u1, 'mask', mask, 'nhood', [1;1;0] );
%      >> slice = 17;
%      >> imshow(fdi(:,:,slice));
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      nhood: a 3x1 array with the radii of the 3-D neighborhood used for
%         local computations (default: [1;1;0] ).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      maxthreads: the maximum number of threads to use in the core mex
%         function (default: Inf, will be adjusted to the actual number of
%         available cores).

% -------------------------------------------------------------------------
assert(ndims(u1)==4,'Argument u1 must be MxNxPx3');
[M,N,P,Q] = size(u1);
assert(Q==3,'Argument u1 must be MxNxPx3');
% -------------------------------------------------------------------------
opt.nhood = [1;1;0];    optchk.nhood = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];
opt.maxthreads = 1e6;   optchk.maxthreads = [true,true];
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Compute the FDI:
fdi = coeffs2localdispersion( ...
    double(permute(u1,[4,1,2,3])), ...
    double(opt.mask), ...
    [1;1;1], ...
    opt.nhood, ...
    opt.maxthreads   );
end
