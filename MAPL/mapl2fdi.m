function fdi = mapl2fdi(mapl,varargin)
%   fdi = mapl2fdi(mapl,varargin)
%
%   Computes the so-called fiber dispersion index (FDI) of a field of MAPL
%   coefficients as computed by atti2mapl, mapl. The FDI is the average 
%   disagreement of the fiber configuration inside each voxel as compared 
%   to its neighbors.
%
%   The mandatory input u1 is returned by atti2mapl.
%
%   INPUTS:
%
%      mapl: a MxNxPxK double array with the MAPL expansion coefficients.
%
%   OUPUTS:
%
%      fdi: a MxNxP volume with the fdi index at each voxel.
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
assert(ndims(mapl)==4,'Argument u1 must be MxNxPxK');
[M,N,P,K] = size(mapl);
% -------------------------------------------------------------------------
opt.nhood = [1;1;0];    optchk.nhood = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];
opt.maxthreads = 1e6;   optchk.maxthreads = [true,true];
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Compute the FDI:
fdi = coeffs2localdispersion( ...
    double(permute(mapl,[4,1,2,3])), ...
    double(opt.mask), ...
    ones(K,1), ...
    opt.nhood, ...
    opt.maxthreads   );
end
