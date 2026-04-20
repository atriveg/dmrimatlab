function fdi = shodf2fdi(sh,varargin)
%   fdi = shodf2fdi(sh,varargin)
%
%   Computes the so-called fiber dispersion index (FDI) of a field of ODFs
%   represented with spherical harmonics, sh. The FDI is the average 
%   disagreement of the fiber configuration inside each voxel as compared 
%   to its neighbors.
%
%   INPUTS:
%
%      sh: a MxNxPxK double array with the SH coefficients of an ODF.
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
assert(ndims(sh)==4,'Argument u1 must be MxNxPxK');
[M,N,P,K] = size(sh);
L = (sqrt(8*(K-1)+9)-3)/2; % SH order
if(abs(L-round(L))>1.0e-9)
    error('Weird size of the SH volume. Its fourth dimension should have size 1, 6, 15, 28, 45, ..., (L+1)(L+2)/2, with even L');
end
% -------------------------------------------------------------------------
opt.nhood = [1;1;0];    optchk.nhood = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];
opt.maxthreads = 1e6;   optchk.maxthreads = [true,true];
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Compute the FDI:
fdi = coeffs2localdispersion( ...
    double(permute(sh,[4,1,2,3])), ...
    double(opt.mask), ...
    ones(K,1), ...
    opt.nhood, ...
    opt.maxthreads   );
end
