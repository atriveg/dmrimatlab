function [paths,pvals,stopconds] = dtiTractography( dti, seeds, varargin )
% function [paths,pvals,stopconds] = dtiTractography( dti, seeds, ...
%                               'option1', value1, 'option2', value2, ... )
%
% Computes standard DTI-based tractography by integrating the vector field
% described by the principal eigenvector of the diffusion tensor at each
% voxel. Starting from a set of seeding points, the algorithm proceeds by
% tracing the streamlines in both directions using a RK4 numerical
% integration scheme until a given stopping criterion (tipically, the FA
% falling below a minimum threshold) is met.
%
%   MANDATORY INPUTS:
%
%      dti: a X x Y x Z x 6 array of doubles representing the unique
%         components of the diffusion tensor as returned by atti2dti.
%      seeds: 3 x N, an array of N seeding points (each column is a point)
%         IN THE XYZ OR ANATOMICAL SPACE, each one is a seed from which
%         the tractography will be started. All of them should be inside
%         the FoV of the image (see optional argument ijk2xyz for a note on
%         this). Note that two streamlines will appear for each seeding
%         point, one per direction (positive or negative) of the main
%         eigenvector at the seeding point.
%
%   OUTPUTS:
%
%      paths: a 2N x 1 cell array, each entry with size 3 x M_n describing
%         a streamline traced from each each point (length might be 0
%         if fiber tracking fails for some reason). Within each 3 x M_n
%         streamline, each column represents a successive point in the
%         streamline (in the XYZ, ANATOMICAL SPACE).
%      pvals: a 2N x 1 cell array, each entry with size 1 x M_n describing
%         the interpolated values of the scalar map used to check the
%         stopping condition (tipically, the FA) reached at each succesive
%         point in the respective streamline.
%      stopconds: a 1 x 2N double array which stores the condition that
%         caused each streamline to stop being fiber tracking, according to
%         the following codes:
%
%           1: normal condition, reached the threshold in the scalar map
%           2: the jump became too small due to interpolation effects
%           3: streamline reached a point outside the FoV
%           4: streamline reached a point outside the mask
%           5: maximum allowed number of points per stream was reached
%           6: maximum streamline length reached
%           7: local curvature radius below allowed minimum 
%
%   OPTIONAL KEY/VALUE ARGUMENTS:
%
%      mask: a X x Y x Z boolean array with a mask (default: all true).
%      scalar: a string either 'fa' or 'cl', telling the algorithm stopping
%         criterion should be based on either the FA (Fractional
%         Anisotropy) or Westin's Cl (Linear Coefficient) (default: fa).
%      threshold: the threshold to apply to the scalar map. If the value of
%         the FA/Cl falls below this threshold, the fiber tracking is
%         stopped (default: 0.4).
%      ijk2xyz: a 4x4 matrix with the homogeneous coordinates transform
%         that converts pixel indices to physical indices. Something like:
%            ijk2xyz = [ R11*s1   R12*s2   R13*s3   o_x ]
%                      [ R21*s1   R22*s2   R23*s3   o_y ]
%                      [ R31*s1   R32*s2   R33*s3   o_z ]
%                      [      0        0        0     1 ]
%         where s1 to S3 are the resolutions (in millimeters) each voxel
%         dimension has, [Rij] is the rotation matrix that aligns the voxel
%         orientations to the anatomical orientations and o_x, o_y, o_z is 
%         the anatomical origin of the image (default: eye(D+1)).
%         NOTE: The ijk indices are assumed to be zero-based (like in any
%         other software for medical imaging), so if you're planning to
%         retrieve pixel indices from the anatomical coordinates returned
%         by this function, something like this:
%               >> streamline = paths{1};
%               >> streamline(4,:) = 1; % Homogeneous coords
%               >> streamlineijk = (ijk2xyz\(streamline'))';
%               >> streamlineijk = streamlineijk(1:3,:);
%         Don't forget to add 1 to each dimension so that these indices
%         match Matlab's one-based indexing:
%               >> streamlineijk = streamlineijk + 1;
%      step: 1 x 1, the (base) step length used for Runge-Kutta method
%         integration of the vector field of directional costs (might be
%         effectively reduced near non-masked or out-of-bounds voxels)
%         (default: 0.5, anatomical units).
%      Mlength: 1 x 1, the maximum length allowed for tracked streamlines
%         (default: 800, in anatomical units).
%      mcurv: 1 x 1, the minimum curvature radius allowed within a
%         streamline before backtracing is aborted (default: 0.8, in
%         anatomical units).
%      prune: a 1 x P vector, with P=0..8 that tells which streamlines
%         have to be pruned from the computation because of certain
%         stopping conditions. For example, if prune = [3,4,7], those
%         streamlines that terminated because either:
%           - A point laid outisde the FoV
%           - A point laid outisde the mask
%           - A point produced a very small curvature radius
%         will be replaced with empty arrays (default: [3,4,8]).
%      maxthreads: the algorithm is run as a multi-threaded mex. This is 
%         the maximum allowed number of threads, which can indeed be 
%         reduced if it exceeds the number of logical cores (default: the 
%         number of logical cores in the machine).

% -------------------------------------------------------------------------
% Check the mandatory input argments:
if(nargin<2)
    error('At lest the dtivolume and the seeds must be supplied');
end
[X,Y,Z,G] = size(dti);
fov  = [X,Y,Z];
assert(G==6,'The fourth dimension of dti should have size 6');
assert( ismatrix(seeds), 'Input seeds should be a 3 x N matrix' );
assert( size(seeds,1)==3, 'Input seeds should be a 3 x N matrix' );
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.mask = true(fov);      optchk.mask = [true,true];
opt.scalar = 'fa';         optchk.scalar = [true,true];
opt.threshold = 0.4;       optchk.threshold = [true,true];
opt.ijk2xyz = eye(4);      optchk.ijk2xyz = [true,true];
opt.step = 0.5;            optchk.step = [true,true];
opt.Mlength = 800;         optchk.Mlength = [true,true];
opt.mcurv = 0.8;           optchk.mcurv = [true,true];
opt.prune = [3,4,8];       optchk.prune = [true,false]; % Variable size
opt.maxthreads = 1e6;      optchk.maxthreads = [true,true];
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -----------------
% Check the "prune" argument
prune = opt.prune;
if(~isempty(prune))
    prune = unique(prune);
    prune = prune(:)';
    assert( all(prune<=7), 'prune argument should contain integer values from 1 to 7 ' );
    assert( all(prune>=1), 'prune argument should contain integer values from 1 to 7 ' );
    assert( all(abs(round(prune)-prune)<=eps(1)), 'prune argument should contain integer values from 1 to 7 ' );
end
% -----------------
% Check the "scalar" argument
scalar = lower(opt.scalar);
switch(scalar)
    case 'fa'
    case 'cl'
    otherwise
        error('Optional argument <scalar> must be either fa or cl');
end
% -------------------------------------------------------------------------
% Compute the spectrum of the tensor volume and the derived measures:
mask = ( double(opt.mask) > 0.5 );
[dirs,~,~,l1,l2,l3] = dti2spectrum( dti, 'mask', opt.mask, ...
    'maxthreads', opt.maxthreads, 'clip', true );
scmap = spectrum2scalar( l1, l2, l3, 'scalar', scalar, 'mask', opt.mask );
% -------------------------------------------------------------------------
% Reshape the directions volume to make it compatible with the mex function:
dirs = permute( dirs, [4,1,2,3] );
% -------------------------------------------------------------------------
% Invert the ijk2xyz matrix, since the mex function needs both:
ijk2xyz = double(opt.ijk2xyz);
xyz2ijk = ijk2xyz\eye(4);
% -------------------------------------------------------------------------
% Directly call the mex implementation:
[paths,pvals,stopconds] = fiberTracking_( double(scmap), double(dirs), ...
    double(mask), double(seeds), ijk2xyz, xyz2ijk, ...
    abs(double(opt.threshold)), ...
    abs(double(opt.step)), abs(double(opt.Mlength)), ...
    abs(double(opt.mcurv)), abs(double(opt.maxthreads)) );
stopconds = stopconds+1;
% -------------------------------------------------------------------------
% Prune the outputs, if necessary:
if(~isempty(prune))
    idxs = [];
    for n=1:length(prune)
        idx  = find(stopconds==prune(n));
        idxs = [idxs,idx]; %#ok<AGROW>
    end
    if(~isempty(idxs))
        paths(idxs,1) = {zeros(3,0)};
        pvals(idxs,1) = {zeros(1,0)};
    end
end
