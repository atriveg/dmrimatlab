function [paths,pcosts,dirs,stopconds] = backTracing( costs, dirs, targets, varargin )
% function [paths,pcosts,dirs,stopconds] = backTracing( costs, dirs, targets, ...
%                               'option1', value1, 'option2', value2, ... )
%
% According to the maps of global costs "costs" and directional costs
% "dirs", as obtained with fastSweeping, traces back streamlines from the
% target points described by "targets" towards the seeding points
% originally used in the fastSweeping function. Note the Field of View and
% the ijk2xyz/xyz2ijk matrices used in both functions must be compatible
% with each other.
%
%   MANDATORY INPUTS:
%
%      costs: a X_1 x X_2 x ... x X_D array of doubles representing the
%         costs of arriving at each point in the FoV from the seeding
%         points, which is returned by the fastSweeping function.
%      dirs: X_1 x X_2 x ... x X_D x D array of doubles with the arrival
%         directions the minimal costs are propagated along. It is roughly
%         the (negative) gradient of costs. Note these dirs correspond to
%         ANATOMICAL, NOT pixel directions. It is also returned by the
%         fastSweeping function.
%      targets: D x N, an array of N target points (each column is a point)
%         IN THE XYZ OR ANATOMICAL SPACE, each one is a target from which
%         backtracing will be started towards the seeding region. All of
%         them should be inside the FoV of the image (see optional argument
%         ijk2xyz for a note on this).
%
%   OUTPUTS:
%
%      paths: a N x 1 cell array, each entry with size D x M_n describing a
%         streamline traced back from each target point (length might be 0
%         if backtraing fails for some reason). Within each D x M_n
%         streamline, each column represents a successive point in the
%         streamline.
%      pcosts: a N x 1 cell array, each entry with size 1 x M_n describing
%         the interpolated costs (from the "costs" input) reached at each
%         successive point in the respective streamline (hence, it should
%         roughly contain monotonically decreasing values).
%      dirs: this is the same as the input argument dirs, but with values
%         corrected at mask boundaries that help to push the streamlines
%         far from non-masked values to prevent early (undesired) stop
%         conditions.
%      stopconds: a 1 x N double array which stores the condition that
%         caused each streamline to stop being traced back, according to
%         the following codes:
%
%           1: normal condition, reached seeding region (low cost reached)
%           2: normal condition, reached seeding region (low jump reached)
%           3: streamline reached a point outside the FoV
%           4: streamline reached a point outside the mask
%           5: maximum allowed number of points per stream was reached
%           6: maximum streamline length reached
%           7: local curvature radius below allowed minimum 
%
%   OPTIONAL KEY/VALUE ARGUMENTS:
%
%      mask: a X_1 x X_2 x ... x X_D boolean array with a mask. Should be 
%         the same previously used with fastSweeping (default: all true).
%      ijk2xyz: a (D+1)x(D+1) matrix with the homogeneous coordinates 
%         transform that converts pixel indices to physical indices. For
%         example, for D=2, it should be something like:
%            ijk2xyz = [ cos(a)*s_i   -sin(a)*s_j   o_x ]
%                      [ sin(a)*s_i    cos(a)*s_j   o_y ]
%                      [          0             0     1 ]
%         where s_i and s_j are the resolutions (in millimeters) each pixel
%         dimension has, a is the rotation of the pixel directions w.r.t
%         the anatomical axes and o_x, o_y is the anatomical origin of the
%         image. Should be the same used with fastSweeping
%         (default: eye(D+1)).
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
%         effectively reduced once we approach seeding points, and also
%         increased under rather unusual conditions) (default: 0.5, anatom-
%         ical units).
%      Mlength: 1 x 1, the maximum length allowed for backtraced
%         streamlines (default: 800, in anatomical units).
%      mcurv: 1 x 1, the minimum curvature radius allowed within a
%         streamline before backtracing is aborted (default: 0.8, in
%         anatomical units).
%      prune: a 1 x P vector, with P=0..7 that tells which streamlines
%         have to be pruned from the computation because of certain
%         stopping conditions. For example, if prune = [3,4,7], those
%         streamlines that terminated because either:
%           - A point laid outisde the FoV
%           - A point laid outisde the mask
%           - A point produced a very small curvature radius
%         will be replaced with empty arrays (default: [3,4]).
%      maxthreads: the algorithm is run as a multi-threaded mex. This is 
%         the maximum allowed number of threads, which can indeed be 
%         reduced if it exceeds the number of logical cores (default: the 
%         number of logical cores in the machine).

% -------------------------------------------------------------------------
% Check the mandatory input argments:
if(nargin<3)
    error('At lest the costs, the directions costs, and the targets must be supplied');
end
ndim = ndims(costs);
fov  = size(costs);
assert( isequal(size(dirs),[fov,ndim]), 'Unmatching sizes of costs and dirs' );
assert( ismatrix(targets), 'Input targets should be a D x N matrix' );
assert( size(targets,1)==ndim, 'The # of rows of variable targets should match the # of dimensions of costs' );
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.mask = true(fov);      optchk.mask = [true,true];
opt.ijk2xyz = eye(ndim+1); optchk.ijk2xyz = [true,true];
opt.step = 0.5;            optchk.step = [true,true];
opt.Mlength = 800;         optchk.Mlength = [true,true];
opt.mcurv = 0.8;           optchk.mcurv = [true,true];
opt.prune = [3,4];         optchk.prune = [true,false]; % Variable size
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
% -------------------------------------------------------------------------
% Tweak the values of the arrival directions at the borders of the mask, so
% that they point towards the region inside the mask and interpolation
% issues do not appear:
mask = double(opt.mask);
if( any(mask(:)<0.5) ) % There is, indeed, a mask
    % -----
    gmask = zeros([fov,ndim]);
    gmask = reshape( gmask, [prod(fov),ndim] );
    % -----
    krn   = [1;0;-1]/2;
    for d=1:ndim
        dperm      = 1:ndim;
        dperm(d)   = 1;
        dperm(1)   = d;
        krnd       = permute( krn, dperm );
        dmask      = convn( mask, krnd, 'same' );
        gmask(:,d) = dmask(:);
    end
    % -----
    gmask = [ gmask, zeros(size(gmask,1),1) ];
    gmask = gmask*(opt.ijk2xyz)';
    gmask = gmask(:,1:end-1);
    % -----
    nrm   = sqrt(sum(gmask.*gmask,2));
    gmask(nrm>eps(1),:)  = gmask(nrm>eps(1),:)./nrm(nrm>eps(1),1);
    gmask(nrm<=eps(1),:) = 0;
    % -----
    dirs                   = reshape( dirs, [prod(fov),ndim] );
    dirs( mask(:)<0.5, : ) = gmask( mask(:)<0.5, : );
    % -----
    dirs = reshape( dirs, [fov,ndim] );
end
% -------------------------------------------------------------------------
% Reshape the arrival directions volume to make it compatible with the mex
% function:
dirs2 = permute( dirs, [ndim+1,1:ndim] );
% -------------------------------------------------------------------------
% Invert the ijk2xyz matrix, since the mex function needs both:
ijk2xyz = double(opt.ijk2xyz);
xyz2ijk = ijk2xyz\eye(ndim+1);
% -------------------------------------------------------------------------
% Directly call the mex implementation:
[paths,pcosts,stopconds] = backTracing_( double(costs), double(dirs2), ...
    mask, double(targets), ijk2xyz, xyz2ijk, ...
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
        paths(idxs,1)  = {zeros(ndim,0)};
        pcosts(idxs,1) = {zeros(1,0)};
    end
end
