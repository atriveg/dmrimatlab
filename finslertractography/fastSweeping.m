function [costs,dirs] = fastSweeping( dcosts, gi, seeds, varargin )
% [costs,dirs] = fastSweeping( ocosts, gi, seeds, ...
%                               'option1', value1, 'option2', value2, ... )
%
%   Implements the "Fast Sweeping" algorithm, based on dynamic programming,
%   to calculate minimal cost trajectories inside a field of directional
%   costs in D-dimensional problems. In precise terms, let phi(x,r) be a
%   directional cost defined at each point x in R^D and for each arrival
%   direction r in R^D : ||r||=1. Let x_0 be a "seeding point" (or a set of 
%   points), for which the cost is 0. The algorithm finds, for each 
%   "target point" x_T in the field of view, the minimum cost of reaching 
%   x_0, defined as the line integral:
%
%     C(x_T) = int_{x_0}^{x_T} phi( x, x'(t)/||x'(t)|| ) ||x'(t)|| dt  [1]
%
%   i.e. the minimum cost reachable for any possible trajectory x(t).
%
%   Additionally the algorithm provides, for each point in the FOV, the
%   optimal "arrival direction" from which the optimal cost was obtained
%   with a single jump, i.e. an approximation to the (negative) gradient of
%   the map of optimal costs.
%
%   MANDATORY INPUTS:
%
%      dcosts: a X_1 x X_2 x ... x X_D x N array of doubles representing a
%         discretization of phi(x,r) at N different points of the manifold
%         r in R^D : ||r||=1 (e.g. a circumference with radius 1 if D=2, or
%         a sphere with radius 1 if D=3). Note the FOV has D dimensions.
%      gi: a N x D array with the unit directions describing the
%         discretization of the orientations space. Note these directions
%         are assumed to be provided in "image coordinates", i.e. a value
%         [0,0,1] points towards the positive direction of the anatomical
%         'z'-axis, NOT the positve direction of the 3-rd pixel coordinate.
%         This implies the 'ijk2xyz' optional parameter should be taken
%         into account almost always.
%      seeds: a X_1 x X_2 x ... x X_D array with integer values describing
%         the seeding points/regions. Typically, those D-dimensional pixels
%         whose value is 1 are used as seeds, though this behavior can be
%         fixed with the optional argument 'seedval'.
%
%   OUTPUTS:
%
%      costs: a X_1 x X_2 x ... x X_D array of doubles with the optimal
%         costs; the value of each pixel minimizes eq. [1] for any possible
%         trajectory x(t). Non-reachable points are assigned an inifinite
%         cost.
%      dirs: X_1 x X_2 x ... x X_D x D array of doubles with the arrival
%         directions the minimal costs are propagated along. It is roughly
%         the (negative) gradient of costs. Note these dirs are always
%         chosen among those provided in gi, so that they correspond to
%         ANATOMICAL, NOT pixel directions.
%
%   OPTIONAL KEY/VALUE ARGUMENTS:
%
%      mask: a X_1 x X_2 x ... x X_D boolean array with a mask. Pixels
%         outside the mask are considered unreachable (default: all true).
%         NOTE: you should make sure that the mask is only made up of
%         connected components, otherwise you will get isolated islands
%         that will never be visited causing the algorithm to iterate over
%         and over until the maximum number of iterations is reached. The
%         following piece of code can remove all but the largest island
%         from a 3-D mask:
%             [objs,nobjs] = bwlabeln(mask,26);
%             szs          = zeros(1,nobjs);
%             for n=1:nobjs
%                tmp    = (objs==n);
%                szs(n) = sum(tmp(:));
%             end
%             [~,mpos] = max(szs);
%             mask     = (objs==mpos);
%      ijk2xyz: a (D+1)x(D+1) matrix with the homogeneous coordinates 
%         transform that converts pixel indices to physical indices. For
%         example, for D=2, it should be something like:
%            ijk2xyz = [ cos(a)*s_i   -sin(a)*s_j   o_x ]
%                      [ sin(a)*s_i    cos(a)*s_j   o_y ]
%                      [          0             0     1 ]
%         where s_i and s_j are the resolutions (in millimeters) each pixel
%         dimension has, a is the rotation of the pixel directions w.r.t
%         the anatomical axes and o_x, o_y is the anatomical origin of the
%         image (default: eye(D+1)).
%      seedval: an integer value indicating which pixels in the "seeds"
%         volume are actually seeds. Note that at least one seeding point
%         must be present in the volume, otherwise the algorithm will not
%         run (default: 1).
%      maxiters: the maximum number of iterations of the Fast-Sweeping
%         algorithm. Each iteration implies processing the entire FOV 2^D
%         times, for all possible combinations of causal and anticausal
%         directions of each image dimension (default: 200).
%      change: if the total change in the optimal costs is 'change' times 
%         the initial one, we consider Fast Sweeping has converged
%         (default: 1.0e-3).
%      verbose: whether or not print information after each iteration of
%         the algorithm (default: false).

% -------------------------------------------------------------------------
% Check the mandatory input argments:
if(nargin<3)
    error('At lest the directional costs volume, the directions table, and the seeds mask must be supplied');
end
ndim = ndims(dcosts)-1; % The number of image dimensions
fov  = size(dcosts);    % The Field of View of the image
N    = fov(end);        % The number of incoming directions to test
if(ndim==1)
    fov = [fov(1:end-1),1];
else
    fov = fov(1:end-1);
end
assert( isequal(size(gi),[N,ndim]), 'Dimensions mismatch for gi' );
assert( isequal(size(seeds),fov), 'The size of seeds should match the FOV of dcosts' );
% Make sure the gi are normalized:
ni = sqrt(sum(gi.*gi,2));
assert(all(ni>1.0e-3),'Some of the gi have 0-norm, cannot process'),
gi = gi./ni; % Broadcasted operation
% Reshape the directional costs volume to make it compatible with the mex
% function:
dcosts = permute( dcosts, [ndim+1,1:ndim] );
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.mask = true(fov);      optchk.mask = [true,true];
opt.ijk2xyz = eye(ndim+1); optchk.ijk2xyz = [true,true];
opt.seedval = 1;           optchk.seedval = [true,true];
opt.maxiters = 200;        optchk.maxiters = [true,true];
opt.change = 1.0e-3;       optchk.change = [true,true];
opt.verbose = false;       optchk.verbose = [true,true];
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Sanity checks:
% Check for bad values
if(any(isnan(dcosts(:))))
    warning('Input volume dcosts contains NaN values');
end
if(any(isinf(dcosts(:))))
    warning('Input volume dcosts contains Inf values');
end
% Make sure all directional costs are non-negative:
dcosts(dcosts<0) = 0.0;
% Initialize the costs to a N-D array of all -1, since this value is
% internally used to represent infinite costs:
costs = -ones(fov);
% Make sure at least one seeding point is present inside the mask:
seeds = ( opt.mask & (seeds==opt.seedval) );
assert( any(seeds(:)), ...
    'Not a single seeding point found within the masked region' );
costs(seeds) = 0;
% -------------------------------------------------------------------------
% Map the physical orientations gi to pixel orientations gi2; these pixel
% orientations are then projected onto the boundaries of a hypercube
% matching the 1x1x1x...x1 radius neighborhood of each pixel:
gi2 = projectGradients(gi,opt.ijk2xyz);
% Now, map back these pixel orientations to NON-normalized physical
% directions + distances, so that the directional costs can be normalized
% by the actual length of the jumps they represent:
dcosts = correctDirectionalCosts(dcosts,gi2,opt.ijk2xyz);
% -------------------------------------------------------------------------
% Compute the interpolation weights and interpolation neighbors:
[neighbors,weights] = computeInterpolationParameters(gi2);
% -------------------------------------------------------------------------
% Convert variable to doubles to avoid conflicts with the mex function:
costs     = double(costs);
dcosts    = double(dcosts);
mask      = double(opt.mask);
neighbors = double(neighbors);
weights   = double(weights);
% -------------------------------------------------------------------------
% Proceed with the iterations:
allvisited = false;
tcost0     = inf;
tcost      = inf;
ccost0     = inf;
ccost      = inf;
stop       = false;
for it=1:opt.maxiters
    % The basic step is carried out by the mex function:
    [costs,odirs] = fastSweepingStep( ...
        costs, ...
        dcosts, ...
        mask, ...
        neighbors-1, ...  % Note the mex file uses C-based indexing!!!
        weights     );
    % Check if all masked voxels have been visited at least once, so that
    % we can begin checking the stop criterion:
    if(~allvisited)
        allvisited = all( costs(opt.mask)>=0 );
    end
    if(allvisited)
        tcost = costs(opt.mask);
        tcost = sum( tcost(:) );
        if(~isinf(tcost0))
            ccost = abs(tcost0-tcost);
            if(~isinf(ccost0))
                % We can actually check
                if( ccost<=ccost0*(opt.change) )
                    stop = true;
                end
            else
                ccost0 = ccost;
            end
        end
        tcost0 = tcost;
    end
    % Print information if needed
    if(opt.verbose)
        fprintf(1,'#%d: [Total cost: %1.5g] [Total change in cost: %1.5g] [All visited: %u] [Stop: %u]\n',it,tcost,ccost,allvisited, stop);
    end
    % Stop if we are done:
    if(stop)
        break;
    end   
end
% -------------------------------------------------------------------------
% Transform the indexed optimal arrival directions to actual D-dimensional
% orientations taken from gi.
opt.mask = opt.mask & (costs>=0); % Unreachable points and seeds are treated in a different way
odirs = odirs(:)+1; % C to matlab indexing
dirs = zeros([prod(fov),ndim]);
for d=1:ndim
    dirs( opt.mask, d ) = gi( round(odirs(opt.mask)), d );
end
dirs( ~opt.mask, : ) = 0;
dirs = reshape(dirs,[fov,ndim]);
% -------------------------------------------------------------------------
% It only remains to fix the proper value to unreachable pixels:
costs(costs<0) = inf;
% -------------------------------------------------------------------------
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gi2 = projectGradients(gi,ijk2xyz)
[N,D] = size(gi);
% Invert the ijk2xyz matrix to get the reverse transform:
xyz2ijk = ijk2xyz\eye(D+1); % (D+1) x (D+1)
% Check where the gi map to:
gi2 = xyz2ijk*[gi';ones(1,N)]; % (D+1) x N
% Check where the origin maps to:
gi0 = xyz2ijk*[zeros(D,1);1]; % (D+1) x 1
% So that pixel displacements are (broadcasted operation):
gi2 = gi2-gi0;     % (D+1) x N
gi2 = gi2(1:D,:)'; % N x D
% Compute the maximum absolute value of the D components of each direction
% in gi2 and normalize each direction by this value. This way, we make sure
% that each row of gi2 has a component equal to 1 and the other two
% components <=1, i.e. we project gi2 in a 1x1x...x1 radius neighborhood:
M   = max(abs(gi2),[],2); % N x 1
gi2 = gi2./M;        % N x D, broadcasted operation
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dcosts = correctDirectionalCosts(dcosts,gi2,ijk2xyz)
[N,D] = size(gi2);
% Check which anatomical position each gi2 maps to:
gi  = ijk2xyz*[gi2';ones(1,N)]; % (D+1) x N
% Check where the origin of the image maps to:
gi0 = ijk2xyz*[zeros(D,1);1];   % (D+1) x 1
% So that anatomical displacements are (broadcasted operation):
gi  = gi - gi0;   % (D+1) x N
gi  = gi(1:D,:)'; % N x D
% The norm of each row of gi stands for the phsyical length each jump
% corresponds to, which is used to correct the units-less dcosts with
% length units:
ni  = sqrt(sum(gi.*gi,2)); % N x 1
dcosts = dcosts.*ni; % Broadcasted operation
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [neighbors,weights] = computeInterpolationParameters(gi2)
[N,D] = size(gi2);
% Start by creating all the available offsets in the 1x1x1x...x1 radius
% neighborhood, which has size 3^D
allneighs = generateDDoffsets(D);
TN = size(allneighs,1);
ct = (TN+1)/2; % This is the position of the center of the neighborhhod
% The interpolation is based on as many neighbors as dimensions
% has the image:
neighbors = zeros(size(gi2));
weights   = zeros(size(gi2));
% Proceed with all rows of gi2, one by one
for n=1:N
    % Pick up one direction:
    g = gi2(n,:);
    % Compute the distances from g to each of the potential neighbors to be
    % used:
    dists = (g-allneighs);       % TN x D, broadcasted
    dists = sum(dists.*dists,2); % TN x 1
    % Make sure the center point of the neighborhood is not used:
    dists(ct) = inf;
    % Sort the distances in ascending order so that the D closest neighbors
    % can be found:
    [sf,sfname,vs] = check_software_platform;
    if(sf==1)
        [~,ptr] = sort(dists,'ascend','MissingPlacement','last');
    else
        [~,ptr] = sort(dists,'ascend');
    end
    ptr = ptr(1:D);
    % Create the interpolation matrix:
    intM = allneighs(ptr,:)';
    intM = intM\eye(D);
    % Compute the interpolation weights:
    weight = intM*(g');
    weight = abs(weight');
    weight = weight/sum(weight);
    % Assign values to the outputs:
    neighbors(n,:) = ptr;
    weights(n,:)   = weight;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function allneighs = generateDDoffsets(D)
% Create all the available offsets in the 1x1x1x...x1 radius neighborhood, 
% which has size 3^D
% Basic offset [-1,0,1], which is replicated as many times as image
% dimensions we have to be used with ndgrid:
offset = cell(1,D); 
for d=1:D
    offset{d} = [-1,0,1];
end
% Define the set of D outputs returned by ndgrid
alloffsets = cell(1,D);
% Call ndgrid to produce the grid. Note that ndgrid, as opposed to
% meshgrid, does not swap the first two dimensions, so that we always get 
% the proper ordering of the neighbors as required by the mex function:
[alloffsets{:}] = ndgrid(offset{:});
% Create a 2-D array to store the result:
allneighs = zeros( numel(alloffsets{1}), D );
for d=1:D
    allneighs(:,d) = alloffsets{d}(:);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
