function seeds = labelmap2seeds( labelmap, varargin )
% function seeds = labelmap2seeds( labelmap, ...
%                               'option1', value1, 'option2', value2, ... )
% 
% From a volume 'labelmap' which is assumed to contain integer values
% between 0 and 65535, randomly generates seeding points to be used with
% tractography algorithms.
%
%   MANDATORY INPUTS:
%
%      labelmap: this is a X x Y x Z array of doubles whose size should
%         match the field of view of the diffusion volumes to be used later
%         on for tractography. It is assumed to contain integers between 0
%         and 65535, so that the seeding points returned will be
%         (approximantely) placed at the locations where the labelmap
%         evaluates to 'label', where 'label' is an integer that can be
%         passed as an optional argument.
%
%   OUTPUTS:
%
%      seeds: this a 3 x nSeeds array of doubles, each column being one of
%         the seeding points created IN THE XYZ OR ANATOMICAL SPACE. The
%         number nSeeds of seeding points can be passed as an optional 
%         argument.
%      
%   OPTIONAL KEY/VALUE ARGUMENTS:
%
%      label: this is the label chosen from the the labelmap image,
%         therefore should be an integer in the range [1,65535] (as 0 is
%         reserved for background voxels). The seeding points will be
%         (approximately) placed only at those image locations where the
%         labelmap evaluates to label (default: 1).
%      nSeeds: this is the number of seeding points requested. If set to
%        -1, then the algorithm will output as many seeds as voxels of the
%        labelmap evaluate to label (default: -1).
%      jitter: this parameter introduces a random jitter in the actual
%         position of each seeding point with respect to the spatial
%         location of the corresponding labelmap point used. It is
%         measured as a fraction of the voxel size, so its value should be
%         roughly in the range [0,1] (default: 0.5).
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
%         NOTE: The ijk indices are assumed to be zero-based.

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<1)
    error('At lest the labelmap must be supplied');
end
if(ndims(labelmap)>3)
    error('Input labelmap must be a 3-D array');
end
labelmap = uint16(round(labelmap));
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.label = 1;             optchk.label = [true,true];
opt.nSeeds = -1;           optchk.nSeeds = [true,true];
opt.jitter = 0.5;          optchk.jitter = [true,true];
opt.ijk2xyz = eye(4);      optchk.ijk2xyz = [true,true];
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Make sure there is at least one seed candidate:
pp = find( labelmap == uint16(round(opt.label)) );
if(isempty(pp))
    error('Unable to find labelmap locations that evaluate to the desired label value');
end
% -------------------------------------------------------------------------
% If we have at least one seed, we can proceed
if(opt.nSeeds>=0) % Otherwise, use all candidates
    if(opt.nSeeds<length(pp))
        % Requested seeds are less than available candidates. Randomly pick
        % up opt.nSeeds of them without repetition
        pp = pp( randperm( length(pp), opt.nSeeds ) );
    else
        % Requested seeds are more than available candidates. Randomly pick
        % up opt.nSeeds of them, with or without repetition:
        pp = pp( ceil( length(pp)*rand(1,opt.nSeeds) ) );
    end
end
% -------------------------------------------------------------------------
% Retrieve ijk-locations from indices:
[ii,jj,kk] = ind2sub( size(labelmap), pp );
% Fix 0-based indexing:
ii = ii-1;
jj = jj-1;
kk = kk-1;
% Add jitter as requested:
ii = ii + 2.0*(rand(size(ii))-0.5)*(opt.jitter);
jj = jj + 2.0*(rand(size(jj))-0.5)*(opt.jitter);
kk = kk + 2.0*(rand(size(kk))-0.5)*(opt.jitter);
% Get homogeneous coordinates:
ijk = [ ii(:)'; jj(:)'; kk(:)'; ones(1,numel(ii)) ];
% Transform to the anatomical space:
seeds = (opt.ijk2xyz)*ijk;
% Erase the all-1 row:
seeds = seeds(1:3,:);
