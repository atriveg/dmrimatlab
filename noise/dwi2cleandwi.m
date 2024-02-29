function DWIfiltered = dwi2cleandwi( DWInoisy, Grads, sigma, varargin )
% function DWIfiltered = dwi2cleandwi( DWInoisy, Grads, sigma, ...
%                             Property1, value1, ...
%                                   ...
%                             Propertyn, valuen )
%
%   NOTE: this is just a re-factorization of the "jaLMMSEDWI.m" function
%   available at: 
%
%       https://es.mathworks.com/matlabcentral/fileexchange/
%          43992-joint-anisotropic-wiener-filter-for-diffusion-weighted-mri
%
%   This function implements a denoising method for DWI (Diffusion Weighted
%   Imaging) volumes contaminated with stationary Rician noise with
%   standard deviation 'sigma' in both the real and imaginary parts of the
%   original Gaussian process. This is a vector LMMSE approach, in the
%   sense that the DWI channels are filtered together using the Wiener
%   filter over the squared magnitude image. The moments estimates are
%   computed as sample means, and to avoid the over-blurring due to this
%   methodology a non-local means scheme is used. The algorihtm is
%   described in detail in the following two references:
%
%   [1] Antonio Tristan-Vega and Santiago Aja-Fernandez, 'DWI filtering
%       using joint information for DTI and HARDI', Medical Image Analysis,
%       Volume 14, Issue 2, Pages 205-218. 2010;
%
%   [2] Antonio Tristan-Vega, Veronique Brion, Gonzalo Vegas-Sanchez-
%       Ferrero, and Santiago Aja-Fernandez, 'Merging squared-magnitude
%       approaches to DWI denoising: An adaptive Wiener filter tuned to the
%       anatomical contents of the image', In Proceedings of IEEE EMBC
%       2013, pp. 507-510. Osaka (Japan). 2013,
%
%   which we ask you to cite in case you use this tool for your research.
%
% ************** ABOUT THE EXECUTION TIME OF THIS SOFTWARE
%
% NOTE 1: Note this software is intended to work with very large data sets,
%   i.e. large 4-D arrays, and it runs a rather complicated non-linear
%   algorithm at each 3-D location. This means the code is potentially
%   very slow, despite we have made an effort to vectorize the operations
%   as much as possible. To have a more clear idea on the execution times
%   you may expect in your computer, you may run the script test_jaLMMSE
%   distributed with this m-file. Note the execution times strongly depend
%   on the parameters you feed to the function as explained below.
%
% NOTE 2: This Matlab implementation does not feature all the possible
%   optimizations and computational enhancements described in the
%   papers cited above. Hence, this program may take several tens of
%   minutes to compute the filtered volume in a typical scenario with real 
%   DWI data. If computational performance is an issue for you, we strongly
%   encourage you to try the ITK/C++ open-source implementation available
%   at the following link (both the source code and multi-platform
%   precompiled binaries may be downloaded):
%
%      http://www.nitrc.org/projects/jalmmse_dwi
%
% NOTE 3: However, since the algorithm is massively parallel, you can take
%   advantage of the Parallel Computing Toolbox in case you own a proper
%   license for it. WE DO NOT USE any parallel computing capabilities in
%   this version, but you may easily modify this m-file to use them in case
%   you have this toolbox installed. Just changing the 'for y=1:Y' loop in
%   this function with a 'parfor (y=1:Y, NUMBEROFMATLABWORKERS )' should do
%   the trick, since the operations at each spatial location are completely
%   independent (we have not tested this point).
%
% **************
%
% MANDATORY INPUT ARGUMENTS:
%
%   DWInoisy: The 4-D matrix with the DWI data set, i.e. a MxNxPx(G+B)
%     array for which the three first dimensions comprise the spatial
%     locations, and the fourth dimension represents each gradient
%     direction or baseline of the data set.
%   Grads: The gradients table, i.e. a (G+B)x3 matrix with the gradient
%     directions corresponding to each channel; if the n-th channel is a
%     gradient image, then G(n,:) is a unit vector with the direction used
%     to acquire it. If the n-th channel is a baseline, G(n,:) should be
%     all zeros.
%   sigma: This is the standard deviation of noise in the complex Gaussian
%     domain, that is assumed to be constant for all spatial locations. You
%     should estimate this parameter by using dwi2noisetd.m.
%
% MANDATORY OUTPUT ARGUMENTS:
%
%   DWIfiltered: The denoised DWI data set, the same size as DWInoisy, and
%     corresponding to the same gradient table.
%
% OPTIONAL PROPERTY/VALUE INPUT ARGUMENTS:
%
%   'rs': This is the 3-D radius of the neighborhood used to estimate
%     sample moments (should be a 3x1 vector or a scalar). In case large
%     radii are used, over-blurring may arise. The non-local means scheme
%     used should keep this artifact not an issue, but the execution speed
%     would drastically worsen as 'rs' grows (Default: [2;2;2]).
%   'rc': The 3-D radius of the patches used to compute nonlocal-means
%     weights. Increasing this value will worsen the adaptiveness of the
%     algorithm, though in this case the execution speed is not compromised
%     (see [2]) (Default: [1;1;1]).
%   'beta': This scalar is also related to the nonlocal means computation
%     of the sample moments. Values far above 1.0 will produce a strong
%     over-blurring, while values much smaller than 1.0 will result in the
%     noise being not removed at all. (Default: 2.0).
%   'Ng': In the normal case, all the gradients are filtered together (and
%     all the baselines are filtered together on their own), which
%     corresponds to Ng=0. For Ng>0, only the Ng closest gradients to the
%     gradient being filtered are mixed together. This limitation may be
%     useful in case strong eddy currents are present in your data set,
%     since very different gradient directions will suffer very different
%     spatial distortions. However, Ng>0 will imply an additional loop to
%     solve the vector Wiener equation for each gradient direction,
%     dramatically increasing the computational time. (Default: 0).
%   'mask': In case it is used, this is a MxNxP binary image used to
%     tell the algorithm to skip background (0-valued in the mask) pixels.
%     In case you can provide a proper mask, the execution speed may be
%     drastically improved. You can use the dwi2otsuthreshold.m function
%     to compute this mask (Default: []).
%   'onlyUNLM': This boolean flag is used to tell the algorithm to skip the
%     Wiener filter correction (in case it is true), so that a pure
%     nonlocal means filter is applied. This will notably accelerate the
%     computations, but the filtering performance will be notably worse,
%     see reference [2]. (Default: false).
%   'filterOutliers': in case pixels with a large signal variability are
%     present, we may either filter them with the UNLM approach (if this
%     boolean flag is true) or keep them untouched (if false). Since this
%     situation will mostly appear in the CSF or background, the impact of
%     this flag in the final outcomes is not truly relevant. (Default:
%     true).

if( nargin<3 )
    error('At least the DWI input, the gradient table, and sigma must be provided');
end

% Check the input size:
[Y,X,Z,N] = size(DWInoisy);
% Parse the input arguments to find the options to use:
opt       = parse_inputs( DWInoisy, Grads, sigma, varargin );
% opt.rs: radii of the estimates window (3x1)
% opt.rc: radii of the patches used to compare pixels (3x1)
% opt.beta: the beta parameter of filter strength (1x1)
% opt.Ng: the number of gradient directions mixed together (1x1)
% opt.mask: either empty or YxXxZ binary image
% opt.gradients: Nx1 binary pointer. Tells which channels are gradients
% opt.baselines: Nx1 binary pointer. Tells which channels are baselines
% opt.Grads: Nx3 matrix; each row is either a 0 vector (baselines) or a
%       vector with norm 1 (gradients).
% opt.GradNeighbors: Either empty (if Ng=0) or a N'xNg matrix of integers,
%       where N' is the number of non-zeros in opt.gradients. The i-th row
%       contains the indices of the Ng closest gradients to the i-th
%       gradient (NOTE this reckoning does not account for baselines).

%Initiallize the output:
DWIfiltered = zeros(Y,X,Z,N);

% Separate gradient images from baseline images:
GRADIENTS   = DWInoisy( :, :, :, opt.gradients );

% Compute the RGB channels:
wR = abs(opt.Grads(opt.gradients,1)); wR = wR./sum(wR); % R weights
wG = abs(opt.Grads(opt.gradients,2)); wG = wG./sum(wG); % G weights
wB = abs(opt.Grads(opt.gradients,3)); wB = wB./sum(wB); % B weights

NWR   = sum(wR.*wR);
NWG   = sum(wG.*wG);
NWB   = sum(wB.*wB);
NWRGB = (NWR+NWG+NWB);

RC = zeros(Y,X,Z);
GC = zeros(Y,X,Z);
BC = zeros(Y,X,Z);

for n=1:length(wR)
    RC = RC + wR(n).*GRADIENTS(:,:,:,n);
    GC = GC + wG(n).*GRADIENTS(:,:,:,n);
    BC = BC + wB(n).*GRADIENTS(:,:,:,n);
end

% Compute the mean values and spatial derivatives of each RGB channel, i.e.
% the local features describing the structure of the image:
[muR,GxR,GyR,GzR,factorsR,hcorrR] = ComputeLocalFeatures3D( RC, opt.rc ); hcorrR = hcorrR*opt.beta; hcorrR = hcorrR*hcorrR;
[muG,GxG,GyG,GzG,factorsG,hcorrG] = ComputeLocalFeatures3D( GC, opt.rc ); hcorrG = hcorrG*opt.beta; hcorrG = hcorrG*hcorrG;
[muB,GxB,GyB,GzB,factorsB,hcorrB] = ComputeLocalFeatures3D( BC, opt.rc ); hcorrB = hcorrB*opt.beta; hcorrB = hcorrB*hcorrB;

% MAIN LOOP 1: Compute the second and fourth order local moments and use
% the Wiener solution to filter the pixel:
for y=1:Y
    for x=1:X
        for z=1:Z
            % Check if we need to perform the computations depending on the
            % mask:
            DOIT = true;
            if( ~isempty(opt.mask) )
                if( ~(opt.mask(y,x,z)) )
                    DOIT = false;
                end
            end
            if( DOIT )
                % We are filtering the pixel (x,y,z). First, create a
                % neighborhood around this pixel checking for out-of-bound
                % indices:
                mx = max( x-opt.rs(1), 1 );
                MX = min( x+opt.rs(1), X );
                my = max( y-opt.rs(2), 1 );
                MY = min( y+opt.rs(2), Y );
                mz = max( z-opt.rs(3), 1 );
                MZ = min( z+opt.rs(3), Z );
                % Keep the center values:
                mu0R = muR(y,x,z); mu0G = muG(y,x,z); mu0B = muB(y,x,z);
                gx0R = GxR(y,x,z); gx0G = GxG(y,x,z); gx0B = GxB(y,x,z);
                gy0R = GyR(y,x,z); gy0G = GyG(y,x,z); gy0B = GyB(y,x,z);
                gz0R = GzR(y,x,z); gz0G = GzG(y,x,z); gz0B = GzB(y,x,z);
                % Get the mean values and gradients of the pixels in the whole
                % search neighborhood:
                muiR = muR(my:MY,mx:MX,mz:MZ); muiG = muG(my:MY,mx:MX,mz:MZ); muiB = muB(my:MY,mx:MX,mz:MZ);
                gxiR = GxR(my:MY,mx:MX,mz:MZ); gxiG = GxG(my:MY,mx:MX,mz:MZ); gxiB = GxB(my:MY,mx:MX,mz:MZ);
                gyiR = GyR(my:MY,mx:MX,mz:MZ); gyiG = GyG(my:MY,mx:MX,mz:MZ); gyiB = GyB(my:MY,mx:MX,mz:MZ);
                gziR = GzR(my:MY,mx:MX,mz:MZ); gziG = GzG(my:MY,mx:MX,mz:MZ); gziB = GzB(my:MY,mx:MX,mz:MZ);
                % Compute the distances from the center to each neighbor:
                distsR = (muiR-mu0R).*(muiR-mu0R) + ...
                    (gxiR-gx0R).*(gxiR-gx0R)*factorsR(1) + ...
                    (gyiR-gy0R).*(gyiR-gy0R)*factorsR(2) + ...
                    (gziR-gz0R).*(gziR-gz0R)*factorsR(3);
                distsG = (muiG-mu0G).*(muiG-mu0G) + ...
                    (gxiG-gx0G).*(gxiG-gx0G)*factorsG(1) + ...
                    (gyiG-gy0G).*(gyiG-gy0G)*factorsG(2) + ...
                    (gziG-gz0G).*(gziG-gz0G)*factorsG(3);
                distsB = (muiB-mu0B).*(muiB-mu0B) + ...
                    (gxiB-gx0B).*(gxiB-gx0B)*factorsB(1) + ...
                    (gyiB-gy0B).*(gyiB-gy0B)*factorsB(2) + ...
                    (gziB-gz0B).*(gziB-gz0B)*factorsB(3);
                % Normalize the distances:
                dists = ( distsR./hcorrR + distsG./hcorrG + distsB./hcorrB ) ./ (NWRGB*sigma*sigma);
                % Compute the weights:
                wis   = exp(-dists);
                % Avoid over-weighting of the central pixel:
                wis(wis>0.367879441171442) = 0.367879441171442;
                % Normalize the weights:
                NORM  = sum(wis(:));
                wis   = wis./NORM;
                % Actually compute the moments:
                vals  = DWInoisy(my:MY,mx:MX,mz:MZ,:);
                vals2 = vals.*vals;
                vals4 = vals2.*vals2;
                wis   = repmat( wis, [1,1,1,N] );
                M2    = sum( sum( sum( vals2.*wis, 1 ), 2 ), 3 );
                M2    = M2(:);
                M4    = sum( sum( sum( vals4.*wis, 1 ), 2 ), 3 );
                M4    = M4(:);
                % Now, compute the moments of the original A from the
                % moments of the Rician-corrupted signal M:
                vals  = DWInoisy(y,x,z,:);
                vals  = vals(:);
                vals2 = vals.*vals;
                DIFF  = vals2 - M2;
                M2    = max( M2-2*sigma*sigma,                    1e-10 );
                M4    = max( M4-8*sigma*sigma.*(M2+sigma*sigma) , 1e-10 );
                if( opt.onlyUNLM ) % In case we want the simple UNLM, we are done
                    DWIfiltered(y,x,z,:) = sqrt(M2);
                else
                    % In this case, we need to actually use the Wiener
                    % scheme to filter the pixel. The exact model used for
                    % the prior distribution is described in [1]. Start by
                    % computing the normalization factor accounting for the
                    % overall signal variability:
                    osv = M4 - M2.*M2;
                    osv = ( osv +  1e-10 )./( M2.*M2 + 1e-10 );
                    osv = mean(osv>0);
                    % If the variability is too small (i.e. the region is
                    % likely to be the background or a very homgeneous
                    % region), we consider there is no signal variation, so
                    % that we keep the UNLM solution. In case all elements
                    % of osv are negative, this is also likely to be the
                    % background, so that we should keep the UNLM solution
                    % also (in this case osv will be NaN).
                    if( isnan(osv) || osv<1e-12 )
                        DWIfiltered(y,x,z,:) = sqrt(M2);
                    elseif( osv>1000 )
                        % If the variability is too large, C_M2M2 is close
                        % to singular and numeric problems may arise. This
                        % outliers will appear in the CSF or the
                        % background, so the exact procedure will not be
                        % that relvant:
                        if( opt.filterOutliers )
                            DWIfiltered(y,x,z,:) = sqrt(M2);
                        else
                            DWIfiltered(y,x,z,:) = vals(:);
                        end
                    else
                        % We need to actually compute the Wiener solution.
                        % This corresponds to lines 605 to 682 in the
                        % ITK/C++ implementation.
                        % Iterations number for the approximation to the
                        % matrix inversion:
                        MAXITER   = 5;
                        % First, filter the baselines They always go
                        % altogether, no matter the value of Ng:
                        pixel     = DWInoisy(y,x,z,:);
                        pixel     = pixel(:);
                        BASELINES = VectorWiener( M2(opt.baselines), DIFF(opt.baselines), osv, sigma, MAXITER );
                        % Now, filter the gradients:
                        if( opt.Ng==0 )
                            % In this case, all the gradients are also
                            % filtered together. This is much faster since
                            % the inversion is performed only once (though
                            % we invert larger matrices) and we don't need
                            % loops. However, this might be problematic is
                            % case strong eddy currents are present in the
                            % data set.
                            GRADIENTS = VectorWiener( M2(opt.gradients), DIFF(opt.gradients), osv, sigma, MAXITER );
                        else
                            % In this case, only the Ng-closest gradients
                            % to the gradient being filtered are used, so
                            % we need to loop along these gradients:
                            vals2     = vals2(opt.gradients);
                            M2        = M2(opt.gradients);
                            DIFF      = DIFF(opt.gradients);
                            GRADIENTS = vals2;
                            for n=1:numel(GRADIENTS)
                                TEMP = VectorWiener( M2(opt.GradNeighbors(n,:)), ...
                                    DIFF(opt.GradNeighbors(n,:)), ...
                                    osv, sigma, MAXITER );
                                GRADIENTS(n,1) = TEMP(1);
                            end
                        end
                        % Put the filtered baselines and gradients in
                        % place:
                        pixel(opt.baselines) = BASELINES;
                        pixel(opt.gradients) = GRADIENTS;
                        DWIfiltered(y,x,z,:) = pixel;
                    end
                end
            end % if( DOIT )
        end % for z=1:Z
    end % for x=1:X
end % for y=1:Y


return;


%--------------------------------------------------------------------------
function filtered = VectorWiener( m2, diffs, osv, sigma, MAXITER )
% function filtered = VectorWiener( m2, diffs, osv, sigma, MAXITER )
%
%    Compute the vector filtering based on the moments and differences. For
%    the inversion in the Wiener formula, we use a series expansion if
%    possible.
%

% Compute the minimum of the second order sample moments in thsi data set.
% If this value is large enough, the power series expansion is convergent,
% otherwise it isn't:
CFACT = 5.0;
minM2 = min(m2(:));
% Pre-whitening of the input:
if( minM2 > CFACT*sigma*sigma )
    % The expansion is convergent:
    bRes = CMMInversionSeries( diffs(:), m2(:), osv, sigma, MAXITER );
else
    % We have to actually invert the matrix:
    bRes = CMMInversion( diffs(:), m2(:), osv, sigma );
end
% Compute the producto with C_A2M2:
dotproduct = (bRes')*m2(:);
% Correct the output value:
filtered   = (1.0+dotproduct*osv) .* m2(:);
% Compute the sqaure root when positive, 0 otherwise:
filtered    = sqrt( max(filtered,0) );

return;

%--------------------------------------------------------------------------
function bRes = CMMInversion( diffs, m2, osv, sigma )
% Create the matrix to invert:
SGM = sigma*sigma;
CMM = (osv.*m2)*m2' + diag( 4*SGM.*(m2+SGM) );
% The output is created by inverting this matrix:
bRes = CMM\diffs;

return;

%--------------------------------------------------------------------------
function bRes = CMMInversionSeries( diffs, m2, osv, sigma, MAXITER )
% This function implements the recurssive formula in [1], eq. (20).
% Fixed parameters:
SGM  = 4*sigma*sigma;
eta  = -1/( SGM*(SGM/osv+sum(m2)) );
ev   = 1./(SGM.*m2);
% Initial estimate:
bRes = diffs;
% Recurssive rule:
for it=1:MAXITER
    bRes = diffs - SGM*(  ev.*bRes + eta*sum(bRes) ); 
end
% Final product by \tilde{C}_{MM}^{-1}
bRes = ev.*bRes + eta*sum(bRes);

return;

%--------------------------------------------------------------------------
function opt = parse_inputs( DWInoisy, Grads, sigma, cell_array )

% Check the number of arguments to parse:
M = length( cell_array );
if( rem(M,2) ~= 0 )
    error('Wrong number of property/value arguments');
end

% Write defaults:
opt.rs             = [2;2;2];
opt.rc             = [1;1;1];
opt.beta           = 2.0;
opt.Ng             = 0;
opt.mask           = [];
opt.onlyUNLM       = false;
opt.filterOutliers = false;

% Parse arguments one by one:
m = 1;
while( m+1 <= M )
    property = cell_array{m};
    value    = cell_array{m+1};
    m        = m+2;
    switch( lower(property) )
        case 'rs'
            rs = value;
            if( ~isnumeric(rs) )
                error('rs must be a 3x1 numeric array or a numeric scalar');
            end
            if( rs ~= round(rs) )
                error('rs must be a 3x1 numeric array of integers or a numeric integer');
            end
            if( numel(rs) == 1 )
                opt.rs = [rs;rs;rs];
            elseif( numel(rs) ~= 3 )
                error('Wrong size of rs, which should be a 3x1 vector or a scalar');
            else
                opt.rs = rs(:);
            end
        case 'rc'
            rc = value;
            if( ~isnumeric(rc) )
                error('rc must be a 3x1 numeric array or a numeric scalar');
            end
            if( rc ~= round(rc) )
                error('rc must be a 3x1 numeric array of integers or a numeric integer');
            end
            if( numel(rc) == 1 )
                opt.rc = [rc;rc;rc];
            elseif( numel(rc) ~= 3 )
                error('Wrong size of rc, which should be a 3x1 vector or a scalar');
            else
                opt.rc = rc(:);
            end
        case 'beta'
            beta = value;
            if( ~isnumeric(beta) )
                error('beta must be a numeric');
            end
            if( numel(beta) == 1 )
                opt.beta = beta;
            else
                error('beta should be scalar, not a vector');
            end
        case 'ng'
            Ng = value;
            if( ~isnumeric(Ng) )
                error('Ng must be a numeric');
            end
            if( Ng ~= round(Ng) )
                error('Ng must be an integer, not a floating point');
            end
            if( numel(Ng) == 1 )
                opt.Ng = Ng;
            else
                error('Ng should be scalar, not a vector');
            end
        case 'onlyunlm'
            onlyunlm = value;
            if( ~islogical(onlyunlm) )
                error('onlyUNLM must be a boolean flag (true or false)');
            end
            if( numel(onlyunlm) == 1 )
                opt.onlyUNLM = onlyunlm;
            else
                error('onlyUNLM should be scalar logical, not a vector');
            end
        case 'filteroutliers'
            filteroutliers = value;
            if( ~islogical(filteroutliers) )
                error('filteroutliers must be a boolean flag (true or false)');
            end
            if( numel(filteroutliers) == 1 )
                opt.filterOutliers = filteroutliers;
            else
                error('filterOutliers should be scalar logical, not a vector');
            end
        case 'mask'
            mask     = value;
            if( ~isempty(mask) )
                if( ~islogical(mask) )
                    error('The mask image should be a 3-D array of logicals');
                end
                opt.mask = mask;
            end
        otherwise
            error(['Unrecognized property name: ',property]);
    end
end

% Sanity checks:
[Y,X,Z,N] = size(DWInoisy);
if( ~isempty(opt.mask) )
    [Y2,X2,Z2,N2] = size(opt.mask);
    if( N2~=1 )
        error('The mask should have size YxXxZ');
    end
    if( X~=X2 || Y~=Y2 || Z~=Z2 )
        error('The mask and the DWI volume should have the same spatial locations');
    end
end
if( size(Grads,1) ~= N )
    error('Wrong gradients table: Grads should have the same number of rows as the fourth dimension of the DWI input');
end
if( size(Grads,2) ~= 3 )
    error('Wrong gradients table: Grads should have 3 columns');
end
if( ~isnumeric(sigma) )
    error('sigma should be numeric');
end
if( numel(sigma)>1 )
    error('sigma should be scalar');
end
if( sigma<=0 )
    error('sigma should be a strictly positive floating point');
end

% Look for gradient and baseline channels:
norms = sqrt(sum(Grads.*Grads,2));
opt.gradients = (norms >  0.01);
opt.baselines = (norms <= 0.01);

% Normalize the gradients table:
Grads(opt.gradients,1) = Grads(opt.gradients,1)./norms(opt.gradients);
Grads(opt.gradients,2) = Grads(opt.gradients,2)./norms(opt.gradients);
Grads(opt.gradients,3) = Grads(opt.gradients,3)./norms(opt.gradients);
Grads(opt.baselines,:) = 0;
opt.Grads              = Grads;

% Check if we are requesting too many gradients to be filtered together:
if( opt.Ng>=sum(double(opt.gradients)) )
    opt.Ng = 0;
end

% In case Ng ~= 0, not all the gradients are filtered together. In this
% case, we need to check which gradients are required to filter the i-th
% gradient. We use a separated function for this task:
opt.GradNeighbors = [];
if( opt.Ng~=0 )
    opt.GradNeighbors = find_grad_neighbors( opt.Grads(opt.gradients,:), opt.Ng );
end

return;

%--------------------------------------------------------------------------
function neighbors = find_grad_neighbors( Grads, Ng )

% Check the total number of gradients:
N = size( Grads, 1 );

% The output will be a NxNg matrix; each row 'i' will comprise the indices
% of the Ng closest gradient directions to the 'i'-th gradient:
neighbors = zeros( N, Ng );
[gx1,gx2] = meshgrid( Grads(:,1), Grads(:,1) );
[gy1,gy2] = meshgrid( Grads(:,2), Grads(:,2) );
[gz1,gz2] = meshgrid( Grads(:,3), Grads(:,3) );
% We cannot compute the Euclidean distance between the gradient directions
% since direction u and direction -u are essentially the same (water
% diffusion is only probed for orientation, not its direction). For this
% reason, we use the opposite of the absolute value of the dot product
% between the gradient directions as a distance metric (assuming all the
% gradients are normalized to unit vectors):
DIST      = -abs( gx1.*gx2 + gy1.*gy2 + gz1.*gz2 );
for n=1:N
    [~,ptr]        = sort( DIST(n,:), 2, 'ascend' );
    neighbors(n,:) = ptr(1:Ng);
end

return;

%--------------------------------------------------------------------------
function [mu,Gx,Gy,Gz,factors,hcorr] = ComputeLocalFeatures3D( I, radii )
% Computes the local mean value and the local gradients of a 3D image.
%
%    I:       the input image
%    radii:   a 3x1 vector of integers with the size of the neighborhood used
%             to compute the local values. Gaussian windows are used
%             generated for each dimension as gausswin(2*radii(d)+1). If not
%             provided, [x=1;y=1;z=1] will be assumed
%    mu:      A 3D image, the same size as I, with local mean.
%    Gx:      A 3D image, the same size as I, with the gradient in the 'x'
%             direction (dimension 2 in matlab).
%    Gy:      A 3D image, the same size as I, with the gradient in the 'y'
%             direction (dimension 1 in matlab).
%    Gz:      A 3D image, the same size as I, with the gradient in the 'z'
%             direction (dimension 3 in matlab).
%    factors: a 3x1 vector with the factors to be applied to each gradient
%             difference to estimate patch distances.
%    hcorr:   the effective reduction in the amount of noise in the
%             distances between patches because of the fitting.

I = double(I);

% Check if the radii where provided:
if( nargin<2 )
    radii = [1;1;1];
else
    if( length(radii) ~= 3 )
        radii = ones(3,1)*radii(1);
    end
end

% Create the gaussian windows for each direction:
gx = gausswin( 2*radii(1) + 1 ); gx = gx./sum(gx);
gy = gausswin( 2*radii(2) + 1 ); gy = gy./sum(gy);
gz = gausswin( 2*radii(3) + 1 ); gz = gz./sum(gz);

% Compute the local mean:
mu = My3DConv( I, gx, gy, gz );

% Create the differentiation kernels:
gdx = (-radii(1):radii(1))';
gdx = (gdx.*gx)./sum(gdx.*gdx.*gx);
gdy = (-radii(2):radii(2))';
gdy = (gdy.*gy)./sum(gdy.*gdy.*gy);
gdz = (-radii(3):radii(3))';
gdz = (gdz.*gz)./sum(gdz.*gdz.*gz);

% Create each gradient image (the minus sign is for consistence with the
% implementation of matlab's 'gradient' function:
Gx  = -My3DConv( I, gdx, gy,  gz  );
Gy  = -My3DConv( I, gx,  gdy, gz  );
Gz  = -My3DConv( I, gx,  gy,  gdz );

% Compute the scaling factors:
factors(1) = sum( (-radii(1):radii(1)).*(-radii(1):radii(1)).*gx' );
factors(2) = sum( (-radii(2):radii(2)).*(-radii(2):radii(2)).*gy' );
factors(3) = sum( (-radii(3):radii(3)).*(-radii(3):radii(3)).*gz' );

% Compute the correction in the h factor. First, compute the 'X' matrix:
[x,y,z]    = meshgrid( -radii(1):radii(1), ...
    -radii(2):radii(2), ...
    -radii(3):radii(3) );
X          = [ ones(size(x(:))), ...
    x(:), y(:), z(:), ...
    x(:).*x(:)/2, y(:).*y(:)/2, z(:).*z(:)/2, ...
    x(:).*y(:), x(:).*z(:), y(:).*z(:) ];
[g1,g2,g3] = meshgrid( gx, gy, gz );
R          = g1(:).*g2(:).*g3(:);
hcorr    = sqrt(trace(diag(R)*X*(X'*X)^(-1)*X'));
return;

%--------------------------------------------------------------------------
function out = My3DConv( I, gx, gy, gz )
% Computes a separable 3D convolution
gx  = gx(:);
gx  = permute(gx,[2,1,3]);
gy  = gy(:);
gz  = gz(:);
gz  = permute(gz,[3,2,1]);
I   = convn( I, gx, 'same' );
I   = convn( I, gy, 'same' );
out = convn( I, gz, 'same' );
return;
