function [mapl,dti,lopt] = atti2mapl( atti, gi, bi, varargin )
% function [mapl,dti,lopt] = atti2mapl( atti, gi, bi, 
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Computes the coefficients of the MAP-MRI signal model as described by
%   Ozarslan et al:
%    
%        Ozarslan E, Koay CG, Shepherd TM, Komlosh ME, Irfanoglu MO,
%        Pierpaoli C, Basser PJ. "Mean apparent propagator (MAP) MRI: a 
%        novel diffusion imaging method for mapping tissue microstructure". 
%        Neuroimage. 2013 Sep; 78:16-32.
%
%   The attenuation signal is represented as a combination of Hermite
%   functions (separable for each space direction) in a transformed
%   "anatomical space" that comes from the rotation and stretching of the
%   original space according to the tensor model previously computed. To
%   fit the model, we use the MAPL (Laplacian-regularized) model described
%   by Fick et al.:
%
%        Fick RHJ, Wassermann D, Caruyer E, Deriche R. "MAPL: Tissue 
%        microstructure estimation using Laplacian-regularized MAP-MRI and 
%        its application to HCP data". Neuroimage. 2016 Jul; 134:365-385. 
%     
%
%   Mandatory inputs:
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradients over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%         These are the scattered q-space samples.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that a multi-shell acquisition is arranged. NOTE: at least 6
%         evenly spaced gradients with bi < bcut (see below) should be
%         present so that the tensor model can be reliably estimated.
%
%   The outputs of the function are as follows:
%
%      mapl: a MxNxPxK array with the K MAPL coefficients computed at each
%         voxel within the field of view. The number of coefficients, K, is
%         directly controlled by the maximum degree of the Hermite
%         polynomials used for the expansion (Nmax):
%             K = (Nmax+2)(Nmax+4)(2*Nmax+3)/24
%      dti: MxNxPx6, a double array with the estimated tensor model  at
%         each voxel. Note this is necessary to interpret the data in mapl,
%         since the coefficients refer to basis functions that are both
%         rotated and scaled dependint on the properties of this dti.
%      lopt: MxNxP, a double array with the regularization parameter used
%         to weight the energy of the Laplacian at each voxel. If the
%         optional parameter 'lambda' (see below) is a positive, finite
%         value, lopt will be filled with this value. In case GCV is used,
%         then lopt stores the optimal value used at each voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   General optional parameters:
%
%      Nmax: the maximum degree of the Hermite polynomials used in the
%         expansion, must be an even integer >=0 (default: 4).
%      lambda: this is a regularization parameter that penalizes highly
%         irregular solutions of the signal. In precise terms, the
%         quadratic deviation from the linear model is added lambda times
%         the energy of the Laplacian of the signal to compute the final 
%         cost function:
%            - If lambda is a finite, positive number, then this fixed
%              value is used for all voxels.
%            - If lambda is empty, or is NaN, or is Inf, or is negative,
%              then the appropriate value is estimated at each voxel using
%              Generalized Cross Validatio. Note this will dramatically
%              increase the processing time, since a large number of matrix
%              inversions have to be performed at each voxel.
%         (default: 0.2, i.e. NO GCV is used).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      bcut: b-values over this value are ignored to estimate the diffusion
%         tensor, which is used to transform the native space to the
%         anatomical space (default: 2000 s/mm^2).
%      ADC0: estimated diffusivity of free water at body temperature. It is
%         used to determine the lower and upper bounds of the eigenvalues
%         of the dti and perform sanity checks (default: 3.0e-3).
%      tau: the effective diffusion time of the acquisition, in seconds, 
%         necessary to determine the regularization matrix (default: 35e-3).
%      const: whether (true) or not (false) constrain the estimated MAPL
%         coefficients to fulfill two conditions:
%             1- E(q) = 1.0
%             2- P(R) >= 0 at a pre-defined set of values of R
%         With these constraints, the estimation has a rigorous physical 
%         meaning but a quadratic problem must be solved at each voxel. 
%         Without it, we just need to solve an unconstrained least squares 
%         problem, which will be pretty faster (default: true, i.e. use 
%         constraints).
%      constRad: if const is set true, inequality constraints of the form
%         P(x_i,y_i,z_i) >= 0, i=1..CI are enforced in the quadratic
%         program. This option controls the value of CI. In brief: two 3-D
%         spherical grids are arranged together in an interleaved scheme so
%         that each one covers the normalized R-space in the long term
%         (R<=10.0) and short term (R<=5.0). The number of constraints
%         equals: CI = 2*45*constRad+1 (default: 7, hence CI=631). Note the
%         number of inequality constraints dramatically impacts the
%         computational load (hence the execution time) of the algorithm.
%
%   Other general options:
%
%      maxthreads: the algorithm is run as a multi-threaded mex. This is 
%         the maximum allowed number of threads, which can indeed be 
%         reduced if it exceeds the number of logical cores (default: the 
%         number of logical cores in the machine).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).

% Check the mandatory input argments:
if(nargin<3)
    error('At lest the atti volume, the gradient table, and the b-values must be supplied');
end
[M,N,P,G] = size(atti);
if(~ismatrix(gi)||~ismatrix(bi))
    error('gi and bi must be 2-d matlab matrixes');
end
if(size(gi,1)~=G)
    error('The number of rows in gi must match the 4-th dimension of atti');
end
if(size(gi,2)~=3)
    error('The gradients table gi must have size Gx3');
end
if(size(bi,1)~=G)
    error('The number of b-values bi must match the 4-th dimension of atti');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end
% Make sure the gi are normalized:
ni = sqrt(sum(gi.*gi,2));
gi(:,1) = gi(:,1)./ni;
gi(:,2) = gi(:,2)./ni;
gi(:,3) = gi(:,3)./ni;
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.Nmax = 4;           optchk.Nmax = [true,true];       % always 1x1 double
opt.lambda = 0.2;       optchk.lambda = [true,true];     % always 1x1 double
opt.tl = 1.0e-7;        optchk.tl = [true,true];         % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];         % always 1x1 double
opt.bcut = 2000;        optchk.bcut = [true,true];       % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.tau = 35.0e-3;      optchk.tau = [true,true];        % always 1x1 double
opt.const = true;       optchk.const = [true,true];      % always 1x1 boolean
opt.constRad = 7;       optchk.constRad = [true,true];   % always 1x1 double
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Sanity checks on Nmax
assert( opt.Nmax>=0, 'Nmax should be postive' );
Nmax = round(opt.Nmax);
assert( abs(Nmax-opt.Nmax)<10*eps, 'Nmax should be an integer' );
Nmax = floor(Nmax/2);
assert( abs(2*Nmax-opt.Nmax)<10*eps, 'Nmax should be even' );
% -------------------------------------------------------------------------
% Sanity checks on lambda
if(isempty(opt.lambda))
    opt.lambda = -1.0;
end
assert( numel(opt.lambda)==1, 'Optional parameter lambda must be a scalar' );
if( isnan(opt.lambda) || isinf(opt.lambda) )
    opt.lambda = -1.0;
end
% -------------------------------------------------------------------------
% Sanity checks and helper variables:
atti(atti<opt.tl) = opt.tl;
atti(atti>opt.tu) = opt.tu;
pdti = (bi<=opt.bcut);
Ndti = length(find(pdti));
assert(Ndti>=6,'At least 6 b-values must be less or equal than bcut to proceed with DTI estimation');
% -------------------------------------------------------------------------
% Estimate the tensor model:
dti = atti(:,:,:,pdti);  % M x N x P x Ndti
dti = atti2dti( dti, gi(pdti,:), bi(pdti), 'mask', opt.mask, ...
    'wls', true, 'nonlinear', true, 'wlsit', 3, 'wsc', 0.01, ...
    'tol', 1.0e-12, 'rcondth', 1.0e-6, 'fixmode', 'a', ...
    'maxthreads', opt.maxthreads, 'unroll', false );
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask  = reshape(opt.mask,[M*N*P,1]);
atti  = reshape(atti,[M*N*P,G]);
atti  = atti(mask,:);
dti2  = reshape(dti,[M*N*P,6]);
dti2  = dti2(mask,:);
% -------------------------------------------------------------------------
% Create constraints, if necessary:
xyz = [];
if(opt.const)
    % We will create two interleaved sphercial meshes with evenly spaced
    % raddi each:
    %   - The first one, with a maximum normalized radius 10.0, for which
    %     the tails of the Hermite functions reasonably vanish.
    %   - The second one, with a maximum normalized radius 4.0, for which
    %     most of the energy of the EAP is contained (this is based on
    %     empirical considerations).
    % In both cases, 45 gradient directions are designed (note that below a
    % normalized radius 4.0 this implies an effective sampling of 30+30=60
    % gradient directions.
    % The number of radial samples is governed by the 'constRad' optional
    % argument:
    MX1 = 10.0;
    MX2 = 4.0;
    NG = 45;
    NP = max(round(opt.constRad),1);
    % The radii are evenly spaced in both cases:
    radii1 = (1:NP)/NP*MX1;
    radii1 = ones(NG,1)*radii1;
    radii1 = radii1(:);
    radii2 = (1:NP)/NP*MX2;
    radii2 = ones(NG,1)*radii2;
    radii2 = radii2(:);
    % The gradient directions of each mesh are desgined together to achieve
    % a proper interleaving:
    [G,sets] = designGradients( [NG,NG], ...
        'weights', [1,1], ...
        'nocomb', true );
    G1 = G(sets==1,:);
    G2 = G(sets==2,:);
    G1 = repmat(G1,[NP,1]);
    G2 = repmat(G2,[NP,1]);
    % Put both meshes together and add the origin:
    x = [0;radii1.*G1(:,1);radii2.*G2(:,1)];
    y = [0;radii1.*G1(:,2);radii2.*G2(:,2)];
    z = [0;radii1.*G1(:,3);radii2.*G2(:,3)];
    % Put all together:
    xyz = [x,y,z];
end
% -------------------------------------------------------------------------
% Call the mex function:
[maplM,loptM] = atti2mapl_( double(atti'), double(dti2'), gi, bi, opt.Nmax, ...
    opt.lambda, opt.ADC0, opt.tau, xyz, opt.maxthreads );
% -------------------------------------------------------------------------
% Reshape back to output the result:
K = size(maplM,1);
% ---
mapl  = zeros(M*N*P,K);
mapl(mask,:) = maplM';
mapl = reshape(mapl,[M,N,P,K]);
% ---
lopt = zeros(M*N*P,1);
lopt(mask,:) = loptM;
lopt = reshape(lopt,[M,N,P]);
% -------------------------------------------------------------------------
end
