function [sh,nit,mu] = micro2shsqrtodf( atti, gi, bi, lpar, lperp, f, varargin )
% function [sh,nit,mu] = micro2shsqrtodf( atti, gi, bi, lpar, lperp, f, ...
%                                'option1', value1, 'option2, value2, ... )
%
%   Computes the SH coefficients of the squared root of a unit mass, 
%   strictly non-negative ODF that best fits the multi-shell attenuation 
%   signal atti (with gradients table gi and b-values bi) according to the 
%   convolutional model:
%
%       atti(u,b) = (1-f)*exp(-b*Diso) 
%              + f*Integral_{S}Phi(v)exp(-b*((lpar-lperp)(u*v)^2+lperp))dv,
%
%   where lpar and lperp parameterize the impulse response as an elemental
%   rank-2 tensor, f stands for the partial volume fraction of water
%   confined in the neural axons, and Diso is the free-water
%   (extra-cellular) diffusivity. Phi(v) is the ODF whose squared root we
%   aim at estimating, so that Phi(v) will be non-negative at any possible
%   point of its continuous domain.
%
%   The SH coefficients of the ODF and the SH coefficients of its squared
%   root are analytically related via Wigner's symbols, and the latter are
%   estimated by solving a bi-quadratic problem with a unique quadratic
%   restriction (if the ODF has unit mass, the coefficients vector of its
%   squared root has norm 1).
%
%   The coefficients of the ODF itself can be retrieved from the
%   coefficients of its squared root (i.e. the output of this function)
%   with a call to sh2squaredsh. NOTE: if the squared root of the ODF
%   spands up to an even order L, then the ODF spans up to an order 2L.
%
%   Mandatory inputs:
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradients over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that a multi-shell acquisition is arranged.
%
%      + The next three inputs are obtained with atti2micro, and are used
%        to compute the convolution kernel:
%
%      lpar: a MxNxP double array with the parallel diffusvity modeling
%         the impulse response (should fulfill 0<lpar<=Diso).
%      lperp: a MxNxP double array with the perpendicular diffusvity
%         modeling the impulse response (should fulfill 0<lerp<lpar).
%      f: a MxNxP double array with the partial volume fraction of
%         intra-cellular water (should fulfill 0<=f<=1). If an empty array
%         is passed, then f=1 for all voxels is assumed, so that
%         ones(M,N,P) has the same effect as [].
%
%   Outputs:
%
%      sh: A MxNxPx(L+1)(L+2)/2 double array with the SH coefficients of
%         the squared root of the ODF at each imaged voxel. NOTE: upon
%         using sh2squaredsh to compute the ODF itself, an array sized
%         MxNxPx(2L+1)(2L+1)/2 will be obtained.
%      nit: MxNxP, the number of iterations it took to converge or -1 if
%         convergence failed.
%      mu: MxNxP, the value of the Lagrange multiplier if Newton-Raphson's
%         method is used.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      L: an even integer with the maximum order of the SH to be used
%         FOR THE SQUARED ROOT OF THE ODF (default: 6, meaning that the
%         ODF itself will be expanded up to order 12).
%      lambda: the Laplace-Beltrami regularization parameter, which applies
%         to the coefficients OF THE ODF itself (default 0.001).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      ADC0: estimated diffusivity of free water at body temperature 
%         (Diso). Should use the same as in atti2micro (default: 3.0e-3).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with isotropic ODFs
%         (default: all trues).
%
%   Sanity checks on the micro-structure model.
%
%      chkmod: wether (true) or not (false) perform sanity checks over lpar
%         and lperp as provided by atti2micro. If true, three corrections
%         are performed:
%            + lpar is ensured to be in the range (ADC0/20,ADC0);
%            + lperp is ensured to be greater than lpar*flperp (see below);
%            + lperp is ensured to be less than lpar*Flperp (see below).
%         (default: true)
%      flperp: if chkmod == true, this parameter provides a lower threshold
%         for lperp as a function of lpar (default: 0.001).
%      Flperp: if chkmod == true, this parameter provides an upper 
%         threshold for lperp as a function of lpar (default: 0.999).
%
%   Advanced parameters:
%
%      recrop: wether (true) or not (false) cropping again the signal to
%         the interval [tl,tu] after the free-water compartment has been
%         substracted (default: false).
%      bth: 1x1, b-values that differ from each other less than this
%         threshold are assumed to belong to the same shell (default: 1).
%      chunksz: 1x1, the block size for chunked operations (default: 1000).
%
%   The following options are directly passed to signal2sqrtshodf, which is
%   the (mex) core function that actually performs the optimization:
%
%      algorithm: One of 'L' for Levenberg-Marquardt's, 'N' for
%         Newton-Raphson's or 'C' for a combination of both (L-M's output
%         is used to initialize N-R's). With L-M's, the DC component of the
%         SH coefficients is written as sqrt(1-||psi_nonDC||^2) and a pure
%         least squares problem is solved in K-1 variables. With N-R's, the
%         unit norm constraint is imposed via a Lagrange multiplier, so
%         that K+1 unknowns are solved for (default: 'C').
%      T: the maximum number of iterations with either of the previous
%         methods (default: 100).
%      maxfails: the algorithm will be assumed to have converged if
%         maxiters iterations are performed without an improvement in the
%         cost function (default: 10).
%      thCost: with L-M's, convergence is assumed if the relative change in
%         the cost function from one iteration to the next one is below
%         this threshold (default: 0.001).
%      thGrad: with N-R's, the convergence threshold for the gradient
%         (default: 0.0001).
%      thCond: with N-R's, the convergence threshold for the unit-mass
%         constraint of the ODF (default: 0.0001);
%      rho0: the initial value of the adaptive damping factor for
%         Levenberg-Marquardt's optimization (default: 1.0).
%      minrcn: minimum reciprocal condition number before a matrix is
%         considered singular and the damping factor is increased
%         (default: 1.0e-5).
%      psi0eps: for L-M's, the minimum value allowed for psi(0,0), the DC
%         component of the squared root of the ODF (default: 1.0e-4).
%      maxthreads: ONLY IN POSIX SYSTEMS the algorithm is run with
%         multiple threads. This is the maximum allowed number of threads,
%         which can indeed be reduced if it exceeds the number of logical
%         cores (default: the number of logical cores in the machine).

% Check the mandatory input arguments:
if(nargin<6)
    error('At lest the atti volume, the gradient table, and the b-values, lpar, lperp, and f must be supplied');
end
[M,N,P,G] = size(atti);
assert(ismatrix(gi),'gi must be a 2-d matlab matrix');
assert(ismatrix(bi),'bi must be a 2-d matlab matrix');
assert(size(gi,1)==G,'The number of rows in gi must match the 4-th dimension of atti');
assert(size(gi,2)==3,'gi must be Gx3');
assert(size(bi,1)==G,'The number of entries in bi must match the 4-th dimension of atti');
assert(size(bi,2)==1,'bi must b a column vector');
assert(isequal(size(lpar),[M,N,P]),sprintf('lpar should have size [%d,%d,%d] for the atti provided',M,N,P));
assert(isequal(size(lperp),[M,N,P]),sprintf('lperp should have size [%d,%d,%d] for the atti provided',M,N,P));
if(~isempty(f))
    assert(isequal(size(f),[M,N,P]),sprintf('f should have size [%d,%d,%d] for the atti provided',M,N,P));
end
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.L = 6;              optchk.L = [true,true];          % always 1x1 double
opt.lambda = 0.001;     optchk.lambda = [true,true];     % always 1x1 double
opt.tl = 1.0e-7;        optchk.tl = [true,true];         % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];         % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
% -------------------------------------------------------------------------
opt.chkmod = true;      optchk.chkmod = [true,true];     % always 1x1 boolean
opt.flperp = 0.001;     optchk.flperp = [true,true];     % always 1x1 double
opt.Flperp = 0.999;     optchk.Flperp = [true,true];     % always 1x1 double
% -------------------------------------------------------------------------
opt.recrop = false;     optchk.recrop = [true,true];     % always 1x1 boolean
opt.bth = 1;            optchk.bth = [true,true];        % always 1x1 double
opt.chunksz = 1000;     optchk.chunksz = [true,true];    % always 1x1 double
% -------------------------------------------------------------------------
opt.algorithm = 'C';    optchk.algorithm = [true,true];  % always a char
opt.T = 100;            optchk.T = [true,true];          % always 1x1 double
opt.maxfails = 10;      optchk.maxfails = [true,true];   % always 1x1 double
opt.thCost = 0.0001;    optchk.thCost = [true,true];     % always 1x1 double
opt.thGrad = 0.0001;    optchk.thGrad = [true,true];     % always 1x1 double
opt.thCond = 0.0001;    optchk.thCond = [true,true];     % always 1x1 double
opt.rho0   = 1.0;       optchk.rho0 = [true,true];       % always 1x1 double
opt.minrcn = 1.0e-5;    optchk.minrcn = [true,true];     % always 1x1 double
opt.psi0eps = 1.0e-4;   optchk.psi0eps = [true,true];    % always 1x1 double
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true]; % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Avoid repeated calls to these functions, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
% Find a suitable initial iteration for the algorithms by using the
% non-constrained function.
% -- Start by computing the SH coefficients of the unconstrained ODF:
sh = micro2shodf( atti, gi, bi, lpar, lperp, f, ...
    'L', opt.L, ...
    'lambda', opt.lambda, ...
    'tl', opt.tl, ...
    'tu', opt.tu, ...
    'ADC0', opt.ADC0, ...
    'mask', opt.mask, ...
    'chkmod', opt.chkmod, ...
    'flperp', opt.flperp, ...
    'Flperp',opt.Flperp, ...
    'optimal', true, ...
    'chunksz', opt.chunksz, ...
    'recrop', opt.recrop, ...
    'bth', opt.bth ...
    ); % M x N x P x K, K=(L+1)(L+2)/2
% -- Use the gradients table passed, gi, to synthesize ODF values:
odf = sh2signal( sh, gi, ...
    'mask', opt.mask, 'chunksz',  opt.chunksz ); % M x N x P x G
% -- Compute the squared root of the ODF (where possible):
odf = sqrt(max(odf,0)); % M x N x P x G
% -- Fit this squared root to the basis of SH:
shc = signal2sh( odf, gi, ...
    'L', opt.L, 'lambda', opt.lambda, ...
    'chunksz',  opt.chunksz, 'mask', opt.mask );  % M x N x P x K
% -------------------------------------------------------------------------
% Make sure the ODF has mass one:
mass            = sqrt(sum(shc.*shc,4));
if(is_broadcast_available)
    shc = shc./mass;
else
    shc = bsxfun( @(x,y)(x./y), shc, mass );
end
shc(isnan(shc)) = 0;
shc(isinf(shc)) = 0;
% -------------------------------------------------------------------------
% Make sure the dwi lay within the proper range:
atti(atti>opt.tu) = opt.tu;
atti(atti<opt.tl) = opt.tl;
% -------------------------------------------------------------------------
% Auto-detect the shells to work with and perform sanity checks:
[bs,ps,Ns] = auto_detect_shells(bi,opt.bth);
p = cell(1,Ns);
for n=1:Ns
    p{n} = find( abs(ps-n)<0.5 );
end
% -------------------------------------------------------------------------
% Reshape things to work comfortably and mask to reduce computations:
K     = (opt.L+1)*(opt.L+2)/2;
mask  = reshape( opt.mask, [M*N*P,1] );
atti  = reshape( atti,     [M*N*P,G] ); atti  = atti(mask,:);  % QxG
lpar  = reshape( lpar,     [M*N*P,1] ); lpar  = lpar(mask,:);  % Qx1
lperp = reshape( lperp,    [M*N*P,1] ); lperp = lperp(mask,:); % Qx1
shc   = reshape( shc,      [M*N*P,K] ); shc   = shc(mask,:);   % QxK
if(~isempty(f))
    f = reshape( f,        [M*N*P,1] ); f     = f(mask,:);     % Qx1
end
% -------------------------------------------------------------------------
% Sanity checks on the diffusion model (if needed)
if(opt.chkmod)
    lpar(lpar>opt.ADC0)    = opt.ADC0;
    lpar(lpar<opt.ADC0/20) = opt.ADC0/20;
    lperp(lperp<lpar*opt.flperp) = opt.flperp.*lpar(lperp<lpar*opt.flperp);
    lperp(lperp>lpar*opt.Flperp) = opt.Flperp.*lpar(lperp>lpar*opt.Flperp);
end
% -------------------------------------------------------------------------
% Now, correct the signal with the free-water compartment at each shell if
% necessary:
if(~isempty(f))
    for n=1:Ns
        b = bs(n);
        % atti(:,:,:,p{n}) has size QxG_n
        % f has size Qx1
        if(is_broadcast_available)
            atti(:,p{n}) = ( atti(:,p{n}) - (1-f)*exp(-b*opt.ADC0) )./f;
        else
            atti(:,p{n}) = bsxfun( @(x,y)(x-y),  atti(:,p{n}), (1-f)*exp(-b*opt.ADC0) );
            atti(:,p{n}) = bsxfun( @(x,y)(x./y), atti(:,p{n}), f                      );
        end
    end
    % ---------------------------------------------------------------------
    % Recrop if necessary:
    if(opt.recrop)
        atti(atti>opt.tu) = opt.tu;L
        atti(atti<opt.tl) = opt.tl;
    end
end
% -------------------------------------------------------------------------
%   bs:    Ns x 1
%   ps:    G x 1
%   lpar:  Q x 1
%   lperp: Q x 1
%   L:     1 x 1
%   atti:  Q x G
%   shc:   Q x K
% Compute the convolution model:
lambda = dmri_compute_convolution_weights_ODF( bs, lpar, lperp, 2*(opt.L) ); % Q x (L+1) x Ns
% Make input sizes compatible with signal2sqrtshodf
atti   = double(atti');                  % G x Q
lambda = permute(lambda,[2,3,1]);        % (L+1) x Ns x Q
shc    = shc';                           % K x Q
Ylm    = GenerateSHMatrix(2*(opt.L),gi); % G x Kp, Kp = (2L+1)(2L+2)/2
% Options to signal2sqrtshodf:
opts.nu         = opt.lambda;
opts.algorithm  = opt.algorithm;
opts.T          = opt.T;
opts.maxfails   = opt.maxfails;
opts.thCost     = opt.thCost;
opts.thGrad     = opt.thGrad;
opts.thCond     = opt.thCond;
opts.rho0       = opt.rho0;
opts.minrcn     = opt.minrcn;
opts.psi0eps    = opt.psi0eps;
opts.maxthreads = opt.maxthreads;
% Call to signal2sqrtshodf:
[shc,nitc,muc,~,~,~] = ...
    signal2sqrtshodf( atti, lambda, ps-1, shc, Ylm, opts ); % K x Q
% Transpose back to work comfortably:
shc = shc'; % Q x K
% -------------------------------------------------------------------------
% Reshape SH and otuput the result:
sh          = zeros(M*N*P,K); % Output buffer
nit         = zeros(M*N*P,1);
mu          = zeros(M*N*P,1);
sh(:,1)     = 1;              % Fill background with isotropic ODFs
sh(mask,:)  = shc;            % Actually computed values
nit(mask,1) = nitc;
mu(mask,1)  = muc;
sh          = reshape(sh,[M,N,P,K]);
nit         = reshape(nit,[M,N,P]);
mu          = reshape(mu,[M,N,P]);
