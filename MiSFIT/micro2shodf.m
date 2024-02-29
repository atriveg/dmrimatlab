function sh = micro2shodf( atti, gi, bi, lpar, lperp, f, varargin )
% function sh = micro2shodf( atti, gi, bi, lpar, lperp, f, ...
%                                'option1', value1, 'option2, value2, ... )
%
%   Computes the SH coefficients of the ODF that best fits the multi-shell
%   attenuation signal atti (with gradients table gi and b-values bi)
%   according to the convolutional model:
%
%       atti(u,b) = (1-f)*exp(-b*Diso) 
%              + f*Integral_{S}Phi(v)exp(-b*((lpar-lperp)(u*v)^2+lperp))dv,
%
%   where lpar and lperp parameterize the impulse response as an elemental
%   rank-2 tensor, f stands for the partial volume fraction of water
%   confined in the neural axons, and Diso is the free-water
%   (extra-cellular) diffusivity. Phi(v) is the ODF we aim at estimating.
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
%      + The next three inputs are obtained with atti2micro:
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
%         the ODF at each imaged voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      L: an even integer with the maximum order of the SH to be used
%         (default: 8).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem (default 0.001).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
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
%         and lperp as provided by atti2micro. This may become especially
%         important if the 'optimal' flag is true, since otherwise many
%         matrix inversions will become singular. If true, three
%         corrections are performed:
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
%      optimal: wether (true) or not (false) use a globally optimal
%         (slower) LLS fitting of the ODF for all shells or a subotimal
%         (faster) LLS fitting shell-by-shell (default: true).
%      chunksz: if optimal == false, the LLS problem reduces to the product
%         of the atti signal with a precomputed matrix. To improve the
%         performance, cunksz voxels are processed together Note that 
%         decrasing this value can dramatically slow down processing in 
%         parallel (default: 1000).
%      recrop: wether (true) or not (false) cropping again the signal to
%         the interval [tl,tu] after the free-water compartment has been
%         substracted (default: false).
%      bth: 1x1, b-values that differ from each other less than this
%         threshold are assumed to belong to the same shell (default: 1).

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
opt.L = 8;              optchk.L = [true,true];        % always 1x1 double
opt.lambda = 0.001;     optchk.lambda = [true,true];   % always 1x1 double
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
opt.optimal = true;     optchk.optimal = [true,true];  % always 1x1 boolean
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
% -------------------------------------------------------------------------
opt.chkmod = true;      optchk.chkmod = [true,true];   % always 1x1 boolean
opt.flperp = 0.001;     optchk.flperp = [true,true];   % always 1x1 double
opt.Flperp = 0.999;     optchk.Flperp = [true,true];   % always 1x1 double
% -------------------------------------------------------------------------
opt.chunksz = 1000;     optchk.chunksz = [true,true];  % always 1x1 double
opt.recrop = false;     optchk.recrop = [true,true];   % always 1x1 boolean
opt.bth = 1;            optchk.bth = [true,true];      % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Avoid repeated calls to these functions, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
% Make sure the dwi lay within the proper range:
atti(atti>opt.tu) = opt.tu;
atti(atti<opt.tl) = opt.tl;
% -------------------------------------------------------------------------
% Auto-detect the shells to work with and perform sanity checks:
[bs,ps,Ns] = auto_detect_shells(bi,opt.bth);
%assert(Ns>1,'At least two different b-values are required for this method');
p = cell(1,Ns);
for n=1:Ns
    p{n} = find( abs(ps-n)<0.5 );
end
% -------------------------------------------------------------------------
% Reshape things to work comfortably and mask to reduce computations:
mask  = reshape( opt.mask, [M*N*P,1] );
atti  = reshape( atti,     [M*N*P,G] ); atti  = atti(mask,:);  % QxG
lpar  = reshape( lpar,     [M*N*P,1] ); lpar  = lpar(mask,:);  % Qx1
lperp = reshape( lperp,    [M*N*P,1] ); lperp = lperp(mask,:); % Qx1
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
        atti(atti>opt.tu) = opt.tu;
        atti(atti<opt.tl) = opt.tl;
    end
end
% -------------------------------------------------------------------------
% Estimate the ODF:
shc = compute_shodf_from_micro(atti,lpar,lperp,bs,gi,p,...
    opt.L,opt.lambda,opt.optimal,opt.chunksz);
% -------------------------------------------------------------------------
% Reshape SH and otuput the result:
K              = size(shc,2)+1;
sh             = zeros(M*N*P,K);
sh(:,1)        = 1/sqrt(4*pi); % Fill backgrpound with isotropic ODFs
sh(mask,2:end) = shc; % Actually computed values
sh             = reshape(sh,[M,N,P,K]);
