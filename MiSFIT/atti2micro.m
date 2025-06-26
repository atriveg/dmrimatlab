function [lambdapar,lambdaperp,f] = atti2micro( atti, gi, bi, varargin )
% function [lambdapar,lambdaperp,f] = atti2micro( atti, gi, bi, 
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Computes the parameters of a micro-structural model for diffusion MRI
%   acquired in multiple shells (i.e. with multiple b-values), defined by 
%   a linear spherical convolution with an ODF as follows:
%
%       atti(u,b) = (1-f)*exp(-b*Diso) 
%              + f*Integral_{S}Phi(v)exp(-b*((lpar-lperp)(u*v)^2+lperp))dv,
%
%   where lpar and lperp parameterize the impulse response as an elemental
%   rank-2 tensor, f stands for the partial volume fraction of water
%   confined in the neural axons, and Diso is the free-water
%   (extra-cellular) diffusivity. Phi(v) is the ODF, which has integral 1
%   in the unit sphere S.
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
%   The three outputs of the function define the micro-structual model:
%
%      f: a MxNxP aray in the range (0,1) with the partial volume fraction
%         of non-free (i.e. fiber-confined) water at each imaged voxel.
%      lambdapar: MxNxP, the elemental diffusivity along the parallel 
%         direction.
%      lambdaperp: MxNxP, the elemental diffusivity along the perpendicular
%         direction.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   General optional parameters:
%
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem that fits SH coeffcients to the signals to
%         compute their orientational averages (default 0.006).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      ADC0: estimated diffusivity of free water at body temperature.
%         (default: 3.0e-3).
%      usef: wether (true) or not (false) computing the free water
%         compartment within each voxel. The optimization is carried out by
%         means of a multi-level greedy search algorithm: a set of values
%         of f within the allowed range is probed and the one that best
%         fits the model is selected. For the next resolution level, the
%         allowed range is narrowed around such f and the process is
%         repeated. NOTE: this estimation dramatically increases the
%         computational load of the method as explained below (default:
%         false).
%
%   For estimating the fiber-impulse response at each voxel:
%
%      mlperp: minimum allowed value for lambdaperp in all cases, so that
%          it will not fall down to zero and the model will not become
%          singular. You can pass just zero to obtain the original
%          non-regularized model (default: 0.01e-3).
%      mu: this is a regularization term imposing a penalty to
%          close-to-singular tensor models, i.e. those voxels with
%          lambdaperp <<< lambdapar. You can pass just zero to obtain the
%          original non-regularized model (default: 5.0e-5).
%      nu: this is a regularization parameter with the exact opposite
%          meaning as the previous one, i.e. it penalizes
%          close-to-spherical tensor models, i.e. those voxels with
%          lambdaperp->lambdapar. You can pass just zero to obtain the
%          original non-regularized model (default: 0.0).
%      nmax: maximum number of iterations in the Newton-Raphson
%         optimization (default: 100).
%      dl: maximum allowed change in both lambda_par and lambda_perp before
%         the Newton Raphson iterations stop (default: 1.0e-8).
%      dC: maximum allowed change in the cost function before the
%         Newton-Raphson iterations stop (default: 1.0e-6).
%      regf: wether (true) or not (false) use the regularization parameter
%         mu while optimizing for f (it makes sense only if the 'usef' flag
%         is on). Otherwise, it will first fit f without regularizing the
%         micro-structure model, then correct the spherical means with such
%         value, then re-compute lamdapar and lambdaperp with the
%         regularization term (default: true).
%      lpar: this parameter makes sense only if:
%              - just 1 shell is acquired.
%              - 2 shells are acquired and the user asks for the estimation
%                of the free water comparment. 
%              - Option 'forcelar' is set true.
%         In all these cases, lambdapar will not be estimated, so it has to
%         be externally fixed. NOTE: the user is responsible of giving lpar
%         a physically meaningful value, since this will not be internally
%         checked. It should be less than ADC0, and an appropriate value
%         could be around 1.7e-3. This parameter may be either a scalar
%         (default: 1.7e-3) or an MxNxP volume.
%      forcelpar: 1x1 logical, this flag may be set true so that the
%         algorithm uses a fixed value for the parallel diffusivity
%         (the one set with option lpar) instead of estimating its actual
%         value (default: false).
%      nolpwarn: in case the function is called in a way that lamdapar
%         has to be fixed to a constant, it will rise a warning unless this
%         flag is set 'true' (default: false).
%
%   For the free-water compartment isolation (only if usef='true'):
%
%      fmin, fmax: the lower and upper limits where f will be looked for at
%         each image voxel, so that fmin should be 0 and fmax should be 1 
%         (or close to 1) (default: 0.01, 1).
%
%   Other general options:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      bth: 1x1, b-values that differ from each other less than this
%         threshold are assumed to belong to the same shell (default: 1).
%      chunksz: To improve the performance of certain operations, cunksz 
%         voxels are processed together. Note that decrasing this value can
%         dramatically slow down processing in parallel (default: 10000).
%      verbose: wether (true) or not (false) show additional information on
%         the progress of the algorithm (default: false).

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

% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.lambda = 0.006;     optchk.lambda = [true,true];    % always 1x1 double
opt.tl = 1.0e-7;        optchk.tl = [true,true];        % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];        % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];      % always 1x1 double
opt.usef = false;       optchk.usef = [true,true];      % always 1x1 boolean
% -------------------------------------------------------------------------
opt.mlperp = 0.01e-3;   optchk.mlperp = [true,true];    % always 1x1 double
opt.mu = 5.0e-5;        optchk.mu = [true,true];        % always 1x1 double
opt.nu = 0.0;           optchk.nu = [true,true];        % always 1x1 double
opt.nmax = 100;         optchk.nmax = [true,true];      % always 1x1 double
opt.dl = 1.0e-8;        optchk.dl = [true,true];        % always 1x1 double
opt.dC = 1.0e-6;        optchk.dC = [true,true];        % always 1x1 double
opt.regf = true;        optchk.regf = [true,true];      % always 1x1 boolean
opt.lpar = 2.0e-3;      optchk.lpar = [true,false];     % variable size double
opt.forcelpar = false;  optchk.forcelpar = [true,true]; % always 1x1 boolean
opt.nolpwarn = false;   optchk.nolpwarn = [true,true];     % always 1x1 boolean
% -------------------------------------------------------------------------
opt.fmin = 0.01;        optchk.fmin = [true,true];      % always 1x1 double
opt.fmax = 1;           optchk.fmax = [true,true];      % always 1x1 double
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];      % boolean the same size as the image field
opt.bth = 1;            optchk.bth = [true,true];       % always 1x1 double
opt.chunksz = 10000;    optchk.chunksz = [true,true];   % always 1x1 double
opt.verbose = false;    optchk.verbose = [true,true];   % always 1x1 boolean
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Sanity check on the optional lpar:
if(~isequal(size(opt.lpar),[M,N,P]))
    if(isscalar(opt.lpar))
        opt.lpar = opt.lpar*ones(M,N,P);
    else
        error('The size of optional argument lpar should match the first three dimensions of atti, or otherwise be scalar');
    end
end
% -------------------------------------------------------------------------
% Auto-detect the shells to work with and perform sanity checks:
[bs,ps,Ns] = auto_detect_shells(bi,opt.bth);
if(opt.verbose)
    fprintf(1,'%d shells detected at b={',Ns);
    for n=1:length(bs)-1
        fprintf('%1.3f, ',bs(n));
    end
    fprintf('%1.3f} s/mm^2\n',bs(end));
end
if(any(bs<opt.bth))
    warning('Found one shell with b near 0. This is probably a baseline and should not be included');
end
% -------------------------------------------------------------------------
% Sanity checks on the number of acquired shells
switch(Ns)
    case 1
        if( opt.usef )
            warning('The free water compartment cannot be estimated from 1 shell. The ''usef'' flag will be ignored');
            if(~opt.nolpwarn)
                warning('The ''lambdapar'' will not be optimized, but fixed to the optional argument ''lpar'' instead');
            end
            opt.usef = false;
        end
    case 2
        if( (opt.usef) && (~opt.forcelpar) )
            if(~opt.nolpwarn)
                warning('To estimate the free water compartment up from two shells, the ''lambdapar'' will not be optimized, but fixed to the optional argument ''lpar'' instead');
            end
        end
    otherwise
end
% -------------------------------------------------------------------------
p = cell(1,Ns);
for n=1:Ns
    p{n} = find( abs(ps-n)<0.5 );
end
% -------------------------------------------------------------------------
% Unroll the DWI image and the mask for easy processing:
atti = reshape(atti,[M*N*P,G]);       % NTxG
mask = reshape(opt.mask,[M*N*P,1]);   % NTx1
atti = atti(mask,:);                  % NMxG, only masked voxels
% -------------------------------------------------------------------------
% Compute the per-shell spherical means of the signal
Smean = compute_spherical_means( atti, gi, Ns, ps, 6, opt.lambda ); % NMxNs
% -------------------------------------------------------------------------
% Remove outliers to keep values in range:
Smean(Smean<opt.tl) = opt.tl;
Smean(Smean>opt.tu) = opt.tu;
% -------------------------------------------------------------------------
% Proceed with the estimation. This will depend on the number of shells and
% the need (or not) to estimate the free water compartment:
NM = size(Smean,1);
if(~opt.usef)
    fm = ones(NM,1);
    if( (Ns==1) || (opt.forcelpar) )
        lpa = reshape(opt.lpar,[M*N*P,1]);
        lpa = lpa(mask);
        [~,lpp,~]  = compute_lambda_perp_from_E( Smean, bs, ...
            opt.mu, opt.nu, ...
            opt.mlperp, lpa, opt.nmax, opt.dl, opt.dC, opt.chunksz );
    else
        [lpa,lpp,~]  = compute_lambdas_from_E( Smean, bs, ...
            opt.mu, opt.nu, ...
            opt.mlperp, opt.nmax, opt.dl, opt.dC, opt.ADC0, opt.chunksz );
    end
else
    [lpa,lpp,fm] = atti2microO1( Smean, bs, opt );
end
% -------------------------------------------------------------------------
% Initialize the outputs:
lambdapar  = zeros(M*N*P,1);
lambdaperp = zeros(M*N*P,1);
f          = zeros(M*N*P,1);
% Fill out masked values:
lambdapar(mask)  = lpa;
lambdaperp(mask) = lpp;
f(mask)          = fm;
% Reshape back the outputs to their proper sizes:
lambdapar  = reshape( lambdapar,  [M,N,P] );
lambdaperp = reshape( lambdaperp, [M,N,P] );
f          = reshape( f,          [M,N,P] );

%%% -----------------------------------------------------------------------
function Smean = compute_spherical_means(atti,gi,Ns,ps,L,lambda)
% atti:   NMxG, the "unrolled" attenuation signal, each column a gradient
% gi:     Gx3, the gradients table
% Ns:     1x1, the number of shells
% ps:     Gx1, tells which shell each gradient belongs to
% L:      1x1, the maximum order of Spherical Harmonics
% lambda: 1x1, Laplace-Beltrami penalti term
% Smean:  NMxNs, the spherical mean for each shell

% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
Smean = zeros(size(atti,1),Ns);
for n=1:Ns % For each shell...
    ptn = abs(ps-n)<1/2;                % Gradients in this shell
    gin = gi(ptn,:);                    % Gnx3
    B   = GenerateSHMatrix( L, gin );   % GnxK, where K=(L+1)(L+2)/2
    LR  = GenerateSHEigMatrix( L );     % KxK
    WLS = (B'*B+lambda*LR^2)\(B');      % (KxK)^(-1) * (KxGn) -> KxGn
    WLS = WLS(1,:)/sqrt(4*pi);          % 1xGn, DC component
    if(is_broadcast_available)
        Smean(:,n) = sum(atti(:,ptn).*WLS,2); % NMx1
    else
        Smean(:,n) = sum( ...
            bsxfun( @(x,y)(x.*y), atti(:,ptn), WLS ), ... % atti(:,ptn).*WLS
            2); % NMx1
    end
end

%%% -----------------------------------------------------------------------
function [lpa,lpp,f] = atti2microO1( Smean, bs, opt )
% -------------
% Smean: NM x Ns, the means for each shell
% atti:  NM x Gn, attenuation signal E(q)
% gi:    Gn x 3, gradients table
% bi:    Gn x 1, bvalues
% bs:    Ns x 1, the b-values at each shell
% p:     1 x Ns cell, a pointer the gradients at each shell
% opt:   optional arguments
% -------------
% lpa:   NM x 1, parallel diffusivities (>=0)
% lpp:   NM x 1, perpendicular diffusivities (>=0)
% fm:    NM x 2, partial volume fraction of non-free (confined) water
% -------------------------------------------------------------------------
% Avoid repeated calls to library utilities:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
R  = size(Smean,1);
% -------------------------------------------------------------------------
% Initalize the partial volume fractions to twice the minimum of its
% physically meaningful range and make sure it is still in bounds:
S0 = exp(-(bs')*opt.ADC0);     % 1 x Ns
if(is_broadcast_available)
    fmin1 = Smean./S0;         % R x Ns
    fmin2 = (1-Smean)./(1-S0); % R x Ns
else
    fmin1 = bsxfun( @(x,y)(x./y), Smean,   S0   ); % R x Ns
    fmin2 = bsxfun( @(x,y)(x./y), 1-Smean, 1-S0 ); % R x Ns
end
fmin = max(1-fmin1,1-fmin2);            % R x Ns
fmin = max( max(fmin,[],2), opt.fmin ); % R x 1
fmax = opt.fmax*ones(R,1);              % R x 1
f    = min(2*fmin,fmax);                % R x 1
% -------------------------------------------------------------------------
% Call the optimization routine for the variable parameters. This will
% depend on the number of available shells:
if( (length(bs)<3) || (opt.forcelpar) )
    % Only two shells are available, so that only lpar and f can be
    % actually optimized; lpar will be fixed throughout:
    lpa = opt.lpar(:);
    lpa = lpa(opt.mask(:));
    [lpa,lpp,f,~] = compute_lambda_perp_and_f_from_E( Smean, bs, ...
        fmin, fmax, f, opt.mu, opt.nu, opt.mlperp, lpa, ...
        opt.nmax, opt.dl, opt.dC, opt.ADC0, opt.chunksz );
else
    % Three or more shells are available, so that the three parameters
    % lpar, lperp, and f can be actually optimized. Run the algorithm with
    % the previous initialization and wituout regularization or
    % restrictions on mlperp (seems to work better):
    if(opt.regf)
        [lpa,lpp,f] = compute_lambdas_and_f_from_E( Smean, bs, fmin, fmax, f,  ...
            opt.mu, opt.nu, opt.mlperp, ...
            opt.nmax, opt.dl, opt.dC, opt.ADC0, opt.chunksz );
    else
        [~,~,f] = compute_lambdas_and_f_from_E( Smean, bs, fmin, fmax, f,  ...
            0, 0, 0, opt.nmax, opt.dl, opt.dC, opt.ADC0, opt.chunksz );
    end
end
% -------------------------------------------------------------------------
if( (~opt.regf) && (length(bs)>2) )
    % Correct the spherical means signal with the value obtained for f:
    Smean = Smean - (1-f)*S0; % R x Ns
    if(is_broadcast_available)
        Smean = Smean./f;                     % R x Ns
    else
        Smean = bsxfun(@(x,y)(x./y),Smean,f); % R x Ns
    end
    % Sanity check (shouldn't be needed):
    Smean(Smean<opt.tl) = opt.tl;
    Smean(Smean>opt.tu) = opt.tu;
    % Recompute the model:
    [lpa,lpp,~]  = compute_lambdas_from_E( Smean, bs, opt.mu, opt.nu, ...
        opt.mlperp, opt.nmax, opt.dl, opt.dC, opt.ADC0, opt.chunksz );
end
% -------------------------------------------------------------------------
