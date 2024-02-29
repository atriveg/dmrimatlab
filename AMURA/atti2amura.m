function mu = atti2amura( atti, gi, bi, varargin )
% function mu = atti2amura( atti, gi, bi, 'option1', value1, ... )
%
%   Computes (apparent) moments (either full, axial, or planar) of 
%   arbitrary orders over the attenuation signal E(q) according to a 
%   mono-exponential model for single-shell acquisitions:
%
%      atti(u,b) = exp(-b*D_0(u)),
%
%   where D_0(u) is the Apparent Diffusion Coefficient measured at b=b0 for
%   each direction 'u' within the unit sphere. In precise terms:
%
%      mu_{full}^{nu}       = Integral_{R^3} ||q||^{nu} E(q) dq^3
%      mu_{axial}^{nu}(u0)  = Integral_{R}   t^{nu} E(t*u0) dt
%      mu_{planar}^{nu}(u0) = Integral_{v \perp u0} ||v||^{nu} E(v) dv^2
%
%   where nu>-3 for full moments, nu>-1 for axial moments, and nu>-2 for
%   planar moments. In all cases, nu can be any real number.
%
%   Most of the common diffusion measurements are described with this
%   scheme: RTOP = mu_{full}^0, RTPP = mu_{axial}^0(u_max), RTAP = 
%   mu_{planar}^0(u_max), QMSD = mu_{full}^2...
%
%   Inputs:
%
%      atti: a MxNxPxG double array with the attenuation signal measured at
%         each voxel within the MxNxP field of view.
%      gi: a Gx3 double array with the gradients table, i.e. each row of gi
%         should be a unit vector describing the corresponding gradient
%         direction.
%      bi: a Gx1 double array with the corresponding b-values of each entry
%         in gi. NOTE: for the AMURA model to make sense, all bi should be
%         similar, otherwise the sampling scheme is a multi-shell and a
%         differente model (such as MiSFIT) should be used. Alternatively,
%         bi can be a single scalar describing the acquired shell.
%
%   Outputs:
%
%      mu: MxNxP double array with the moment requested.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      type: 1x1 string array, one of 'f', 'a', or 'p', corresponding to
%          either full, axial, or planar. Alternatively, you can specify 
%          either of 'rtop', 'rtpp', 'rtap', or 'qmsd' (default: 'rtop').
%      nu: 1x1 double (or []) with the order of the moment to be computed.
%         If you specified a particular measure in 'type' (e.g. 'rtop'),
%         this is not used and may be left empty (default: not used).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).
%      u0: only used for axial and planar moments; MxNxPx3 double array
%         with the directions u0 for which axial and planar moments are
%         computed. If left empty, [], the direction of maximum diffusivity
%         will be internally computed and used (default: []).
%
%   Parameters related to SH computations:
%
%      L: an even integer with the maximum order of the SH to be used
%         (default: 8).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem (default 0.001).
%
%   Sanity checks on the attenuation signal:
%
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%
%   Advanced parameters:
%
%      tau: 1x1, the effective diffusion time of the dMRI sequence in
%         miliseconds (default: 70.0e-3).
%      chunksz: the evaluation of SH at desired directions is done by 
%         repeatedly calling GenerateSHMatrix. This is done chunk-by-chunk
%         for efficiency (default: 256).
%      clean: 1x1 double in the range [0,100]. This is a simple outlier 
%         rejection parameter to avoid out-of-range values: 0 means no
%         outlier rejection is applied; >0 means outlier rejection is
%         applied, the closer to 100 the more agressive (default: 50).

% Check the mandatory input arguments:
if(nargin<3)
    error('At least atti, gi, and bi must be supplied');
end
[M,N,P,G] = size(atti);
assert(ismatrix(gi),'gi must be a Gx3 matrix');
assert(size(gi,2)==3,'The second dimension of gi must have size 3');
assert(size(gi,1)==G,'The first dimension of gi must match the fourth dimension of atti');
if(isscalar(bi))
    bi = repmat(bi,[G,1]);
else
    assert(ismatrix(bi),'bi must be either Gx1 or 1x1');
    assert(size(bi,2)==1,'bi should be a column vector');
    assert(size(bi,1)==G,'The first dimension of bi should be either singleton or match the fourth dimension of atti');
end
% -------------------------------------------------------------------------
opt.type = 'rtop';      optchk.type = [true,false];    % Variable length string
opt.nu = 0;             optchk.nu = [true,false];      % Double, either 1x1 or empty
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
opt.u0 = [];            optchk.u0 = [true,false];      % Double, either MxNxPx3 or empty
% -------------------------------------------------------------------------
opt.L = 8;              optchk.L = [true,true];        % always 1x1 double
opt.lambda = 0.001;     optchk.lambda = [true,true];   % always 1x1 double
% -------------------------------------------------------------------------
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
% -------------------------------------------------------------------------
opt.tau = 70.0e-3;      optchk.tau = [true,true];      % always 1x1 double
opt.chunksz = 256;      optchk.chunksz = [true,true];  % always 1x1 double
opt.clean = 50;         optchk.clean = [true,true];    % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
assert(abs(opt.L-floor(opt.L))<100*eps,'Optional argument ''L'' must be integer');
assert(abs(opt.L/2-floor(opt.L/2))<100*eps,'Optional argument ''L'' must be even');
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
% Check consistency of outlier rejection:
assert( opt.clean>=0 && opt.clean<=100, '''clean'' must be in the range [0,100]' );
% -------------------------------------------------------------------------
% Determine the kind and the order of the moment to compute:
[type,nu] = micro2moments_type(opt.type,opt.nu);
% -------------------------------------------------------------------------
% If the moment is axial or planar, u0 must be checked, too
if( type>0 )
    if(~isempty(opt.u0))
        [M2,N2,P2,K2] = size(opt.u0);
        assert(isequal([M,N,P],[M2,N2,P2]),'The first three dimensions of u0 should match those of atti');
        assert(K2==3,'u0 must have size MxNxPx3, so that the last dimension stands for a 3x1 unit vector');
    else
        sh = atti2shadc( atti, gi, bi, 'L', 2, 'lambda', 0.001, ...
            'tl', opt.tl, 'tu', opt.tu, 'mask', opt.mask, ...
            'chunksz', opt.chunksz );                 % MxNxPx6
        tens   = shadc2dti(sh,'mask',opt.mask,...
            'chunksz',opt.chunksz,'unroll',false );   % MxNxPx6
        opt.u0 = dti2xyz( tens, 'mask', opt.mask );   % MxNxPx3
    end
end
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably:
mask  = reshape( opt.mask, [M*N*P,1] ); % M*N*P x 1
atti  = reshape( atti, [M*N*P,G] );     % M*N*P x G
atti  = atti(mask,:);                   % Q x G
if(type>0)
    u0 = reshape( opt.u0, [M*N*P,3] );  % M*N*P x 3
    u0 = u0(mask,:);                    % Q x 3
end
% -------------------------------------------------------------------------
% Compute the ADC from the attenuation signal:
atti(atti<opt.tl) = opt.tl;
atti(atti>opt.tu) = opt.tu;
atti = -log(atti); % Q x G
if(is_broadcast_available)
    atti = atti./(bi');                       % Q x G
else
    atti = bsxfun( @(x,y)(x./y), atti, bi' ); % Q x G
end
% Now atti is D_0!
% -------------------------------------------------------------------------
% SH decompositions:
B  = GenerateSHMatrix( opt.L, gi );    % GxK, where K=(L+1)(L+2)/2
LR = GenerateSHEigMatrix( opt.L );     % KxK
B  = (B'*B+(opt.lambda).*LR.*LR)\(B'); % (KxK)^(-1) * (KxG) -> KxG
% -------------------------------------------------------------------------
% Effectively compute the moments
if(type>0) % Axial or planar
    % ---------------------------------------------------------
    Q   = size(atti,1);
    mu0 = zeros(Q,1);
    % ---------------------------------------------------------
    if(type==1) % Axial
        % ------------------------------
        B   = B';                       % G x K
        % ------------------------------
        atti = atti.^(-(nu+1)/2);
        % ------------------------------
        cnst = gamma((nu+1)/2);
        cnst = cnst/(4*pi*pi*opt.tau)^((nu+1)/2);
        % ------------------------------
    else        % Planar
        % ------------------------------
        FRT = GenerateFRTMatrix(opt.L); % K x K
        B   = FRT*B;                    % K x G
        B   = B';                       % G x K
        % ------------------------------
        atti = atti.^(-(nu+2)/2);
        % ------------------------------
        cnst = gamma((nu+2)/2);
        cnst = cnst/(4*pi*pi*opt.tau)^((nu+2)/2)/2;
        % ------------------------------
    end
    % ---------------------------------------------------------
    % Work chunk-by-chunk to avoid memory blowing up
    for ck=1:ceil(Q/opt.chunksz)
        % ---------
        idi = (ck-1)*opt.chunksz+1;  % From...
        idf = min(ck*opt.chunksz,Q); % ... to
        % ---------
        % Compute atti SH for this chunk:
        sh  = atti(idi:idf,:)*B; % ( Qc x G ) x ( G x K ) -> Qc x K
        % ---------
        % SH matrix describing the present chunk:
        B2 = GenerateSHMatrix( opt.L, u0(idi:idf,:) ); % Qc x K
        % ---------
        % Evaluate at the desired direction:
        mu0(idi:idf,1) = sum(B2.*sh,2);     % Qc x 1
        % ---------
    end
    % ---------------------------------------------------------
    % Normalization:
    mu0 = mu0*cnst;
    % ---------------------------------------------------------
else % Full
    B  = B(1,:);                                 % 1 x G, DC coefficient
    atti = atti.^(-(3+nu)/2);                    % Q x G
    if(is_broadcast_available)
        mu0 = sum( atti.*B, 2 );                 % Q x 1
    else
        mu0 = sum( ...
            bsxfun( @(x,y)(x.*y), atti, B ), ... % Q x G
            2 );                                 % Q x 1
    end
    mu0 = mu0*(gamma((nu+3)/2)*sqrt(pi)/(4*pi*pi*opt.tau)^((nu+3)/2));
end
% -------------------------------------------------------------------------
% Remove outliers if necessary:
if(opt.clean>100*eps)
    mu0 = clean_moments_outliers(mu0,opt.clean);
end
% -------------------------------------------------------------------------
% Assign the output
mu = zeros(M*N*P,1);
mu(mask,:) = mu0;
mu = reshape(mu,[M,N,P]);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [type,nu] = micro2moments_type(type,nu)
if(length(type)==1)
    assert(~isempty(nu),'You must provide the moment order if a standard dMRI measurement name is not used');
    assert(isscalar(nu),'The moment order, nu, must be 1x1');
    switch(lower(type))
        case 'f'
            type = 0;
            assert(nu>-3+1e9*eps,'For full moments, nu must be > -3');
        case 'a'
            type = 1;
            assert(nu>-1+1e9*eps,'For axial moments, nu must be > -1');
        case 'p'
            type = 2;
            assert(nu>-2+1e9*eps,'For planar moments, nu must be > -2');
        otherwise
            error(['Cannot parse ',type,' as a valid moment (type: f|a|p)']);
    end
else
    switch(lower(type))
        case 'rtop'
            type   = 0;
            nu     = 0;
        case 'rtpp'
            type   = 1;
            nu     = 0;
        case 'rtap'
            type   = 2;
            nu     = 0;
        case 'qmsd'
            type   = 0;
            nu     = 2;
        otherwise
            error(['Cannot recognize ',type,' as a standar, E(q)-based dMRI measurement']);
    end
end
