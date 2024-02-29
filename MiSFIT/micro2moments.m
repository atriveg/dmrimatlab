function mu = micro2moments( sh, lpar, lperp, f, varargin )
% function mu = micro2moments( sh, lpar, lperp, f, 'option1', value1, ... )
%
%   Computes moments (either full, axial, or planar) of arbitrary orders
%   over either the attenuation signal E(q) or the diffusion propagator
%   P(R) according to a linear convolutional model:
%
%      atti(u,b) = Integral_{S}Phi(v)exp(-b*((lpar-lperp)(u*v)^2+lperp))dv,
%
%   where lpar and lperp parameterize the impulse response as an elemental
%   rank-2 tensor and the ODF, Phi(v), is modeled through its SH expansion
%   over the unit sphere. Either type of moments are defined as:
%
%      mu_{full,E}^{nu}       = Integral_{R^3} ||q||^{nu} E(q) dq^3
%      mu_{axial,E}^{nu}(u0)  = Integral_{R}   t^{nu} E(t*u0) dt
%      mu_{planar,E}^{nu}(u0) = Integral_{v \perp u0} ||v||^{nu} E(v) dv^2
%      mu_{full,P}^{nu}       = Integral_{R^3} ||R||^{nu} P(R) dR^3
%      mu_{axial,P}^{nu}(r0)  = Integral_{R^3} t^{nu} P(t*r0) dt
%      mu_{planar,P}^{ru}     = Integral_{s \perp r0} ||s||^{nu} P(s) ds^2
%
%   where nu>-3 for full moments, nu>-1 for axial moments, and nu>-2 for
%   planar moments. In all cases, nu can be any real number, though the
%   computation of integer moments is computationally more efficient and
%   reliable.
%
%   Most of the common diffusion measurements are described with this
%   scheme: RTOP = mu_{full,E}^0, RTPP = mu_{axial,E}^0(u_max), RTAP = 
%   mu_{planar,E}^0(u_max), QMSD = mu_{full,E}^2, MSD = mu_{full,P}^2...
%
%   Inputs:
%
%      sh: a MxNxPxK, with K=(L+1)*(L+2)/2 and L>0 even, double array with
%         the coefficients of the ODF obtained with micro2shodf.
%         Alternatively, this argument may be an empty array, [], if the
%         moment to be computed is a 'full' one.
%      lpar: a MxNxP double array with the parallel diffusivity modeling
%         the impulse response (should fulfill 0<lpar<=Diso). This is
%         obtained with atti2micro.
%      lperp: a MxNxP double array with the perpendicular diffusvity
%         modeling the impulse response (should fulfill 0<lerp<lpar). This
%         is obtained with atti2micro.
%      f: a MxNxP double array with the partial volume fraction of
%         intra-cellular water (should fulfill 0<=f<=1). If an empty array
%         is passed, then f=1 for all voxels is assumed, so that
%         ones(M,N,P) has the same effect as [].
%
%   Outputs:
%
%      mu: MxNxP double array with the moment requested.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      type: 1x2 string array, one of 'Ef', 'Ea', 'Ep', 'Pf', 'Pa', Pp',
%         the first character indicating the signal the moment is computed
%         over and the second one if it is either full, axial, or planar.
%         Alternatively, you can specify either of 'rtop', 'rtpp', 'rtap',
%         'qmsd', or 'msd' (default: 'rtop').
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
%   Sanity checks on the micro-structure model:
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
%      ADC0: estimated diffusivity of free water at body temperature 
%         (Diso). Should use the same as in atti2micro (default: 3.0e-3).
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
%         applied, the closer to 100 the more agressive. **You should use
%         this flag just for visualization purposes, but not for actual 
%         quantitative analyses** (default: 0).

% Check the mandatory input arguments:
if(nargin<4)
    error('At least the sh volume, lpar, lperp and f must be supplied');
end
[M,N,P] = size(lpar);
assert(isequal(size(lpar),size(lperp)),'lpar and lperp must be the same size');
if(~isempty(f))
    assert(isequal(size(f),[M,N,P]),'lpar and f must be the same size');
end
% SH will be checked later on, since it might be left empty depending on
% the type of moment to compute
% -------------------------------------------------------------------------
opt.type = 'rtop';      optchk.type = [true,false];    % Variable length string
opt.nu = 0;             optchk.nu = [true,false];      % Double, either 1x1 or empty
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
opt.u0 = [];            optchk.u0 = [true,false];      % Double, either MxNxPx3 or empty
% -------------------------------------------------------------------------
opt.chkmod = true;      optchk.chkmod = [true,true];   % always 1x1 boolean
opt.flperp = 0.001;     optchk.flperp = [true,true];   % always 1x1 double
opt.Flperp = 0.999;     optchk.Flperp = [true,true];   % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
% -------------------------------------------------------------------------
opt.tau = 70.0e-3;      optchk.tau = [true,true];      % always 1x1 double
opt.chunksz = 256;      optchk.chunksz = [true,true];  % always 1x1 double
opt.clean = 0;          optchk.clean = [true,true];    % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Check consistency of outlier rejection:
assert( opt.clean>=0 && opt.clean<=100, '''clean'' must be in the range [0,100]' );
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
% Determine the kind and the order of the moment to compute:
[signal,type,nu] = micro2moments_type(opt.type,opt.nu);
% -------------------------------------------------------------------------
% Check now sh if necessary
if(type>0) % Axial or planar, we really need to check
    if(isempty(sh))
        error('An actual SH volume is required to compute the requested type of moment, so you cannot leave empty input sh');
    else
        [M2,N2,P2,K] = size(sh);
        assert(isequal([M,N,P],[M2,N2,P2]),'The first three dimensions of sh should match those of lpar and lperp');
        L = (-3+sqrt(1+8*K))/2;
        assert( abs(round(L)-L)<1000*eps, 'This is a weird size for the fourth dimension of sh; make sure it is a SH volume' );
        assert( L>=2, 'This method makes no sense with trivial SH volumes with L=0' );
    end
end
% -------------------------------------------------------------------------
% If the moment is axial or planar, u0 must be checked, too
if( type>0 )
    if(~isempty(opt.u0))
        [M2,N2,P2,K2] = size(opt.u0);
        assert(isequal([M,N,P],[M2,N2,P2]),'The first three dimensions of u0 should match those of lpar and lperp');
        assert(K2==3,'u0 must have size MxNxPx3, so that the last dimension stands for a 3x1 unit vector');
    else
        tens   = shadc2dti(sh,'mask',opt.mask,...
            'chunksz',opt.chunksz,'unroll',false );   % MxNxPx6
        opt.u0 = dti2xyz( tens, 'mask', opt.mask );   % MxNxPx3
    end
end
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably:
mask  = reshape( opt.mask, [M*N*P,1] );
lpar  = reshape( lpar,     [M*N*P,1] ); lpar  = lpar(mask,:);  % Qx1
lperp = reshape( lperp,    [M*N*P,1] ); lperp = lperp(mask,:); % Qx1
if(type>0)
    sh = reshape( sh, [M*N*P,K] );     sh = sh(mask,:); % Qx((L+1)(L+2)/2)
    u0 = reshape( opt.u0, [M*N*P,3] ); u0 = u0(mask,:); % Qx3
end
% -------------------------------------------------------------------------
% Time for sanity checks on the micro-structural model:
if(opt.chkmod)
    lpar(lpar>opt.ADC0)    = opt.ADC0;
    lpar(lpar<opt.ADC0/20) = opt.ADC0/20;
    lperp(lperp<lpar*opt.flperp) = opt.flperp.*lpar(lperp<lpar*opt.flperp);
    lperp(lperp>lpar*opt.Flperp) = opt.Flperp.*lpar(lperp>lpar*opt.Flperp);
end
% -------------------------------------------------------------------------
% Compute the integrals we will need:
if(signal>0) % P(R)
    switch(type)
        case 0 % Full moment
            ints = dmri_compute_Pkernel_integrals(lpar,lperp,nu+3,0);                                  % Qx1
            CONST = gamma((nu+3)/2)*(4*opt.tau*lperp).^(nu/2).*sqrt(lperp)./sqrt(pi*lpar);             % Qx1
        case 1 % Axial moment
            ints  = dmri_compute_Pkernel_integrals(lpar,lperp,nu+1,L/2);                               % Qx(L/2+1)
            CONST = gamma((nu+1)/2)*(4*opt.tau*lperp).^((nu-2)/2).*sqrt(lperp)./sqrt(pi*pi*pi*lpar);   % Qx1
        case 2 % Planar moment
            ints  = dmri_compute_Pkernel_integrals(lpar,lperp,nu+2,L/2);                               % Qx(L/2+1)
            CONST = gamma((nu+2)/2)*(4*opt.tau*lperp).^((nu-1)/2).*sqrt(lperp)./sqrt(pi*pi*pi*lpar)/2; % Qx1
    end
else % E(q)
    switch(type)
        case 0 % Full moment
            ints  = dmri_compute_Ekernel_integrals(lpar,lperp,nu+3,0);       % Qx1
            CONST = gamma((nu+3)/2)*pi./(4*pi*pi*opt.tau*lperp).^((nu+3)/2); % Qx1
        case 1 % Axial moment
            ints  = dmri_compute_Ekernel_integrals(lpar,lperp,nu+1,L/2);     % Qx(L/2+1)
            CONST = gamma((nu+1)/2)./(4*pi*pi*opt.tau*lperp).^((nu+1)/2);    % Qx1
        case 2 % Planar moment
            ints = dmri_compute_Ekernel_integrals(lpar,lperp,nu+2,L/2);      % Qx(L/2+1)
            CONST = gamma((nu+2)/2)./(4*pi*pi*opt.tau*lperp).^((nu+2)/2)/2;  % Qx1
    end
end
% For axial and planar moments, combine into integrals with Legendre
% polynomials:
if(type>0)
    ints = dmri_compute_Legendre_projection(ints); % Qx(L/2+1)
end
% -------------------------------------------------------------------------
% Effectively compute the moments
if(type>0) % Axial or planar
    % ---------------------------------------------------------
    % Correct SH coefficients with the convolution coefficients
    sh = sh.*ints(:,micro2moments_expandF(L)); % Qx(L+1)(L+2)/2
    % ---------------------------------------------------------
    % If neccesary, apply the FRT:
    if(type>1)
        frt = diag(GenerateFRTMatrix(L))';    % 1x(L+1)(L+2)/2
        if(is_broadcast_available)
            sh = sh.*frt;                     % Qx(L+1)(L+2)/2
        else
            sh = bsxfun(@(x,y)(x.*y),sh,frt); % Qx(L+1)(L+2)/2
        end
    end
    % ---------------------------------------------------------
    % Evaluate at desired directions
    Q   = size(sh,1);
    mu0 = zeros(Q,1);
    for ck=1:ceil(Q/opt.chunksz)
        % ------------------------------------------
        % Characterize the present chunk:
        idi = (ck-1)*opt.chunksz+1;
        idf = min(ck*opt.chunksz,Q);
        B   = GenerateSHMatrix( L, u0(idi:idf,:) ); % Qcx(L+1)(L+2)/2
        mu0(idi:idf,1) = sum(sh(idi:idf,:).*B,2);
    end
    % ---------------------------------------------------------
    % Correct with the appropriate constant:
    mu0 = CONST.*mu0;   % Qx1
    % ---------------------------------------------------------
else % Full
    mu0 = CONST.*ints;  % Qx1
end
% -------------------------------------------------------------------------
% In case an actual volume with the free-water compartment was provided,
% use it here:
if(~isempty(f))
    f   = reshape( f, [M*N*P,1] );                            % (M*N*P)x1
    f   = f(mask,:);                                          % Qx1
    muc = isotropic_moment(signal,type,nu,opt.ADC0,opt.tau);  % 1x1
    mu0 = f.*mu0 + (1-f)*muc;                                 % Qx1                          
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
function [signal,type,nu] = micro2moments_type(type,nu)
if(length(type)==2)
    assert(~isempty(nu),'You must provide the moment order if a standard dMRI measurement name is not used');
    assert(isscalar(nu),'The moment order, nu, must be 1x1');
    switch(upper(type(1)))
        case 'E'
            signal = 0;
        case 'P'
            signal = 1;
        otherwise
            error(['Cannot parse first character of ',type,', ',type(1),' as a valid signal (type E|P)']);
    end
    switch(lower(type(2)))
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
            error(['Cannot parse second character of ',type,', ',type(2),' as a valid moment (type f|a|p)']);
    end
else
    switch(lower(type))
        case 'rtop'
            signal = 0;
            type   = 0;
            nu     = 0;
        case 'rtpp'
            signal = 0;
            type   = 1;
            nu     = 0;
        case 'rtap'
            signal = 0;
            type   = 2;
            nu     = 0;
        case 'msd'
            signal = 1;
            type   = 0;
            nu     = 2;
        case 'qmsd'
            signal = 0;
            type   = 0;
            nu     = 2;
        otherwise
            error(['Cannot recognize ',type,' as a standar dMRI measurement']);
    end
end

% -------------------------------------------------------------------------
function ptr  = micro2moments_expandF( L )
% Expands (l,0) SH coefficients to full (l,m), m=-l..l coefficients
ptr    = zeros(1,(L+1)*(L+2)/2);
ptr(1) = 1;
pos    = 1;
for l=2:2:L
    nl  = 2*l+1;
    ptr(1,pos+1:pos+nl) = l/2+1;
    pos = pos +nl;
end

% -------------------------------------------------------------------------
function mu = isotropic_moment( signal, type, nu, ADC0, tau )

K = 4*tau*ADC0;
if(signal>0) % P(R)
    switch(type)
        case 0 % Full moment
            mu = 2*K^(nu/2)*gamma((nu+3)/2) / sqrt(pi);
        case 1 % Axial moment
            mu = K^(nu/2-1)*gamma((nu+1)/2) / sqrt(pi^3);
        case 2 % Planar moment
            mu = K^((nu-1)/2)*gamma((nu+2)/2) / sqrt(pi);
    end
else         % E(q)
    switch(type)
        case 0 % Full moment
            mu = 2*gamma((nu+3)/2) / (pi^(nu+2)*K^((nu+3)/2));
        case 1 % Axial moment
            mu = gamma((nu+1)/2) / (pi^(nu+1)*K^((nu+1)/2));
        case 2 % Planar moment
            mu = gamma((nu+2)/2) / (pi^(nu+1)*K^((nu+2)/2));
    end
end
