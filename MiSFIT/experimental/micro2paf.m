function pa = micro2paf( sh, lpar, lperp, f, varargin )
% function pa = micro2pa( sh, lpar, lperp, 'option1', value1, ... )
%
%   Computes the so-called Propagator Anisotropy of diffusion at each voxel
%   according to a linear convolutional model:
%
%      atti(u,b) = (1-f)exp(-bDo)+f*Integral_{S}Phi(v)exp(-b*((lpar-lperp)(u*v)^2+lperp))dv.
%
%   The PA is defined as the sine of the angle between the EAP, P(R), and
%   its isotropic counterpart defined as the spherical average of P(R). The
%   angle between two signals is defined in the common way in terms of
%   their inner product.
%   Inputs:
%
%      sh: a MxNxPxK, with K=(L+1)*(L+2)/2 and L>0 even, double array with
%         the coefficients of the ODF obtained with micro2shodf.
%         Alternatively, this argument may be an empty array, [], if the
%         user requests the 'micro-PA'.
%      lpar: a MxNxP double array with the parallel diffusivity modeling
%         the impulse response (should fulfill 0<lpar<=Diso). This is
%         obtained with atti2micro.
%      lperp: a MxNxP double array with the perpendicular diffusvity
%         modeling the impulse response (should fulfill 0<lerp<lpar). This
%         is obtained with atti2micro.
%      f: a MxNxP double array with the partial volume fraction of
%           intra-cellular water (should fulfill 0<=f<=1). If an empty
%           array is passed, then f=1 for all voxels is assumed, so that
%           ones(M,N,P) has the same effect as [].
%
%   Outputs:
%
%      pa: MxNxP double array with PA at each voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).
%      micro: wether (true) or not (false) compute the microscopic version
%         of the PA, which depends only on lpar and lperp but not on sh.
%         Accordingly, sh input may be left empty in this case (default:
%         false).
%      epsilon: to improve the contrast of the raw PA, a gamma correction
%         is usually performed with the form:
%              pa = pa0^(3*epsilon)/(1-3*pa0^epsilon+3*pa0^(2*epsilon));
%         Use empty brackets, [], to avoid this correction and work with
%         the raw PA (default: 0.4).
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

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<4)
    error('At least the sh volume, lpar, lperp and f must be supplied');
end
[M,N,P] = size(lpar);
assert(isequal(size(lpar),size(lperp)),'lpar and lperp must be the same size');
if(~isempty(f))
    assert(isequal(size(f),size(lpar)),'lpar and f must be the same size');
else
    f = ones(M,N,P);
end
% SH will be checked later on, since it might be left empty depending on
% the type of moment to compute
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
opt.micro = false;      optchk.micro = [true,true];    % always 1x1 boolean
opt.epsilon = 0.4;      optchk.epsilon = [true,false]; % might be an empty array
% -------------------------------------------------------------------------
opt.chkmod = true;      optchk.chkmod = [true,true];   % always 1x1 boolean
opt.flperp = 0.001;     optchk.flperp = [true,true];   % always 1x1 double
opt.Flperp = 0.999;     optchk.Flperp = [true,true];   % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Check the sh input argument if necessary:
if(~opt.micro) % Regular PA
    % ------------------------------
    % Check if sh has a proper size:
    [M2,N2,P2,K] = size(sh);
    assert(isequal([M2,N2,P2],[M,N,P]),'If the ''micro'' option is false, the first 3 dimensions of sh must match those of lpar and lperp');
    % ------------------------------
    % Check if the fourth dimension of sh is correct, and look for the SH
    % expansion order:
    L = (-3+sqrt(1+8*K))/2;
    assert( abs(round(L)-L)<1000*eps, 'This is a weird size for the fourth dimension of sh; make sure it is a SH volume' );
    assert( L>=2, 'This method makes no sense with trivial SH volumes with L=0' );
end
% Check optional input argument epsilon:
if( ~isempty(opt.epsilon) )
    assert(numel(opt.epsilon)==1,'Optional argument ''epsilon'' must be either an empty array or have size 1x1');
end
% -------------------------------------------------------------------------
% Everything seems correct. Try loading the pre-computed weights if
% necessary:
if(~opt.micro)
    if(exist('dmri_PA_weights.mat','file')==2)
        S = load('dmri_PA_weights.mat');
        if(S.L<L)
            fprintf(1,'Your sh volume contain SH coefficients up to order %d but only\n',L);
            fprintf(1,'those factors up to order %d are present in dmri_PA_weights.mat.\n',S.L);
            fprintf(1,'I am now running dmri_compute_PA_weights.m to compute the required\n');
            fprintf(1,'weights. This will take a while, but you will not get this message\n');
            fprintf(1,'again unless you request weights for orders L larger than those\n');
            fprintf(1,'stored in dmri_PA_weights.mat or you explicitly delete this file.\n');
            fprintf(1,'Proceeding... ');
            dmri_compute_PA_weights(L);
            fprintf(1,'[DONE]\n');
            S = load('dmri_PA_weights.mat');
        end
    else % Weights not already computed!
        fprintf(1,'I cannot find file dmri_PA_weights.mat with the precomputed\n');
        fprintf(1,'weigths. I am now running dmri_compute_PA_weights.m to\n');
        fprintf(1,'create it, which may take a while. You will not get this\n');
        fprintf(1,'message again unless you request weights for orders L\n');
        fprintf(1,'larger than those stored in dmri_PA_weights.mat or you\n');
        fprintf(1,'explicitly delete this file. Proceeding... ');
        dmri_compute_PA_weights( max(L,12) );
        fprintf(1,'[DONE]\n');
        S = load('dmri_PA_weights.mat');
    end
    % Now structure S ha all the information required
end
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
use_parallel           = use_parallel_test;
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably:
mask  = reshape( opt.mask, [M*N*P,1] );
lpar  = reshape( lpar,     [M*N*P,1] ); lpar  = lpar(mask,:);  % Qx1
lperp = reshape( lperp,    [M*N*P,1] ); lperp = lperp(mask,:); % Qx1
if(~opt.micro)
    sh = reshape( sh, [M*N*P,K] );   % (M*N*P)x((L+1)(L+2)/2)
    sh = sh(mask,:);                 % Qx((L+1)(L+2)/2)
end
Q     = size(lpar,1);
% -------------------------------------------------------------------------
% Time for sanity checks on the micro-structural model:
if(opt.chkmod)
    lpar(lpar>opt.ADC0)    = opt.ADC0;
    lpar(lpar<opt.ADC0/20) = opt.ADC0/20;
    lperp(lperp<lpar*opt.flperp) = opt.flperp.*lpar(lperp<lpar*opt.flperp);
    lperp(lperp>lpar*opt.Flperp) = opt.Flperp.*lpar(lperp>lpar*opt.Flperp);
end
% -------------------------------------------------------------------------
% Compute the basic value of the PA in each case
if(opt.micro) % Micro-PA
    rho = (lpar-lperp)/2./sqrt(lpar.*lperp);
    pa0  = sqrt( 1 - min(atan(rho)./rho,1) );
else
    % --------------------------------------------------------------
    % Compute rho, which is the micro-structure descriptor:
    rho   = lperp./(lpar-lperp); % Q x 1
    rhol  = log(rho);            % Q x 1
    % --------------------------------------------------------------
    % Interpolate the weigths in the logarithmic domain, where they are
    % smooth enough:
    rholi = S.rhol;              % 1 x NI
    rhol(rhol<rholi(1)) = rholi(1);
    wPAli = S.wPAl;              % (LI/2+1) x NI
    wPA   = zeros(Q,L/2+1);      % Q x (L/2+1)
    if(use_parallel)
        parfor l=1:L/2+1
            wPA(:,l) = exp( interp1(rholi,wPAli(l,:),rhol) ); % Q x 1
        end
    else
        for l=1:L/2+1
            wPA(:,l) = exp( interp1(rholi,wPAli(l,:),rhol) ); % Q x 1
        end
    end
    wPA(isinf(wPA)) = 0;
    wPA(isnan(wPA)) = 0;

    % --------------------------------------------------------------
    % Apply these weights to the (squared) SH volume:
    sh = sh.*sh;                    % Q x K
    
    % --------------------------------------------------------------
    % Free-water Fraction

    f2      = f.*f;
    f2      = reshape( f2, [M*N*P,1] );
    f2      = f2(mask, :);


    fiso    = (1-f).*(1-f);
    fiso    = reshape( fiso, [M*N*P,1] );
    fiso    = fiso(mask, :);
    iso     = fiso/(2*opt.ADC0)^(3/2);


    ptr     = dmri_sh_expand_coeffs(L);
    delta   = (lpar-lperp);
    aniso   = pi*f2.*sum(sh.*wPA(:, ptr), 2)./sqrt(delta.*delta.*delta);
    aniso0  = pi*f2.*sh(:,1).*wPA(:,1)./sqrt(delta.*delta.*delta);

    fs      = f.*(1-f);
    fs      = reshape( fs, [M*N*P,1] );
    fs      = fs(mask, :);
    mix_num = 4*sqrt(pi)*fs.*sh(:,1);
    mix_den = (opt.ADC0+lperp).*sqrt(opt.ADC0+lpar);
    mix     = mix_num./mix_den;

    num     = iso + mix + aniso0;
    den     = iso + mix + aniso;
    pa0      = sqrt(1-num./den);
    
    
%     Old code!    
%     num = sh(:,1).*wPA(:,1);        % Q x 1, numerator
%     ptr = dmri_sh_expand_coeffs(L); % 1 x K
%     den = sum( sh.*wPA(:,ptr), 2 ); % Q x 1
%     pa0 = sqrt(1-num./den);         % Q x 1


    % --------------------------------------------------------------
    % Sanity checks for weird values of rho:
    pa0(rho>S.rho(end)) = 0;
    pa0(isnan(pa0))     = 0;
    % --------------------------------------------------------------
end
% -------------------------------------------------------------------------
% Compute the gamma-correction if necessary:
if(~isempty(opt.epsilon))
    pa0 = pa0.^(3*opt.epsilon)./ ...
        ( 1 - 3*pa0.^(opt.epsilon) + 3*pa0.^(2*opt.epsilon));
end
% -------------------------------------------------------------------------
% Reshape things back and exit:
pa = zeros(M*N*P,1);
pa(mask,:) = pa0;
pa = reshape(pa,[M,N,P]);
% -------------------------------------------------------------------------
