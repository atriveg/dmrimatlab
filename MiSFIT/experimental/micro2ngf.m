function ng = micro2ngf( sh, lpar, lperp, f, dti, varargin )
% function ng = micro2ng( sh, lpar, lperp, dti, 'option1', value1, ... )
%
%   Computes the so-called Non-Gaussianity (NG) of diffusion at each voxel
%   according to a linear convolutional model:
%
%      atti(u,b) = Integral_{S}Phi(v)exp(-b*((lpar-lperp)(u*v)^2+lperp))dv.
%
%   The NG is defined as the sine of the angle between the EAP, P(R), and
%   its Gaussian counterpart defined as the diffusion tensor model that 
%   best fits the multi-shell measurements. The angle between two signals 
%   is defined in the common way in terms of their inner product.
%
%   Inputs:
%
%      sh: a MxNxPxK, with K=(L+1)*(L+2)/2 and L>0 even, double array with
%         the coefficients of the ODF obtained with micro2shodf
%      lpar: a MxNxP double array with the parallel diffusivity modeling
%         the impulse response (should fulfill 0<lpar<=Diso). This is
%         obtained with atti2micro
%      lperp: a MxNxP double array with the perpendicular diffusvity
%         modeling the impulse response (should fulfill 0<lerp<lpar). This
%         is obtained with atti2micro
%      f: a MxNxP double array with the partial volume fraction of
%           intra-cellular water (should fulfill 0<=f<=1). If an empty
%           array is passed, then f=1 for all voxels is assumed, so that
%           ones(M,N,P) has the same effect as [].
%      dti: a MxNxPx6 double array with the tensor volume to compute the
%         distance to. The proper way to compute this tensor is:
%         1- Use micro2atti to generate a synthetic attenuation signal from
%            you model, something like:
%            >> atti2 = micro2atti( sh, lpar, lperp, ...
%                   gi(bi<=BMAX,:), bi(bi<=BMAX,:), 'mask', mask );
%            Note only those shells below BMAX should be used (for BMAX
%            something like 2000) since otherwise the tensor model will not
%            be physically consistent.
%         2- Use either signal2dti or:
%            >> sh2 = atti2shadc( atti2, gi(bi<=BMAX,:), ...
%                   bi(bi<=BMAX,:), ...
%                   'mask', mask, 'L', 2 , 'lambda', 0 );
%            >> dti = shadc2dti(sh2, 'mask', mask ,'unroll', false);
%            to estimate the desired tensor.
%
%   Outputs:
%
%      ng: MxNxP double array with the NG at each voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).
%      epsilon: to improve the contrast of the raw NG, a gamma correction
%         can be performed with the form:
%              ng = ng0^(3*epsilon)/(1-3*ng0^epsilon+3*ng0^(2*epsilon));
%         Use empty brackets, [], to avoid this correction and work with
%         the raw NG (default: 0.4).
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
%   Other patameters:
%
%      lambda: the algorithm is based on fitting SH to the (inverse of the
%         squared root of) the determinant of the sum of two matrices. This
%         is the Laplace-Beltrami penalty parameter for this fitting
%         (default: 1.0e-6).
%      chunksz: This parameter is directly passed to the signal2sh function
%         to fit the signal in the basis of SH chunk-by-chunk (see the help
%         therein) (default: 100).

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<5)
    error('At least the sh volume, lpar, lperp, f and dti must be supplied');
end
%%% ----------------
[M,N,P] = size(lpar);
assert(isequal(size(lpar),size(lperp)),'lpar and lperp must be the same size');
if(~isempty(f))
    assert(isequal(size(f),size(lpar)),'lpar and f must be the same size');
else 
    f = ones(M,N,P);
end
%%% ----------------
[M2,N2,P2,K] = size(sh);
assert(isequal([M2,N2,P2],[M,N,P]),'The first 3 dimensions of sh must match those of lpar and lperp');
L = (-3+sqrt(1+8*K))/2;
assert( abs(round(L)-L)<1000*eps, 'This is a weird size for the fourth dimension of sh; make sure it is a SH volume' );
assert( L>=2, 'This method makes no sense with trivial SH volumes with L=0' );
L = round(L);
%%% ----------------
[M2,N2,P2,K2] = size(dti);
assert(isequal([M2,N2,P2],[M,N,P]),'The first 3 dimensions of dti must match those of lpar and lperp');
assert(K2==6,'The fourth dimension of dti must have size 6. Perhaps you missed the ''unroll'', false option in signal2dti?');
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
opt.epsilon = 0.4;      optchk.epsilon = [true,false]; % might be an empty array
% -------------------------------------------------------------------------
opt.chkmod = true;      optchk.chkmod = [true,true];   % always 1x1 boolean
opt.flperp = 0.001;     optchk.flperp = [true,true];   % always 1x1 double
opt.Flperp = 0.999;     optchk.Flperp = [true,true];   % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
% -------------------------------------------------------------------------
opt.lambda = 1.0e-6;    optchk.lambda = [true,true];   % always 1x1 double
opt.chunksz = 100;      optchk.chunksz = [true,true];  % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
if( ~isempty(opt.epsilon) )
    assert(numel(opt.epsilon)==1,'Optional argument ''epsilon'' must be either an empty array or have size 1x1');
end
% -------------------------------------------------------------------------
% Everything seems correct. Try loading the pre-computed weights if
% necessary:
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
% Now structure S has all the information required
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
use_parallel           = use_parallel_test;
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
% Compute the spectrum of the tensor volume:
[u1,u2,u3,l1,l2,l3] = dti2spectrum( dti, 'mask', opt.mask );
l1(l1<0) = 0; l1(l1>opt.ADC0) = opt.ADC0;
l2(l2<0) = 0; l2(l2>opt.ADC0) = opt.ADC0;
l3(l3<0) = 0; l3(l3>opt.ADC0) = opt.ADC0;
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably:
mask  = reshape( opt.mask, [M*N*P,1] ); % M*N*P  x 1
% ---
lpar  = reshape( lpar,     [M*N*P,1] ); % M*N*P  x 1
lpar  = lpar(mask,:);                   % Q x 1
lperp = reshape( lperp,    [M*N*P,1] ); % M*N*P  x 1
lperp = lperp(mask,:);                  % Q x 1
% ---
sh    = reshape( sh, [M*N*P,K] );       % (M*N*P)x((L+1)(L+2)/2)
sh    = sh(mask,:);                     % Qx((L+1)(L+2)/2)
% ---
u1 = reshape(u1,[M*N*P,3]);             % M*N*P x 3
u2 = reshape(u2,[M*N*P,3]);             % M*N*P x 3
u3 = reshape(u3,[M*N*P,3]);             % M*N*P x 3
u1 = u1(mask,:);                        % Q x 3
u2 = u2(mask,:);                        % Q x 3
u3 = u3(mask,:);                        % Q x 3
% ---
l1 = reshape(l1,[M*N*P,1]);             % M*N*P x 1
l2 = reshape(l2,[M*N*P,1]);             % M*N*P x 1
l3 = reshape(l3,[M*N*P,1]);             % M*N*P x 1
l1 = l1(mask);                          % Q x 1
l2 = l2(mask);                          % Q x 1
l3 = l3(mask);                          % Q x 1
% ---
Q = size(lpar,1);
% -------------------------------------------------------------------------
% Time for sanity checks on the micro-structural model:
if(opt.chkmod)
    lpar(lpar>opt.ADC0)    = opt.ADC0;
    lpar(lpar<opt.ADC0/20) = opt.ADC0/20;
    lperp(lperp<lpar*opt.flperp) = opt.flperp.*lpar(lperp<lpar*opt.flperp);
    lperp(lperp>lpar*opt.Flperp) = opt.Flperp.*lpar(lperp>lpar*opt.Flperp);
end
% -------------------------------------------------------------------------
% Compute integral weights by interpolation:
delta = (lpar-lperp);
rho   = lperp./delta;    % Q x 1
rhol  = log(rho);        % Q x 1
rholi = S.rhol;          % 1 x NI
rhol(rhol<rholi(1)) = rholi(1);
wPAli = S.wPAl;          % (LI/2+1) x NI
wPA   = zeros(Q,L/2+1);  % Q x (L/2+1)
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
% -------------------------------------------------------------------------
% Apply these weights to the (squared) SH volume to compute the squared 
% module of the EAP. The term sqrt(4*pi*tau)^(3/2) is simplified
% everywhere:
ptr   = dmri_sh_expand_coeffs(L);           % 1 x K
modE  = sum( sh.*sh.*wPA(:,ptr), 2 );       % Q x 1
modE  = pi*modE./sqrt(delta.*delta.*delta); % Q x 1
% -------------------------------------------------------------------------
% Compute the (inverse of) the squared module of the dti, The term 
% sqrt(4*pi*tau)^(3/2) is simplified everywhere:


    %! EP 
    imodG = sqrt(8*l1.*l2.*l3);                 % Q x 1
    imodGiso = sqrt((opt.ADC0 + l1).*(opt.ADC0 + l2).*(opt.ADC0 + l3));                 % Q x 1


% -------------------------------------------------------------------------
% Sample the determinant of dti plus the elemental tensor in the unit
% sphere with sufficient accuracy.
% 1. Use at least twice as samples as the number of coefficients to fit:
lev = log2(max((2*K-16)/15,1)) + 2; % See the help on icosamplesSphere
vi  = icosamplesSphere(ceil(lev),...
    'O1',true,'iters',20,'verbose',false); % G x 3
% Sample the function:
l1   = l1+lperp; % Q x 1
l2   = l2+lperp; % Q x 1
l3   = l3+lperp; % Q x 1
prj1 = ( u1(:,1)*vi(:,1)' + u1(:,2)*vi(:,2)' + u1(:,3)*vi(:,3)' ); % Q x G
prj2 = ( u2(:,1)*vi(:,1)' + u2(:,2)*vi(:,2)' + u2(:,3)*vi(:,3)' ); % Q x G
prj3 = ( u3(:,1)*vi(:,1)' + u3(:,2)*vi(:,2)' + u3(:,3)*vi(:,3)' ); % Q x G
if(is_broadcast_available)
    prj1 = prj1.*prj1./l1; % Q x G
    prj2 = prj2.*prj2./l2; % Q x G
    prj3 = prj3.*prj3./l3; % Q x G
else
    prj1 = bsxfun( @(x,y)(x./y), prj1.*prj1, l1 ); % Q x G
    prj2 = bsxfun( @(x,y)(x./y), prj2.*prj2, l2 ); % Q x G
    prj3 = bsxfun( @(x,y)(x./y), prj3.*prj3, l3 ); % Q x G
end
detp  = prj1 + prj2 + prj3;    % Q x G
if(is_broadcast_available)
    detp = delta.*detp + 1;    % Q x G
    detp = detp.*(l1.*l2.*l3); % Q x G
else
    detp = bsxfun(@(x,y)(x.*y),delta,detp) + 1;  % Q x G
    detp = bsxfun(@(x,y)(x.*y),detp,l1.*l2.*l3); % Q x G
end
detp = 1./sqrt(detp);              % Q x G
% -------------------------------------------------------------------------
% Compute the SH for this spherical function:
sh2 = signal2sh( reshape(detp,[Q,1,1,size(vi,1)]), vi, ...
    'L', L, 'lambda', opt.lambda, 'chunksz', opt.chunksz );
sh2 = reshape(sh2,[Q,K]);

f = reshape( f, [M*N*P,1] );
f = f(mask, :);

% Isotropic Fraction 
fiso = (1-f).*(1-f);
iso  =  fiso/(2*opt.ADC0)^(3/2);

% Anisotropic Fraction 
ptr     = dmri_sh_expand_coeffs(L);
f2      = f.*f;
aniso   = pi*f2.*sum(sh.*sh.*wPA(:, ptr), 2)./sqrt(delta.*delta.*delta);

% Mix Fraction 
fs      = f.*(1-f);
mix_num = 4*sqrt(pi)*fs.*sh(:,1);
mix_den = (opt.ADC0+lperp).*sqrt(opt.ADC0+lpar);
mix     = mix_num./mix_den;

denMOD     = iso + mix + aniso;


% -------------------------------------------------------------------------
% Compute the cross-product and the NG:
cp  = sum(sh.*sh2,2);                            % Q x 1
%ng0 = sqrt(min(max(1-imodG.*cp.*cp./modE,0),1)); % Q x 1


numMOD = ((1-f)./imodGiso + f.*cp);
numMOD = numMOD.*numMOD;
numMOD = numMOD.*imodG;
%numMOD = numMOD./modE;
ng0 = sqrt(min(max(1-numMOD./denMOD,0),1)); % Q x 1


% -------------------------------------------------------------------------
% Compute the gamma-correction if necessary:
if(~isempty(opt.epsilon))
    warning('The gamma-coontrast enhancement should not be applied to the NG');
    ng0 = ng0.^(3*opt.epsilon)./ ...
        ( 1 - 3*ng0.^(opt.epsilon) + 3*ng0.^(2*opt.epsilon));
end
% -------------------------------------------------------------------------
% Reshape things back and exit:
ng = zeros(M*N*P,1);
ng(mask,:) = ng0;
ng = reshape(ng,[M,N,P]);
% -------------------------------------------------------------------------
