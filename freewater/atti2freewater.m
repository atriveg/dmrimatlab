function [f,SH] = atti2freewater( atti, gi, bi, varargin )
% function [f,SH] = atti2freewater( atti, gi, bi, 'opt1', value1, 'opt2', value2, ... )
%
%   Fits the mono-exponential ADC and free water fraction estimated from a
%   DWI model following a dumb exhaustive search/least squares fit
%   procedure similar the one described by Koay et all (this generalizes
%   the tensor model to a general mono-exponential model):
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradient over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that multi-shell acquisitions are allowed.
%
%   The signal S_i is assumed to follow a two-compartment mono-exponential
%   model:
%
%      S_i(gi,bi) = f*S_0*exp(-bi*ADC(gi)) + (1-f)*S_0*exp(-bi*ADC0),
%
%   that is further regularized with the additional condition:
%
%      sqrt(tau)*(f-f0) = 0;
%
%   where f is the partial volume fraction of non-free (intra-cellular)
%   water and ADC0 is the diffusivity of free water at body temperature,
%   which has been empirically fixed near 3.0e-3 mm^2/s.  ADC(gi) is the 
%   apparent diffusion coefficient of non-free water, a function defined 
%   over the unit sphere to be fit at each image voxel in the basis of
%   even spherical harmonics up to order L to yield a total of
%   M=(L+1)*(L+2)/2 SH coefficients. f0 is an initial estimate of the
%   actual value of f, and tau is a Tikhonov regularization parameter to
%   fit the model. The value of f0 are subsequently updated as iterations
%   proceed, so that f0 at iteration n+1 equals f at iteration n.
%
%   Then:
%
%         (S_i/S_0) * (1-f*(S_0/S_i)*exp(-bi*ADC_free)) / (1-f)
%                     = exp(-bi*ADC(gi))
%      -> [-log(dwi*(1-f*exp(-bi*ADC)/dwi)/(1-f))]/bi
%                     = sum_{m=1}^{M} SH(m)*Y_m(gi) [Linear Least Squares]
%
%
%   which is once again regularized with the usual Tikhonov redundant
%   condition applied to the energy og the Laplacian of the ADC:
%
%      sqrt(lambda)*laplacian_eigenvalue(m)*SH(m) = 0, m=1...M
%
%   The two outputs of the function are:
%
%      f: a MxNxP aray in the range (0,1) with the partial volume fraction
%         of non-free (i.e. fiber-confined) water at each imaged voxel.
%      SH: a MxNxPx((L+1)*(L+2)/2) array with the coefficients computed at 
%         each voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   For the representation of the DWI in the basis of Spherical Harmonics:
%      L: an even integer with the maximum order of the SH to be used
%         (default: 6).
%      lambda: the Tikhonov regularization parameter for the linear least
%         squares problem that fits SH coeffcients to the without-free-
%         water ADC. It penalizes the energy of the Laplacian
%         (default 0.006).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      ADC0: estimated diffusivity of free water at body temperature. Do
%         not change this value unless youy have a good reason to do so
%         (default: 3.0e-3).
%      pos: wether (true) or not (false) avoiding estimates of f leading to
%        values of the attenuation signal outside the range [tl,tu]
%        (default: false).
%   For the multi-resolution greedy search algorithm:
%      O1: wether (true) or not (false) perform a first step of searching
%         over the whole domain of (discretized) f for a first solution. If
%         this is not done, a constant value f0=0.5 will be used as the
%         first iteration for O2 and/or O3 (default: true).
%      fmin, fmax: the lower and upper limits where f will be looked for at
%         each image voxel, so that fmin should be 0 and fmax should be 1 
%         (or close to 1) (default: 0.1, 1).
%      fbins: how many evenly spaced values of f are probed within its
%         search range in search of the optimal one. This is interpreted as
%         a radius, so that a total of 2*fbins+1 will be probed, fbins 
%         before the current value, the current value, and fbins after the
%         current value (default: 5 -> 11 values are probed).
%      fdepth: how many times the search range of f is refined around the
%         optimal value found at the previous depth level (default: 2).
%   For the Newton-Raphson optimziation algorithm:
%      O2: wether (true) or not (false) use the additional Newton-Raphson
%          optimization to further fit the signal model described above.
%          Note this step may be very, very slow (default: false).
%      tau: the Laplace-Beltrami penalty for the square root of the ADC
%          (default: 0.001)
%      fnriters: the maximum number of Newton-Raphson iterations to wait
%          for convergence until premature exit (default: 100).
%      fnrth: iterations stop (i.e. convergence is detected) if the
%          absolute difference between the previous value of f and the
%          current one is below this threshold (default: 1.0e-4).
%      t0: penalty term in the log-barrier method to ensure f lays within
%          the range (0,1) (default: 0.01). This is the intial value, it
%          decreases as the iteration proceed.
%      rcm: minimum allowed condition number for the matrix to be inverted
%          at each Newton-Raphson (or alike) step (default: 1.0e-6).
%      lmb: initial Levenberg-Marquardt parameter for each processed voxel.
%          If the Newton-Raphson step is successful (the cost decreases),
%          it will be divided by 2 so that iterations become more alike to
%          Newton-Raphson. If not (the cost does not decrease) it is
%          multiplied by 2 so taht iterations become more alike to
%          gradient-descent (default: 1.0).
%   Other general options:
%      chunksz: the LLS problem reduces to the product of the dwi signal
%         by an inverse matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the LLS
%         inverse at eaxh step, hence taking advantage of matlab's
%         capabilities (default: 100).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      verbose: wether (true) or not (false) show a sort of progress bar in
%         the command line (default: false).

% Check the mandatory input argments:
if(nargin<3)
    error('At lest the dwi volume, the gradient table, and the b-values must be supplied');
end
[M,N,P,G] = size(atti);
if(~ismatrix(gi)||~ismatrix(bi))
    error('gi and bi must be 2-d matlab matrixes');
end
if(size(gi,1)~=G)
    error('The number of rows in gi must match the 4-th dimension of dwi');
end
if(size(gi,2)~=3)
    error('The gradients table gi must have size Gx3');
end
if(size(bi,1)~=G)
    error('The number of b-values bi must match the 4-th dimension of dwi');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end

% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.L = 6;              optchk.L = [true,true];        % always 1x1 double
opt.lambda = 0.05;      optchk.lambda = [true,true];   % always 1x1 double
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
opt.pos = false;        optchk.pos  = [true,true];     % always 1x1 boolean
% -------------------------------------------------------------------------
opt.O1 = true;          optchk.O1 = [true,true];       % always 1x1 boolean
opt.fmin = 0.01;        optchk.fmin = [true,true];     % always 1x1 double
opt.fmax = 1;           optchk.fmax = [true,true];     % always 1x1 double
opt.fbins = 5;          optchk.fbins = [true,true];    % always 1x1 double
opt.fdepth = 2;         optchk.fdepth = [true,true];   % always 1x1 double
% -------------------------------------------------------------------------
opt.O2 = false;         optchk.O2 = [true,true];       % always 1x1 boolean
opt.tau = 0.001;        optchk.tau = [true,true];      % always 1x1 double
opt.fnriters = 100;     optchk.fnriters = [true,true]; % always 1x1 double
opt.fnrth = 1.0e-4;     optchk.fnrth = [true,true];    % always 1x1 double
opt.t0 = 0.01;          optchk.t0 = [true,true];       % always 1x1 double
opt.rcm = 1.0e-6;       optchk.rcm = [true,true];      % always 1x1 double
opt.lmb = 1.0;          optchk.lmb = [true,true];      % always 1x1 double
% -------------------------------------------------------------------------
opt.chunksz = 100;      optchk.chunksz = [true,true];  % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the size as the image field
opt.verbose = false;    optchk.verbose = [true,true];  % always 1x1 boolean
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Sanity check over the dwi volume:
atti(atti<opt.tl) = opt.tl;
atti(atti>opt.tu) = opt.tu;

% Unroll the DWI image, the output, the SH coefficients, and the mask for 
% easy processing:
atti = reshape(atti,[M*N*P,G]);
mask = reshape(opt.mask,[M*N*P,1]);
f    = zeros(M*N*P,1);
SH   = zeros(M*N*P,(opt.L+1)*(opt.L+2)/2);

if(opt.O1)
    f(mask) = atti2freewaterO1( atti(mask,:), gi, bi, opt );
else
    f(mask) = ((opt.fmax+opt.fmin)/2)*ones(imagesz);
end

adci  = compute_corrected_atti(atti(mask,:),f(mask), ...
    bi,opt.tu,opt.tl,opt.ADC0);
SHtmp = signal2sh( reshape(adci,[size(adci,1),1,1,G]), gi, ...
    'L', opt.L, 'lambda', opt.lambda, 'chunksz', opt.chunksz );
SH(mask,:) = reshape(SHtmp,[size(SHtmp,1),size(SHtmp,4)]);

if(opt.O2)
    [f,SH] = atti2freewaterO2( atti, gi, bi, opt, f, SH, mask );
end

% Reshape back the outputs to their proper sizes:
SH = reshape(SH,[M,N,P,(opt.L+1)*(opt.L+2)/2]);
f  = reshape(f, [M,N,P]);

%%% -----------------------------------------------------------------------
function f = atti2freewaterO1( atti, gi, bi, opt )
% atti: RxG (R<=M*N*P, all voxels where mask=true)
% gi:   Gx3
% bi:   Gx1
% f:    Rx1
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
R = size(atti,1);
G = size(atti,2);

% Initialize the partial volume fractions:
f   = ((opt.fmax+opt.fmin)/2)*ones(R,1);      % Rx1
res = zeros([R,2*opt.fbins+1]);               % Rx(2*bins+1)
df  = (opt.fmax-opt.fmin);                    % bins spacing
% Iterate:
if(opt.verbose)
    txtlength = 0;
end
for l=1:opt.fdepth
    df = df/(2*opt.fbins);
    for n=-opt.fbins:1:opt.fbins
        % -- Curent value of f at each voxel:
        fn   = max(   min( f+n*df, opt.fmax ),   opt.fmin   ); % Rx1
        % -- Compute the corrected version of the DWI:
        [adci,impossible] = compute_corrected_atti(atti,fn,bi,opt.tu,opt.tl,opt.ADC0);
        % -- Compute the SH coefficients of the ADC:
        SH = signal2sh( reshape(adci,[R,1,1,G]), gi, ...
            'L', opt.L, 'lambda', opt.lambda, ...
            'chunksz', opt.chunksz );
        % -- Reconstruct the ADC:
        ADC  = sh2signal( SH, gi, 'chunksz', opt.chunksz ); % Rx1x1xG
        ADC  = reshape(ADC,[R,G]);
        ADC  = ADC.*ADC;
        % -- Compute the error and store:
        if(is_broadcast_available)
            ADC = ADC.*reshape(bi,[1,G]); % RxG
            ADC = exp(-ADC);              % RxG
            ADC = ADC.*fn + (1-fn).*exp(-(opt.ADC0)*reshape(bi,[1,G]));
        else
            % bi has size Gx1
            ADC = bsxfun( @(x,y)(x.*y), ADC, reshape(bi,[1,G]) );
            ADC = exp(-ADC);
            ADC = bsxfun( @(x,y)(x.*y), ADC, fn );
            ADC = ADC + exp(-(opt.ADC0))*bsxfun( @(x,y)(x.*y), (1-fn), reshape(bi,[1,G]) );
        end
        resn = sum((ADC-atti).*(ADC-atti),2);
        if(opt.pos)
            resn(impossible) = inf;
        end
        res(:,n+opt.fbins+1) = resn;
        if(opt.verbose)
            fprintf(1,repmat('\b',[1,txtlength]));
            msg = sprintf('depth %d of %d; tested value %d of %d',l,opt.fdepth,n+opt.fbins+1,2*opt.fbins+1);
            txtlength = length(msg);
            fprintf(1,msg);
        end
    end
    [~,minIdx] = min(res,[],2); % Rx1
    f = f + (minIdx-opt.fbins-1)*df;
    f = max(min(f,opt.fmax),opt.fmin);
end
if(opt.verbose)
    fprintf(1,repmat('\b',[1,txtlength]));
end

%%% -----------------------------------------------------------------------
function [f,SH] = atti2freewaterO2( atti, gi, bi, opt, f, SH, mask )
% atti: RxG (R<=M*N*P, all voxels where mask=true)
% gi:   Gx3
% bi:   Gx1
% f:    Rx1
% SH:   RxK
Ym  = (GenerateSHMatrix(opt.L, gi))';
LM  = GenerateSHEigMatrix(opt.L);
it  = 0;
t   = opt.t0*ones(size(f));
lmb = opt.lmb*ones(size(f));
% Make sure f is not exactly 0 or 1 to avoid numerical issues
epsf = 1.0e-3;
f(f<epsf)     = epsf;
f(f>1.0-epsf) = 1.0-epsf;
while( any(mask) && (it<opt.fnriters) ) % While there is at least one voxel to process
    % ---------------------------------------------------------------------
    it    = it+1;
    prc   = find(mask);
    N     = length(prc);
    f0    = f;
    fails = false(size(f));
    % ---------------------------------------------------------------------
    for ck=1:ceil(N/opt.chunksz)
        % -----------------------------------------
        idi = (ck-1)*opt.chunksz+1;
        idf = min(ck*opt.chunksz,N);
        % -----------------------------------------
        [fc,SHc,fail,lm] = newton_raphson_step( ...
            atti(prc(idi:idf),:), ...
            f(prc(idi:idf),1), ...
            SH(prc(idi:idf),:), ...
            bi, Ym, opt.ADC0, ...
            opt.tau, ...
            LM, ...
            t(prc(idi:idf),1), ...
            opt.rcm, ...
            lmb(prc(idi:idf),1), ...
            epsf ...
            ); 
        % -----------------------------------------
        f(prc(idi:idf),:)     = fc;
        SH(prc(idi:idf),:)    = SHc;
        fails(prc(idi:idf),:) = fail;
        lmb(prc(idi:idf),:)   = lm;
        % -----------------------------------------
    end
    % ---------------------------------------------------------------------
    t(~fails) = 0.9*t(~fails);
    mask      = mask & ( (abs(f-f0)>opt.fnrth) | (t>1.0e-3) | fails );
    % ---------------------------------------------------------------------
end

%%% -----------------------------------------------------------------------
function [adci,impossible] = compute_corrected_atti(atti,f,bi,tu,tl,ADC0)
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
G  = size(atti,2);                 % atti is RxG
Ci = exp(-ADC0*reshape(bi,[1,G])); % 1xG

if(is_broadcast_available)
    % f has size Rx1; Ci is 1xG: f.*bi is size RxG as atti
    atti = atti - (1-f).*Ci; % RxG
    atti = atti./f;          % RxG
else
    atti = atti - bsxfun( @(x,y)(x.*y), (1-f), Ci );
    atti = bsxfun( @(x,y)(x./y), atti, f );
end

impossible = any( (atti>tu) | (atti<tl), 2 ); % Rx1
atti(atti>tu) = tu;
atti(atti<tl) = tl;
adci = -log(atti); % RxG

if(is_broadcast_available)
    adci = adci./reshape(bi,[1,G]);
else
    adci = bsxfun( @(x,y)(x./y), adci, reshape(bi,[1,G]) );
end
adci = sqrt(max(adci,0));

%%% -----------------------------------------------------------------------
function [fc,SHc,fail,lmb] = newton_raphson_step( Si, f, SH, bi, Ym, ADC0, tau, LM, t, rcm, lmb, epsf )
% Si:   RxG, meaured signal (R is a chunk, not the whole volume)
% f:    Rx1, current estimation of non-free water
% SH:   RxK, SH coefficients
% bi:   Gx1, b-values
% Ym:   KxG, values of the Spherical Harmonic functions. Ym(k,g) is the
%            value of the k-th SH at gradient direction number g
% ADC0: 1x1, free water diffusitivity
% tau:  1x1, Laplace-Beltrami penalty term for sqrt(ADC)
% LM:   KxK, eigenvalues of SH functions
% t:    Rx1, log-barrier parameter
% rcm:  1x1, minimum allowed condition number
% lmb:  Rx1, Levenberg-Marquardts' parameter
% fc:   Rx1, updated estimation of non-free water
% SHc:  RxK, updated SH coefficients
% fail: Rx1 (logical) failed iterations
% lmb:  Rx1, corrected Levenberg-Marquardts' parameter
% epsf: 1x1, threshold to avoid t getting too close to 0 or 1
Q0   = compute_optimization_cost(Si,f,SH,bi,Ym,ADC0,tau,LM,zeros(size(t))); % Rx1
t    = Q0.*t;
Q0   = compute_optimization_cost(Si,f,SH,bi,Ym,ADC0,tau,LM,t);              % Rx1
g    = compute_optimization_gradient(Si,f,SH,bi,Ym,ADC0,tau,LM,t);          % Rx(K+1)
h    = compute_optimization_hessian(Si,f,SH,bi,Ym,ADC0,tau,LM,t);           % Rx(K+1)x(K+1)
fc   = f;
SHc  = SH;
fail = false(size(Si,1),1);
for r=1:size(Si,1)
    % No other way than computing voxel by voxel, since matrix inversion
    % cannot be programmed vector-wise
    J  = g(r,:)';
    H0 = permute(h(r,:,:),[2,3,1]);
    % Regular Newton-Raphson step should look like:
    %  xn = xp - H^(-1)J
    % but we will use a Levenberg-Marquardt-like term instead:
    %  xn = xp - (H'H+lambda I_N)^(-1)H'J
    H  = H0;                         % H is now positive definite
    sz = sqrt(trace(H*H)/size(H,1)); % 1x1, RMS size of the Hessian
    H  = H0 + lmb(r,1)*sz*eye(size(H,1));
    % H can be singular if Levenberg-Marquardt's parameter is small enough:
    if( isnan(rcond(H)) || isinf(rcond(H)) )
        H  = lmb(r,1)*sz*eye(size(H,1));
    end
    if(rcond(H)<rcm)
        H  = lmb(r,1)*sz*eye(size(H,1));
    end
    delta = -H\J;
    % Make sure f is still in bounds, otherwise reduce the step:
    fest  = f(r,1) + delta(1,1);
    if(fest<=epsf)
        delta = delta*( (f(r,1)-epsf) / (f(r,1)-fest) );
    elseif(fest>=1.0-epsf)
        delta = delta*( (1.0-epsf-f(r,1)) / (fest-f(r,1)) );
    end
    % Update the values:
    fc(r,1)  = f(r,1)  + delta(1,1);
    SHc(r,:) = SH(r,:) + delta(2:end,1)';
end
Q    = compute_optimization_cost(Si,fc,SHc,bi,Ym,ADC0,tau,LM,t); % Rx1
undo = (Q>Q0);
undo = ( undo | isnan(fc) | isinf(fc) | ...
    any(isinf(SHc),2) | any(isnan(SHc),2) );
fc(undo,1)   = f(undo,1);
SHc(undo,:)  = SH(undo,:);
fail(undo,:) = true;
lmb(undo,:)  = lmb(undo,:)*2;
lmb(~undo,:) = lmb(~undo,:)/2;

%%% -----------------------------------------------------------------------
function Q = compute_optimization_cost(Si,f,SH,bi,Ym,ADC0,tau,LM,t,epst)
% Si:   RxG, measured signal
% f:    Rx1, fraction of non-free water
% SH:   RxK, SH coefficients
% bi:   Gx1, b-values
% Ym:   KxG, values of the Spherical Harmonic functions. Ym(k,g) is the
%       value of the k-th SH at gradient direction number g
% ADC0: 1x1, free water diffusitivity
% tau:  1x1, Laplace-Beltrami penalty term for sqrt(ADC)
% LM:   KxK, eigenvalues of SH functions
% t:    Rx1, log-barrier parameter
% Q:    Rx1, the computed cost
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
if(nargin<10)
    epst = eps;
end
ADCi = SH*Ym;      % RxG
ADCi = ADCi.*ADCi; % RxG
if(is_broadcast_available)
    ADCi = ADCi.*(bi'); % RxG
else
    ADCi = bsxfun( @(x,y)(x.*y), ADCi, bi' ); % RxG
end
atti = exp(-ADCi); % RxG
if(is_broadcast_available)
    atti = atti.*f + (1-f).*exp(-(bi')*ADC0); % RxG
else
    atti = bsxfun( @(x,y)(x.*y), atti, f ) + ...
        bsxfun( @(x,y)(x.*y), 1-f, exp(-(bi')*ADC0) ); % RxG
end
barrier = t.*log(f.*(1-f)); % Rx1
barrier(t<epst) = 0;
Q  = sum((atti-Si).*(atti-Si),2)/2 - barrier;  % Rx1
SH = SH*LM;
Q  = Q + (tau/2)*sum(SH.*SH,2);

%%% -----------------------------------------------------------------------
function g = compute_optimization_gradient(Si,f,SH,bi,Ym,ADC0,tau,LM,t)
% Si:   RxG, measured signal
% f:    Rx1, fraction of non-free water
% SH:   RxK, SH coefficients
% bi:   Gx1, b-values
% Ym:   KxG, values of the Spherical Harmonic functions. Ym(k,g) is the
%       value of the k-th SH at gradient direction number g
% ADC0: 1x1, free water diffusitivity
% tau:  1x1, Laplace-Beltrami penalty term for sqrt(ADC)
% LM:   KxK, eigenvalues of SH functions
% t:    Rx1, log-barrier parameter
% g:    Rx(K+1), the computed gradient
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
SHDi = SH*Ym;      % RxG
ADCi = SHDi.*SHDi; % RxG
if(is_broadcast_available)
    ADCi = ADCi.*(bi'); % RxG
else
    ADCi = bsxfun( @(x,y)(x.*y), ADCi, bi' ); % RxG
end
Ai   = exp(-ADCi);       % RxG
A0i  = exp(-(bi')*ADC0); % 1xG
if(is_broadcast_available)
    deltai = f.*Ai + (1-f).*A0i - Si; % RxG
    g      = deltai.*(Ai-A0i);        % RxG
else
    deltai = bsxfun( @(x,y)(x.*y), f, Ai ) + ...
        bsxfun( @(x,y)(x.*y), (1-f), A0i ) - Si;   % RxG
    g      = deltai.*bsxfun( @(x,y)(x-y),Ai,A0i);  % RxG
end
g = sum(g,2) - t.*(1-2*f)./(f.*(1-f)); % Rx1, just the first of K+1 columns
deltai = -2*deltai.*Ai;                % RxG
if(is_broadcast_available)
    deltai = (deltai.*f).*(bi'); % RxG
else
    deltai = bsxfun( @(x,y)(x.*y), bsxfun( @(x,y)(x.*y), deltai, f ), bi' ); % RxG
end
deltai = deltai.*SHDi;
if(is_broadcast_available)
    deltai = deltai.*permute(Ym,[3,2,1]); % RxGxK
else
    deltai = bsxfun( @(x,y)(x.*y), deltai, permute(Ym,[3,2,1]) ); % RxGxK
end
g  = [ g, permute(sum(deltai,2),[1,3,2]) ]; % Rx(K+1)
g2 = tau*SH*(LM*LM);                        % RxK
g(:,2:end) = g(:,2:end) + g2;               % Rx(K+1)

%%% -----------------------------------------------------------------------
function h = compute_optimization_hessian(Si,f,SH,bi,Ym,ADC0,tau,LM,t)
% Si:   RxG, measured signal
% f:    Rx1, fraction of non-free water
% SH:   RxK, SH coefficients
% bi:   Gx1, b-values
% Ym:   KxG, values of the Spherical Harmonic functions. Ym(k,g) is the
%       value of the k-th SH at gradient direction number g
% ADC0: 1x1, free water diffusitivity
% tau:  1x1, Laplace-Beltrami penalty term for sqrt(ADC)
% LM:   KxK, eigenvalues of SH functions
% t:    Rx1, log-barrier parameter
% h:    Rx(K+1)x(K+1), the computed hessian
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
SHDi = SH*Ym;      % RxG
ADCi = SHDi.*SHDi; % RxG
if(is_broadcast_available)
    ADCi = ADCi.*(bi'); % RxG
else
    ADCi = bsxfun( @(x,y)(x.*y), ADCi, bi' ); % RxG
end
Ai   = exp(-ADCi);       % RxG
A0i  = exp(-(bi')*ADC0); % 1xG
if(is_broadcast_available)
    deltai = f.*Ai + (1-f).*A0i - Si; % RxG
    h      = (Ai-A0i);                % RxG
else
    deltai = bsxfun( @(x,y)(x.*y), f, Ai ) + ...
        bsxfun( @(x,y)(x.*y), (1-f), A0i ) - Si;   % RxG
    h      = bsxfun( @(x,y)(x-y), Ai, A0i );       % RxG
end
h = sum(h.*h,2) - t.*( (-1+2*f-2*f.*f)./(f.*f.*(1-f).*(1-f)) ); % Rx1, just the (1,1) element of the Hessian
if(is_broadcast_available)
    h2 = f.*(Ai-A0i); % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), f, bsxfun( @(x,y)(x-y), Ai, A0i ) ); % RxG
end
h2 = -2*(h2+deltai).*Ai; % RxG
if(is_broadcast_available)
    h2 = h2.*(bi'); % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), h2, bi' ); % RxG
end
h2 = h2.*SHDi; % RxG
if(is_broadcast_available)
    h2 = h2.*permute(Ym,[3,2,1]); % RxGxK
else
    h2 = bsxfun( @(x,y)(x.*y), h2, permute(Ym,[3,2,1]) ); % RxGxK
end
h = [ h, permute(sum(h2,2),[1,3,2]) ]; % Rx(K+1), just the first of K+1 slices
if(is_broadcast_available)
    h2 = f.*Ai + deltai;                         % RxG
    h2 = h2.*(bi');                              % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), f, Ai ) + deltai; % RxG
    h2 = bsxfun( @(x,y)(x.*y), h2, bi' );        % RxG
end
h2 = Ai.*(2*(SHDi.*SHDi).*h2-deltai); % RxG
if(is_broadcast_available)
    h2 = h2.*(bi');                        % RxG
    h2 = 2*f.*h2;                          % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), h2, bi' );  % RxG
    h2 = bsxfun( @(x,y)(x.*y), 2*f, h2 );  % RxG
end
Yn = permute(Ym,[3,2,1,4]); % 1xGxKx1
Yp = permute(Ym,[3,2,4,1]); % 1xGx1xK
if(is_broadcast_available)
    h2 = h2.*(Yn.*Yp);  % RxGxKxK
else
    h2 = bsxfun( @(x,y)(x.*y), h2, bsxfun( @(x,y)(x.*y), Yn, Yp ) ); % RxGxKxK
end
h2 = permute(sum(h2,2),[1,3,4,2]); % RxKxK
K  = size(h2,2);
h(:,2:K+1,2:K+1) = h2;                            % Rx(K+1)x(K+1)
h(:,1,2:end)     = permute(h(:,2:end,1),[1,3,2]); % Rx(K+1)x(K+1)
h2 = [zeros(1,size(LM,2)+1);zeros(size(LM,1),1),(LM*LM)];
h2 = tau*h2;                                      % (K+1)x(K+1), sym
h2 = permute(h2,[3,1,2]);                         % 1x(K+1)x(K+1)
if(is_broadcast_available)
    h = h + h2;                                   % Rx(K+1)x(K+1)
else
    h = bsxfun( @(x,y)(x+y), h, h2 );             % Rx(K+1)x(K+1)
end

%%% -----------------------------------------------------------------------
function derivatives_tester %#ok<DEFNU>
% -------------------------------------------------------------------------
% ---
rng(13);
sg = 0.05;
% ---
N  = 51;
Gi = randn(N,3);
Gi = Gi./repmat(sqrt(sum(Gi.*Gi,2)),[1,3]);
bi = rand(N,1)*2000 + 1000;
% ---
ADC0 = 3.0e-3;
% ---
f1  = 0.4;
u11 = [1;1;0]/sqrt(2);   l11 = 1.2e-3;
u12 = [1;-1;0]/sqrt(2);  l12 = 0.4e-3;
u13 = [0;0;1];           l13 = 0.3e-3;
DT1 = l11*(u11*u11') + l12*(u12*u12') + l13*(u13*u13');
% ---
f2  = 0.3;
u21 = [0;1;1]/sqrt(2);   l21 = 1.1e-3;
u22 = [0;-1;1]/sqrt(2);  l22 = 0.3e-3;
u23 = [1;0;0];           l23 = 0.3e-3;
DT2 = l21*(u21*u21') + l22*(u22*u22') + l23*(u23*u23');
% ---
ds1 = Gi(:,1).*Gi(:,1)*DT1(1,1) + Gi(:,2).*Gi(:,2)*DT1(2,2) + ...
    Gi(:,3).*Gi(:,3)*DT1(3,3) + 2*Gi(:,1).*Gi(:,2)*DT1(1,2) + ...
    2*Gi(:,1).*Gi(:,3)*DT1(1,3) + 2*Gi(:,2).*Gi(:,3)*DT1(2,3);
ds2 = Gi(:,1).*Gi(:,1)*DT2(1,1) + Gi(:,2).*Gi(:,2)*DT2(2,2) + ...
    Gi(:,3).*Gi(:,3)*DT2(3,3) + 2*Gi(:,1).*Gi(:,2)*DT2(1,2) + ...
    2*Gi(:,1).*Gi(:,3)*DT2(1,3) + 2*Gi(:,2).*Gi(:,3)*DT2(2,3);
Si = f1*exp(-bi.*ds1) + f2*exp(-bi.*ds2) + (1-f1-f2)*exp(-bi*ADC0);
Si = abs( Si + sg*randn(N,1) + sg*1i*randn(N,1) );
% ---
L      = 2;
lambda = 0.03;
B      = GenerateSHMatrix( L, Gi );
LR     = GenerateSHEigMatrix( L );
WLS    = (B'*B+(lambda).*(LR*LR))\(B');
% ---
Ai = Si/(f1+f2);
Ai = log(Ai);
Ai = -Ai./bi;
Ai(Ai<0) = 0;
Ai = sqrt(Ai);
SH = WLS*Ai;
% ---
f0 = (1-f1-f2)+(rand-0.5)*0.2;
t  = 0.01;
Q  = compute_optimization_cost(Si',f0,SH',bi,B',ADC0,t);
g  = compute_optimization_gradient(Si',f0,SH',bi,B',ADC0,t);
h  = compute_optimization_hessian(Si',f0,SH',bi,B',ADC0,t);
h  = permute(h,[2,3,1]);
% ---
% -------------------------------------------------------------------------
% ---
cf = [f0;SH];
df = abs(cf)*0.001;
% ---
Ain = B*SH;
Ain = exp(-bi.*(Ain.*Ain));
Ain = f0*Ain + (1-f0)*exp(-bi*ADC0);
Ain = Ain - Si;
Qn  = sum(Ain.*Ain)/2 - t*log(f0*(1-f0));
% ---
gn    = zeros(1,size(B,2)+1);
gn(1) = ( ...
    compute_optimization_cost(Si',f0+df(1)/2,SH',bi,B',ADC0,t) - ...
    compute_optimization_cost(Si',f0-df(1)/2,SH',bi,B',ADC0,t) ...
    )/df(1);
for k=1:size(B,2)
    SHn = SH; SHn(k) = SH(k) + df(k+1)/2;
    SHp = SH; SHp(k) = SH(k) - df(k+1)/2;
    gn(k+1) = ( ...
        compute_optimization_cost(Si',f0,SHn',bi,B',ADC0,t) - ...
        compute_optimization_cost(Si',f0,SHp',bi,B',ADC0,t) ...
        )/df(k+1);
end
% ---
hn    = zeros(size(B,2)+1,size(B,2)+1);
hn(1,:) = ( ...
    compute_optimization_gradient(Si',f0+df(1)/2,SH',bi,B',ADC0,t) - ...
    compute_optimization_gradient(Si',f0-df(1)/2,SH',bi,B',ADC0,t) ...
    )/df(1);
for k=1:size(B,2)
    SHn = SH; SHn(k) = SH(k) + df(k+1)/2;
    SHp = SH; SHp(k) = SH(k) - df(k+1)/2;
    hn(k+1,:) = ( ...
        compute_optimization_gradient(Si',f0,SHn',bi,B',ADC0,t) - ...
        compute_optimization_gradient(Si',f0,SHp',bi,B',ADC0,t) ...
        )/df(k+1);
end
% -------------------------------------------------------------------------
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'**** Costs: (%1.7f, %1.7f, %1.4g)\n',Q,Qn,abs(Q-Qn)/Qn);
fprintf(1,'**** Gradients: [Max. rel. err.: %1.4g]\n', max(abs(g-gn)./abs(gn)) );
for k=1:size(B,2)+1
    fprintf('(%1.4f, %1.4f, %1.4g), ', g(k), gn(k), abs(g(k)-gn(k))/abs(gn(k)) );
end
fprintf(1,'\n');
fprintf(1,'**** Hessians: [Max. rel. err.: %1.4f]\n',max(abs(h(:)-hn(:))./abs(hn(:))));
for k=1:size(B,2)+1
    for l=k:size(B,2)+1
        fprintf('(%1.4f, %1.4f, %1.4g), ', h(k,l), hn(l,k), abs(h(k,l)-hn(l,k))/abs(hn(k,l)) );
    end
end
fprintf('\n');
fprintf(1,'-----------------------------------------------------------\n');
keyboard;

% -------------------------------------------------------------------------
function newton_raphson_tester %#ok<DEFNU>
% -------------------------------------------------------------------------
% ---
rng(13);
PSNR = 5;
sg   = 1/PSNR;
% ---
N  = 51;
Gi = randn(N,3);
Gi = Gi./repmat(sqrt(sum(Gi.*Gi,2)),[1,3]);
bi = rand(N,1)*2000 + 1000;
% ---
ADC0 = 3.0e-3;
% ---
f1  = 0.4;
u11 = [1;1;0]/sqrt(2);   l11 = 1.2e-3;
u12 = [1;-1;0]/sqrt(2);  l12 = 0.4e-3;
u13 = [0;0;1];           l13 = 0.3e-3;
DT1 = l11*(u11*u11') + l12*(u12*u12') + l13*(u13*u13');
% ---
f2  = 0.3;
u21 = [0;1;1]/sqrt(2);   l21 = 1.1e-3;
u22 = [0;-1;1]/sqrt(2);  l22 = 0.3e-3;
u23 = [1;0;0];           l23 = 0.3e-3;
DT2 = l21*(u21*u21') + l22*(u22*u22') + l23*(u23*u23');
% ---
ds1 = Gi(:,1).*Gi(:,1)*DT1(1,1) + Gi(:,2).*Gi(:,2)*DT1(2,2) + ...
    Gi(:,3).*Gi(:,3)*DT1(3,3) + 2*Gi(:,1).*Gi(:,2)*DT1(1,2) + ...
    2*Gi(:,1).*Gi(:,3)*DT1(1,3) + 2*Gi(:,2).*Gi(:,3)*DT1(2,3);
ds2 = Gi(:,1).*Gi(:,1)*DT2(1,1) + Gi(:,2).*Gi(:,2)*DT2(2,2) + ...
    Gi(:,3).*Gi(:,3)*DT2(3,3) + 2*Gi(:,1).*Gi(:,2)*DT2(1,2) + ...
    2*Gi(:,1).*Gi(:,3)*DT2(1,3) + 2*Gi(:,2).*Gi(:,3)*DT2(2,3);
Si = f1*exp(-bi.*ds1) + f2*exp(-bi.*ds2) + (1-f1-f2)*exp(-bi*ADC0);
Si = abs( Si + sg*randn(N,1) + sg*1i*randn(N,1) );
% ---
L      = 4;
lambda = 0.03;
B      = GenerateSHMatrix( L, Gi );
LR     = GenerateSHEigMatrix( L );
WLS    = (B'*B+(lambda).*(LR*LR))\(B');
% ---
Ai = Si/(f1+f2);
Ai = log(Ai);
Ai = -Ai./bi;
Ai(Ai<0) = 0;
Ai = sqrt(Ai);
SH = WLS*Ai;
% ---
f0   = (1-f1-f2)+(rand-0.5)*0.2;
fact = f1+f2;
% -------------------------------------------------------------------------
Ym   = (GenerateSHMatrix(L, Gi))';
it   = 0;
t    = 0.01;
rcm  = 1.0e-6;
epsf = 1.0e-3;
lmb  = 0.1;
doit = true;
f    = f0;
SH   = SH';
Si   = Si';
fprintf(1,'\n');
while( doit && (it<1000) )
    % ---------------------------------------------------------------------
    it   = it+1;
    f0   = f;
    % ---------------------------------------------------------------------
    Q0 = compute_optimization_cost(Si,f,SH,bi,Ym,ADC0,t);
    [f,SH,fail,lmb] = newton_raphson_step( ...
        Si, ...
        f, ...
        SH, ...
        bi, Ym, ADC0, ...
        t, ...
        rcm, ...
        lmb, ...
        epsf ...
        );
    Qf = compute_optimization_cost(Si,f,SH,bi,Ym,ADC0,t);
    % ---------------------------------------------------------------------
    doit      = (abs(f-f0)>0.000001) || fail;
    if(~fail)
        t = 0.9*t;
    end
    % ---------------------------------------------------------------------
    if(fail)
        fprintf(1,'it: %d:: Q: %1.4f -> %1.4f (failed, t=%1.4g, lambda=%1.4g), f: %1.4f -> %1.4f\n',it,Q0,Qf,t,lmb,f0,f);
    else
        fprintf(1,'it: %d:: Q: %1.4f -> %1.4f (succeeded, t=%1.4g, lambda=%1.4g), f: %1.4f -> %1.4f\n',it,Q0,Qf,t,lmb,f0,f);
    end
    fc(it)   = f0; %#ok<AGROW>
    fc(it+1) = f;  %#ok<AGROW>
    Qc(it)   = Q0; %#ok<AGROW>
    Qc(it+1) = Qf; %#ok<AGROW>
end
figure(1);
subplot(1,2,1);
hold('on');
plot(1:length(fc),fc,'Color',[0,0,0],'LineWidth',2);
plot([1,length(fc)],[fact,fact],'Color',[0,0,0.5],'LineWidth',2);
grid('on'); xlabel('it'); ylabel('f');
subplot(1,2,2);
plot(1:length(Qc),Qc,'Color',[0,0,0],'LineWidth',2);
grid('on'); xlabel('it'); ylabel('Q');
keyboard;
