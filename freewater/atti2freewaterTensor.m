function [f,tens] = atti2freewaterTensor( atti, gi, bi, varargin )
% function [f,tens] = atti2freewaterTensor( atti, gi, bi, 'opt1', value1, 'opt2', value2, ... )
%
%   Fits the DT-MRI and free water fraction estimated from a
%   DWI model following a dumb exhaustive search/least squares fit
%   procedure similar the one described by Koay et all:
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
%      S_i(gi,bi) = f*S_0*exp(-bi*ADC(gi)) + (1-f)*S_0*exp(-bi*gi^T·D·gi),
%
%   that is further regularized with the additional condition:
%
%      sqrt(tau)*(f-f0) = 0;
%
%   where f is the partial volume fraction of non-free (intra-cellular)
%   water and ADC0 is the diffusivity of free water at body temperature,
%   which has been empirically fixed near 3.0e-3 mm^2/s. gi^T·D·gi is the 
%   DT-MRI non-free water representation. f0 is an initial estimate of the
%   actual value of f, and tau is a Tikhonov regularization parameter to
%   fit the model. The value of f0 are subsequently updated as iterations
%   proceed, so that f0 at iteration n+1 equals f at iteration n.
%
%   Then:
%
%         (S_i/S_0) * (1-f*(S_0/S_i)*exp(-bi*ADC_free)) / (1-f)
%                     = exp(-bi*gi^T·D·gi)
%      -> [-log(dwi*(1-f*exp(-bi*ADC)/dwi)/(1-f))]/bi
%                     = gi^T·D·gi [Linear Least Squares]
%
%
%
%   The two outputs of the function are:
%
%      f: a MxNxP aray in the range (0,1) with the partial volume fraction
%         of non-free (i.e. fiber-confined) water at each imaged voxel.
%      tens: a MxNxPx6 array with the coefficients computed at 
%         each voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   For the representation of the DWI in the basis of Spherical Harmonics:
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      ADC0: estimated diffusivity of free water at body temperature. Do
%         not change this value unless youy have a good reason to do so
%         (default: 3.0e-3).
%      pos: wether (true) or not (false) avoiding estimates of f leading to
%        values of the attenuation signal outside the range [tl,tu]
%        (default: false).
%      flb: wether (true) or not (false) using a dumb estimate of a lower
%        bound of f based on the third eigenvalue of the raw signal
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
%         optimal value found at the previous depth level (default: 3).
%   For the dumb iterative algorithm:
%      O2: wether (true) or not (false) use the additional iterative step
%         to look for f (default: false).
%      fiters: the algorithm alternates between optimizing the SH
%         coefficients in a logarithmic domain for a fixed f and optimizing
%         f in the natural domain for fixed SH coefficients. A max. number
%         of fiters iterations is used (default: 20).
%      fth: iterations stop if the absolute differente between the previous
%         value of f and the current one is below this threshold (default:
%         0.01).
%   For the Newton-Raphson optimization algorithm:
%      O3: wether (true) or not (false) use the additional Newton-Raphson
%          optimization to further fit the signal model described above.
%          Note this step may be very, very slow (default: false).
%      reinit: wether (true) or not (false) re-initialize the estimation
%          provided by O1 in case it seems suspicious (i.e. large mean
%          diffusivity of the tensor compartment) (default: true).
%      fnriters: the maximum number of Newton-Raphson iterations to wait
%          for convergence until premature exit (default: 100).
%      fnrth: iterations stop (i.e. convergence is detected) if the
%          absolute difference between the previous value of f and the
%          current one is below this threshold (default: 1.0e-4).
%      t0: penalty term in the log-barrier method to ensure f lays within
%          the range (0,1) (default: 0.01). This is the intial value, it
%          decreases as the iteration proceed.
%      rcm: minimum allowed condition number for the matrix to be inverted
%          at each Newton-Raphson (or alike) step (default: 1.0e-9).
%      lmb: [UNUSED] Included only for backwards compatibility.
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
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
opt.pos = false;        optchk.pos  = [true,true];     % always 1x1 boolean
opt.flb = false;        optchk.flb = [true,true];      % always 1x1 boolean
% -------------------------------------------------------------------------
opt.O1 = true;          optchk.O1 = [true,true];       % always 1x1 boolean
opt.fmin = 0.01;        optchk.fmin = [true,true];     % always 1x1 double
opt.fmax = 1;           optchk.fmax = [true,true];     % always 1x1 double
opt.fbins = 5;          optchk.fbins = [true,true];    % always 1x1 double
opt.fdepth = 3;         optchk.fdepth = [true,true];   % always 1x1 double
% -------------------------------------------------------------------------
opt.O2 = false;         optchk.O2 = [true,true];       % always 1x1 boolean
opt.fiters = 20;        optchk.fiters = [true,true];   % always 1x1 double
opt.fth = 0.01;         optchk.fth = [true,true];      % always 1x1 double
% -------------------------------------------------------------------------
opt.O3 = false;         optchk.O3 = [true,true];       % always 1x1 boolean
opt.reinit = true;      optchk.reinit = [true,true];   % always 1x1 boolean
opt.fnriters = 100;     optchk.fnriters = [true,true]; % always 1x1 double
opt.fnrth = 1.0e-4;     optchk.fnrth = [true,true];    % always 1x1 double
opt.t0 = 0.01;          optchk.t0 = [true,true];       % always 1x1 double
opt.rcm = 1.0e-9;       optchk.rcm = [true,true];      % always 1x1 double
opt.lmb = 1.0;          optchk.lmb = [true,true];      % always 1x1 double
% -------------------------------------------------------------------------
opt.chunksz = 100;      optchk.chunksz = [true,true];  % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the size as the image field
opt.verbose = false;    optchk.verbose = [true,true];  % always 1x1 boolean
opt.test = false;       optchk.test = [true,true];     % always 1x1 boolean
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
f(mask) = 1;
tens = zeros(M*N*P,6);

% If necessary, pre-compute a tensor fitting of the raw signal:
flb  = -inf(M*N*P,1);
if(opt.flb)
    % This is based on the following paper:
    %        Paul M. Macey, M. Albert Thomas, Luke A. Anderson,
    %        "DTI-based upper limit of voxel free water fraction"
    %        Heliyon 4 (2018), e00700.
    % Doesn't seem to be great deal.
    adci  = compute_corrected_atti( atti(mask,:),f(mask), ...
        bi,opt.tu,opt.tl,opt.ADC0);
    tenstmp = signal2dti( reshape(adci,[size(adci,1),1,1,G]), gi, [], ...
        'wls', false, 'unroll', false, 'chunksz', opt.chunksz ); % M*N*P x 1 x 1 x 6
    [~,~,~,~,~,l3] = dti2spectrum( tenstmp ); % M*N*P x 1 x 1
    flb(mask,:) = min(max(1-l3/opt.ADC0,0),1);
end

if(opt.O1)
    f(mask) = atti2freewaterO1( atti(mask,:), gi, bi, opt );
else
    f(mask) = ((opt.fmax+opt.fmin)/2)*ones(imagesz);
end

if(opt.O2)
    [f,tens] = atti2freewaterO2( atti, gi, bi, opt, f );
else
    adci  = compute_corrected_atti(atti(mask,:),f(mask), ...
        bi,opt.tu,opt.tl,opt.ADC0);
    tenstmp = signal2dti( reshape(adci,[size(adci,1),1,1,G]), gi, [], ...
        'wls', false, 'unroll', false, 'chunksz', opt.chunksz ); % M*N*P x 1 x 1 x 6
    tens(mask,:) = reshape(tenstmp,[size(tenstmp,1),6]); % M*N*P x 6
end
if(opt.O3)
    if(opt.reinit)
        % The following hack is suggested by Hoy et al. in order to avoid model
        % ambiguities for mostly isotropic voxels:
        MD  = (tens(:,1)+tens(:,4)+tens(:,6))/3;
        bad = ( mask & (MD>opt.ADC0/2) );
        f(bad,1)    = 0.5;
        tens(bad,:) = tens(bad,:)/2;
    end
    if(opt.fnriters>0)
        % mum = compute_mum_from_gi(gi); %#ok<NASGU>
        [f,tens] = atti2freewaterO3( atti, gi, bi, opt, f, tens, mask );
    end
end

% Finally, remove non-feasible values of f
bad = ( (f<flb) & mask );
if(any(bad))
    f(f<flb) = flb(f<flb);
    adci  = compute_corrected_atti( atti(bad,:),f(bad), ...
        bi,opt.tu,opt.tl,opt.ADC0);
    tenstmp = signal2dti( reshape(adci,[size(adci,1),1,1,G]), gi, [], ...
        'wls', false, 'unroll', false, 'chunksz', opt.chunksz ); % M*N*P x 1 x 1 x 6
    tens(bad,:) = reshape(tenstmp,[size(tenstmp,1),6]); % M*N*P x 6
end

% Reshape back the outputs to their proper sizes:
tens = reshape(tens,[M,N,P,6]);
f  = reshape(f, [M,N,P]);

%%% -----------------------------------------------------------------------
function mum = compute_mum_from_gi(gi) %#ok<DEFNU>
mum = [ gi(:,1).*gi(:,1), 2*gi(:,1).*gi(:,2), 2*gi(:,1).*gi(:,3), ...
    gi(:,2).*gi(:,2), 2*gi(:,2).*gi(:,3), gi(:,3).*gi(:,3) ]'; % 6xG

%%% -----------------------------------------------------------------------
function f = atti2freewaterO1( atti, gi, bi, opt ) % UPDATED
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
        % -- Current value of f at each voxel:
        fn   = max(   min( f+n*df, opt.fmax ),   opt.fmin   ); % Rx1
        % -- Compute the corrected version of the DWI:
        [adci,impossible] = compute_corrected_atti(atti,fn,bi,opt.tu,opt.tl,opt.ADC0);
        % -- Compute the tensor coefficients of the ADC:
        tensor = signal2dti( reshape(adci,[R,1,1,G]), gi, [], ...
            'wls', false, 'unroll', false, 'chunksz', opt.chunksz );
        % -- Reconstruct the ADC:
        ADC  = dti2signal( tensor, gi, 'chunksz', opt.chunksz ); % Rx1x1xG
        ADC  = reshape(ADC,[R,G]);
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
function [f,SH] = atti2freewaterO2( atti, gi, bi, opt, f ) %#ok<STOUT,INUSD> % UPDATED
error('This mode is deprecated. Please switch the ''O2'' flag off in the main function');

%%% -----------------------------------------------------------------------
function [f2,tensor2] = atti2freewaterO3( atti, gi, bi, opt, f, tensor, mask ) % UPDATED
% atti: RxG (R<=M*N*P, all voxels where mask=true)
% gi:     Gx3
% bi:     Gx1
% f:      Rx1
% tensor: Rx6
mum = [ gi(:,1).*gi(:,1), 2*gi(:,1).*gi(:,2), 2*gi(:,1).*gi(:,3), ...
    gi(:,2).*gi(:,2), 2*gi(:,2).*gi(:,3), gi(:,3).*gi(:,3) ]'; % 6xG
% Make sure f is not exactly 0 or 1 to avoid numerical issues
epsf = 1.0e-3;
f(f<epsf)     = epsf;
f(f>1.0-epsf) = 1.0-epsf;
% ------------------------------------
mask0   = mask;
Q0      = -ones(size(f));
Q       = -ones(size(f));
f2      = f;
tensor2 = tensor;
% ------------------------------------
it  = 0;
while( any(mask) && (it<opt.fnriters) ) % While there is at least one voxel to process
    % ---------------------------------------------------------------------
    it      = it+1;
    prc     = find(mask);
    N       = length(prc);
    f0      = f;
    tensor0 = tensor;
    fails   = false(size(f));
    % ---------------------------------------------------------------------
    for ck=1:ceil(N/opt.chunksz)
        % -----------------------------------------
        idi = (ck-1)*opt.chunksz+1;
        idf = min(ck*opt.chunksz,N);
        % -----------------------------------------
        [fc,tensorc,fail,Q0c,Qc] = newton_raphson_step( ...
            atti(prc(idi:idf),:), ...
            f(prc(idi:idf),1), ...
            tensor(prc(idi:idf),:), ...
            bi, mum, opt.ADC0, ...
            opt.rcm ); 
        % -----------------------------------------
        if(it==1)
            Q0(prc(idi:idf),:) = Q0c;
        end
        Q(prc(idi:idf),:)      = Qc;
        f(prc(idi:idf),:)      = fc;
        tensor(prc(idi:idf),:) = tensorc;
        fails(prc(idi:idf),:)  = fail;
        % -----------------------------------------
    end
    % ---------------------------------------------------------------------
    mask = mask & ( (abs(f-f0)>opt.fnrth) | (max(abs(tensor-tensor0),[],2)>opt.fnrth) ) & (~fails);
    % ---------------------------------------------------------------------
end
mask0            = mask0 & (Q<Q0);
f2(mask0,:)      = f(mask0,:);
tensor2(mask0,:) = tensor(mask0,:);

%%% -----------------------------------------------------------------------
function [adci,impossible] = compute_corrected_atti(atti,f,bi,tu,tl,ADC0) % UPDATED
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

%%% -----------------------------------------------------------------------
function [fc,tensorc,fail,Q0,Q] = newton_raphson_step( Si, f, tensor, bi, mum, ADC0, rcm ) % UPDATED
% Si:      RxG, meaured signal (R is a chunk, not the whole volume)
% f:       Rx1, current estimation of non-free water
% tensor:  Rx6, tensor coefficients
% bi:      Gx1, b-values
% mum:     6xG, values of the tensor products of gradients
% ADC0:    1x1, free water diffusitivity
% rcm:     1x1, minimum allowed condition number
% lmb:     Rx1, Levenberg-Marquardts' parameter
% fc:      Rx1, updated estimation of non-free water
% tensorc: Rx6, updated tensor coefficients
% fail:    Rx1 (logical) failed iterations
% lmb:     Rx1, corrected Levenberg-Marquardts' parameter
Q0      = compute_optimization_cost(Si,f,tensor,bi,mum,ADC0);     % Rx1
g       = compute_optimization_gradient(Si,f,tensor,bi,mum,ADC0); % Rx7
h       = compute_optimization_hessian(Si,f,tensor,bi,mum,ADC0);  % Rx7x7
fc      = f;
tensorc = tensor;
fail    = false(size(Si,1),1);
for r=1:size(Si,1)
    % No other way than computing voxel by voxel, since matrix inversion
    % cannot be programmed vector-wise
    J     = g(r,:)';                   % 7x1
    H     = permute(h(r,:,:),[2,3,1]); % 7x7
    % -----------
    rcc = rcond(H);
    if( isnan(rcc) || isinf(rcc) || (rcc<rcm) )
        fail(r) = true;
        fc(r,1)      = f(r,1);
        tensorc(r,:) = tensor(r,:);
    else
        step  = H\J;
        step  = step';
        % -----------
        % Update the values:
        fc(r,1)      = f(r,1)      - step(1,1);
        tensorc(r,:) = tensor(r,:) - step(1,2:end);
    end
end
Q    = compute_optimization_cost(Si,fc,tensorc,bi,mum,ADC0); % Rx1
fail = ( fail | isnan(Q) | isinf(Q) );
Q(fail)         = Q0(fail);
tensorc(fail,:) = tensor(fail,:);
fc(fail,:)      = f(fail,:);


%%% -----------------------------------------------------------------------
function Q = compute_optimization_cost(Si,f,tensor,bi,mum,ADC0) % UPDATED
% Si:     RxG, measured signal
% f:      Rx1, fraction of non-free water
% tensor: Rx6, tensor coefficients
% bi:     Gx1, b-values
% mum:    6xG, values of the gradients products
% ADC0:   1x1, free water diffusitivity
% Q:      Rx1, the computed cost
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
ADCi = tensor*mum; % RxG
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
Q = sum((atti-Si).*(atti-Si),2)/2;  % Rx1

%%% -----------------------------------------------------------------------
function g = compute_optimization_gradient(Si,f,tensor,bi,mum,ADC0) % UPDATED
% Si:     RxG, measured signal
% f:      Rx1, fraction of non-free water
% tensor: Rx6, tensor coefficients
% bi:     Gx1, b-values
% mum:    6xG, values of the gradients products
% ADC0:   1x1, free water diffusitivity
% g:      Rx7, the computed gradient
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
ADCi = tensor*mum; % RxG
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
g = sum(g,2);          % Rx1, just the first of 7 columns
deltai = -deltai.*Ai;  % RxG
if(is_broadcast_available)
    deltai = (deltai.*f).*(bi'); % RxG
else
    deltai = bsxfun( @(x,y)(x.*y), bsxfun( @(x,y)(x.*y), deltai, f ), bi' ); % RxG
end
if(is_broadcast_available)
    deltai = deltai.*permute(mum,[3,2,1]); % RxGx6
else
    deltai = bsxfun( @(x,y)(x.*y), deltai, permute(mum,[3,2,1]) ); % RxGx6
end
g = [ g, permute(sum(deltai,2),[1,3,2]) ]; % Rx7

%%% -----------------------------------------------------------------------
function h = compute_optimization_hessian(Si,f,tensor,bi,mum,ADC0) % UPDATED
% Si:     RxG, measured signal
% f:      Rx1, fraction of non-free water
% tensor: Rx6, tensor coefficients
% bi:     Gx1, b-values
% mum:    6xG, values of the gradients products
% ADC0:   1x1, free water diffusitivity
% h:      Rx7x7, the computed hessian
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
ADCi = tensor*mum; % RxG
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
h = sum(h.*h,2); % Rx1, just the (1,1) element of the Hessian
if(is_broadcast_available)
    h2 = f.*(Ai-A0i); % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), f, bsxfun( @(x,y)(x-y), Ai, A0i ) ); % RxG
end
h2 = -(h2+deltai).*Ai; % RxG
if(is_broadcast_available)
    h2 = h2.*(bi');    % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), h2, bi' ); % RxG
end
if(is_broadcast_available)
    h2 = h2.*permute(mum,[3,2,1]); % RxGxK
else
    h2 = bsxfun( @(x,y)(x.*y), h2, permute(mum,[3,2,1]) ); % RxGxK
end
h = [ h, permute(sum(h2,2),[1,3,2]) ]; % Rx7, just the first of 7 slices
if(is_broadcast_available)
    h2 = f.*Ai + deltai;                         % RxG
    h2 = h2.*((bi.*bi)');                        % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), f, Ai ) + deltai; % RxG
    h2 = bsxfun( @(x,y)(x.*y), h2, (bi.*bi)' );  % RxG
end
h2 = Ai.*h2; % RxG
if(is_broadcast_available)
    h2 = f.*h2;                          % RxG
else
    h2 = bsxfun( @(x,y)(x.*y), f, h2 );  % RxG
end
mun = permute(mum,[3,2,1,4]); % 1xGx6x1
mup = permute(mum,[3,2,4,1]); % 1xGx1x6
if(is_broadcast_available)
    h2 = h2.*(mun.*mup);  % RxGx6x6
else
    h2 = bsxfun( @(x,y)(x.*y), h2, bsxfun( @(x,y)(x.*y), mun, mup ) ); % RxGx6x6
end
h2 = permute(sum(h2,2),[1,3,4,2]); % Rx6x6
h(:,2:7,2:7) = h2;   % Rx7x7
h(:,1,2:end) = permute(h(:,2:end,1),[1,3,2]);

%%% -----------------------------------------------------------------------
function derivatives_tester %#ok<DEFNU> % UPDATED
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
GR  = [ Gi(:,1).*Gi(:,1), 2*Gi(:,1).*Gi(:,2), 2*Gi(:,1).*Gi(:,3), ...
    Gi(:,2).*Gi(:,2), 2*Gi(:,2).*Gi(:,3), Gi(:,3).*Gi(:,3) ];
WLS = ((GR')*GR)\(GR'); % 6xG
% ---
Ai = Si/(f1+f2);
Ai = log(Ai);
Ai = -Ai./bi;
Ai(Ai<0) = 0;
SH = WLS*Ai;
% ---
f0 = (1-f1-f2)+(rand-0.5)*0.2;
Q  = compute_optimization_cost(Si',f0,SH',bi,GR',ADC0);
g  = compute_optimization_gradient(Si',f0,SH',bi,GR',ADC0);
h  = compute_optimization_hessian(Si',f0,SH',bi,GR',ADC0);
h  = permute(h,[2,3,1]);
% ---
% -------------------------------------------------------------------------
% ---
cf = [f0;SH];
df = abs(cf)*0.001;
% ---
Ain = GR*SH;
Ain = exp(-bi.*Ain);
Ain = f0*Ain + (1-f0)*exp(-bi*ADC0);
Ain = Ain - Si;
Qn  = sum(Ain.*Ain)/2;
% ---
gn    = zeros(1,size(GR,2)+1);
gn(1) = ( ...
    compute_optimization_cost(Si',f0+df(1)/2,SH',bi,GR',ADC0) - ...
    compute_optimization_cost(Si',f0-df(1)/2,SH',bi,GR',ADC0) ...
    )/df(1);
for k=1:size(GR,2)
    SHn = SH; SHn(k) = SH(k) + df(k+1)/2;
    SHp = SH; SHp(k) = SH(k) - df(k+1)/2;
    gn(k+1) = ( ...
        compute_optimization_cost(Si',f0,SHn',bi,GR',ADC0) - ...
        compute_optimization_cost(Si',f0,SHp',bi,GR',ADC0) ...
        )/df(k+1);
end
% ---
hn    = zeros(size(GR,2)+1,size(GR,2)+1);
hn(1,:) = ( ...
    compute_optimization_gradient(Si',f0+df(1)/2,SH',bi,GR',ADC0) - ...
    compute_optimization_gradient(Si',f0-df(1)/2,SH',bi,GR',ADC0) ...
    )/df(1);
for k=1:size(GR,2)
    SHn = SH; SHn(k) = SH(k) + df(k+1)/2;
    SHp = SH; SHp(k) = SH(k) - df(k+1)/2;
    hn(k+1,:) = ( ...
        compute_optimization_gradient(Si',f0,SHn',bi,GR',ADC0) - ...
        compute_optimization_gradient(Si',f0,SHp',bi,GR',ADC0) ...
        )/df(k+1);
end
% -------------------------------------------------------------------------
fprintf(1,'-----------------------------------------------------------\n');
fprintf(1,'**** Costs: (%1.7f, %1.7f, %1.4g)\n',Q,Qn,abs(Q-Qn)/Qn);
fprintf(1,'**** Gradients: [Max. rel. err.: %1.4g]\n', max(abs(g-gn)./abs(gn)) );
for k=1:size(GR,2)+1
    fprintf('(%1.4f, %1.4f, %1.4g), ', g(k), gn(k), abs(g(k)-gn(k))/abs(gn(k)) );
end
fprintf(1,'\n');
fprintf(1,'**** Hessians: [Max. rel. err.: %1.4f]\n',max(abs(h(:)-hn(:))./abs(hn(:))));
for k=1:size(GR,2)+1
    for l=k:size(GR,2)+1
        fprintf('(%1.4f, %1.4f, %1.4g), ', h(k,l), hn(l,k), abs(h(k,l)-hn(l,k))/abs(hn(k,l)) );
    end
end
fprintf('\n');
fprintf(1,'-----------------------------------------------------------\n');

% -------------------------------------------------------------------------
function newton_raphson_tester %#ok<DEFNU>
% -------------------------------------------------------------------------
% ---
rng(19);
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
B   = [ Gi(:,1).*Gi(:,1), 2*Gi(:,1).*Gi(:,2), 2*Gi(:,1).*Gi(:,3), ...
    Gi(:,2).*Gi(:,2), 2*Gi(:,2).*Gi(:,3), Gi(:,3).*Gi(:,3) ];
WLS = ((B')*B)\(B'); % 6xG
% ---
Ai       = Si/(f1+f2);
Ai       = log(Ai);
Ai       = -Ai./bi;
Ai(Ai<0) = 0;
tensor   = WLS*Ai;
% ---
f0   = (f1+f2)+(rand-0.5)*0.2;
fact = f1+f2;
% -------------------------------------------------------------------------
Ym     = B';
it     = 0;
rcm    = 1.0e-9;
doit   = true;
f      = f0;
tensor = tensor';
Si     = Si';
fprintf(1,'\n');
while( doit && (it<1000) )
    % ---------------------------------------------------------------------
    it   = it+1;
    f0   = f;
    % ---------------------------------------------------------------------
    Q0 = compute_optimization_cost(Si,f,tensor,bi,Ym,ADC0);
    [f,tensorc,fail,~,~] = newton_raphson_step( ...
        Si, ...
        f, ...
        tensor, ...
        bi, Ym, ADC0, ...
        rcm );
    Qf = compute_optimization_cost(Si,f,tensorc,bi,Ym,ADC0);
    % ---------------------------------------------------------------------
    doit   = (abs(f-f0)>0.00000001) || (norm(tensorc-tensor)>0.00000001) || fail;
    tensor = tensorc;
    % ---------------------------------------------------------------------
    if(fail)
        fprintf(1,'it: %d:: Q: %1.4f -> %1.4f (failed), f: %1.4f -> %1.4f\n',it,Q0,Qf,f0,f);
        break;
    else
        fprintf(1,'it: %d:: Q: %1.4f -> %1.4f (succeeded), f: %1.4f -> %1.4f\n',it,Q0,Qf,f0,f);
    end
    fc(it)   = f0; %#ok<AGROW>
    fc(it+1) = f;  %#ok<AGROW>
    Qc(it)   = Q0; %#ok<AGROW>
    Qc(it+1) = Qf; %#ok<AGROW>
end
close(figure(1));
figure(1);
subplot(1,2,1);
hold('on');
plot(1:length(fc),fc,'Color',[0,0,0],'LineWidth',2);
plot([1,length(fc)],[fact,fact],'Color',[0,0,0.5],'LineWidth',2);
grid('on'); xlabel('it'); ylabel('f');
subplot(1,2,2);
plot(1:length(Qc),Qc,'Color',[0,0,0],'LineWidth',2);
grid('on'); xlabel('it'); ylabel('Q');
