function [f,lperp,nu,lpar] = atti2freewater( atti, gi, bi, varargin )
% function [f,lperp,nu,lpar] = atti2freewater( atti, gi, bi, 
%                                    'opt1', value1, 'opt2', value2, ... )
% function     [~,~,nu,lpar] = atti2freewater( [], gi, bi, 
%                                    'opt1', value1, 'opt2', value2, ... )
%
%   This is a wrapper function to atti2micro, intended to replicate the
%   method described in:
%
%       Antonio Tristan-Vega; Guillem Paris; Rodrigo de Luis-Garcia; 
%       Santiago Aja-Fernandez. A"ccurate free-water estimation in white 
%       matter from fast diffusion MRI acquisitions using the spherical 
%       means technique". Magnetic Resonance in Medicine 87(2), 
%       pp. 1028–1035. Wiley, 2022.
% 
%   to estimate the so-called Free-Water volume fraction from multi-shell
%   acquisitions (or, at least from "one shell and a half"). It computes
%   the spherical means of the signal at each shell and fits a signal
%   model that mixes up the FW compartment plus a signal coming from a
%   Gaussian, rotation invariant impulse response, for which the main 
%   eigenvalue is fixed (and the two secondary eigenvalues are equal). The
%   computation of the spherical means makes the actual distribution of
%   fibers (i.e. the fiber ODF) inside the voxel irrelvant.
%
%   MANDATORY INPUTS: 
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradient over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%         NOTE: atti may be as well an empty array if query mode is used
%         (i.e. if the function is called just to retrieve the optimal
%         regularization parameters nu and lpar).
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that multi-shell acquisitions are allowed.
%
%   OUPUTS:
%
%      f: a MxNxP aray in the range (0,1) with the partial volume fraction
%         of non-free, i.e. confined, water (consistent with atti2micro and
%         atti2freewaterTensor). If atti=[] is passed, then f is also an 
%         empty array [].
%      lperp: MxNxP, the secondary eigenvalue of the Gaussian kernel that
%         describes non-free water (equals the third one due to rotational
%         invariance). If atti=[] is passed, then f is also an empty 
%         array [].
%      nu: 1x1, the regularization parameter used to promote prolate
%         kernels, i.e. to avoid nearly-spherical convolution kernels that
%         would turn the problem ill-posed. 
%      lpar: 1x1, the fixed value of the main eigenvalue of the 
%         convolution kernel, the same for all voxels in the volume, 
%         usually in the range [2.0e-3,2.1e-3] mm^2/s.
%      NOTE: If the values of nu and lpar are directly passed as optional 
%         arguments, then the same values are returned. If no values, or 
%         negative values, or Infs, or Nans are passed for one of them,
%         then the algorithm is set to query mode and their optimal values
%         are found and returned.
%      NOTE: the optimal value of nu will depend just on the gradients
%         table, i.e. gi and bi. This means that (1) you can use the
%         function setting atti to [] to find these parameters and (2) for 
%         batch processing of a whole database you can use query search 
%         just for the first volume, then fix the returned values of nu and
%         lpar for the remaining ones.
%
%   OPTIONAL arguments may be passed as name/value pairs in the regular
%   Matlab style:
%
%      nu: the regularization parameter used to promote prolate kernels.
%         Defaults to -1, so that query mode is forced.
%      lpar: the fixed value (for all voxels) of the parallel diffusivity.
%         Defaults to 2.1e-3. Set -1, inf or nan to force query mode.
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      ADC0: estimated diffusivity of free water at body temperature. Do
%         not change this value unless you have a good reason to do so
%         (default: 3.0e-3).
%      plot: whether (true) or not (false) plot a graph of
%         "actual-vs-estimated" 1-f (i.e. non-FW volume fraction) based on
%         synthetic data. This may be used to check the consistency of the
%         regularization parameters nu and lambda (default: false).
%
%   Other general options:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% Check the mandatory input argments:
if(nargin<3)
    error('At least the atti volume, the gradients table, and the b-values must be supplied');
end
if(~isempty(atti))
    [M,N,P,G] = size(atti);
else
    M = 0; N=0; P=0;
    if(ismatrix(gi))
        G = size(gi,1);
    else
        G = -1;
    end
end
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
    error('The number of b-values bi must match the number of entries of gi');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.nu = -1;            optchk.nu = [true,true];
opt.lpar = 2.1e-3;      optchk.lpar = [true,true];
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
opt.plot = false;       optchk.plot = [true,true];
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
if( isnan(opt.nu) || isinf(opt.nu) || (opt.nu<0)  )
    opt.nu = -1;
end
if( isnan(opt.lpar) || isinf(opt.lpar) || (opt.lpar<0) )
    opt.lpar = -1;
end
if( (opt.nu<0) || (opt.lpar<0) )
    [nu,lpar] = find_optimal_nu_lpar(gi,bi,opt.ADC0,opt.nu,opt.lpar);
else
    nu   = opt.nu;
    lpar = opt.lpar;
end
if(opt.plot)
    popts.SAMPLING = 0;
    popts.METHOD   = 'micro-fixed';
    popts.Nnb      = 4;
    popts.PSNR     = 18;
    popts.nu       = nu;
    popts.lpar     = lpar;
    popts.gi       = gi;
    popts.bi       = bi;
    test_atti2freewater(popts);
end
if(isempty(atti)) % Only query search
    lperp = [];
    f     = [];
    return;
end
% -------------------------------------------------------------------------
% Resolve with a proper call to atti2micro:
[~,lperp,f] = atti2micro( atti, gi, bi, ...
    'lambda', 0.001, ...
    'mu', 0.00, ...
    'nu', nu, ...
    'bth', 100, ...
    'forcelpar', true, ...
    'lpar', lpar, ...
    'verbose', false, ...
    'usef', true, ...
    'regf', true, ...
    'fmin', 0.001, 'fmax', 1, ...
    'nmax', 100, 'dl', 1.0e-7, 'dC', 1.0e-6, ...
    'nolpwarn', true, ...
    'ADC0', opt.ADC0, ...
    'mask', opt.mask ...
    );
end

%%% -----------------------------------------------------------
function [nu,lpar] = find_optimal_nu_lpar(gi,bi,ADC0,nu,lpar)
% All f-values to be tested:
P   = 20;
f   = linspace(0.5,1.0,P);
% Generate N voxels with random configurations of 2
% crossing fibers contaminated with Rician noise:
N    = 200;
PSNR = 18;
atti = arrange_synthetic_random_signals(N,gi,bi,ADC0,f,PSNR); % NxPx1xG
% Greedy search in nu and lpar if necessary
M = 40;
Q = 5;
if(nu<0)
    nuv = exp( linspace(log(1e-5),log(5),M) );
else
    nuv = nu;
end
if(lpar<0)
    lparv = linspace(1.95e-3,2.15e-3,Q);
else
    lparv = lpar;
end
% Compute the bias of the estimation for each pair
% of values (lpar,nu):
fom = zeros( length(lparv), length(nuv) );
for l=1:length(lparv)
    for n=1:length(nuv)
        fest = estimate_f_config_lnu( atti, gi, bi, ADC0, ...
            lparv(l), nuv(n) ); % N x P
        fom(l,n) = mean(abs(median(fest,1)-f));
    end
end
% Find the minimum fom:
[~,pos] = min(fom(:));
[rw,cl] = ind2sub( size(fom), pos );
nu      = nuv(cl);
lpar    = lparv(rw);
end

%%% -----------------------------------------------------------
function fest = estimate_f_config_lnu( atti, gi, bi, ADC0, lpar, nu )
[~,~,fest] = atti2micro( atti, gi, bi, ...
    'lambda', 0.001, ...
    'mu', 0.00, ...
    'nu', nu, ...
    'bth', 100, ...
    'forcelpar', true, ...
    'lpar', lpar, ...
    'verbose', false, ...
    'usef', true, ...
    'regf', true, ...
    'fmin', 0.001, 'fmax', 1, ...
    'nmax', 100, 'dl', 1.0e-7, 'dC', 1.0e-6, ...
    'nolpwarn', true, ...
    'ADC0', ADC0 ...
    ); % NxPx1x1 -> NxP
end

%%% -----------------------------------------------------------
function atti = arrange_synthetic_random_signals(N,gi,bi,ADC0,f,PSNR)
% Random rotations:
q1  = randn(N,4);
nq1 = sqrt(sum(q1.*q1,2));
q1  = q1./nq1;
ph2 = rand(N,1)*(pi/2);
q2  = [ cos(ph2/2), zeros(N,1), zeros(N,1), sin(ph2/2) ];
q2  = multiply_quatrot(q1,q2);
% Random PVF per compartment:
fc1  = 1/2 + 0.4*(rand(1,N)-1/2); % 1 x N
fc2  = 1 - fc1;                   % 1 x N
% Random eigenvectors:
ii   = [1,0,0];
jj   = [0,1,0];
kk   = [0,0,1];
xx1  = apply_quatrot( q1, ii )';
yy1  = apply_quatrot( q1, jj )';
zz1  = apply_quatrot( q1, kk )';
xx2  = apply_quatrot( q2, ii )';
yy2  = apply_quatrot( q2, jj )';
zz2  = apply_quatrot( q2, kk )';
% Random eigenvalues:
l11 = 1.34e-3 + 0.27e-3*randn(1,N);
l21 = 0.39e-3 + 0.10e-3*randn(1,N);
l31 = 0.23e-3 + 0.08e-3*randn(1,N);
l11 = min( max(l11,0), ADC0 );
l21 = min( max(l21,0), ADC0 );
l31 = min( max(l31,0), ADC0 );
l21 = min(l11,l21);
l31 = min(l21,l31);
% ---
l12 = 1.34e-3 + 0.27e-3*randn(1,N);
l22 = 0.39e-3 + 0.10e-3*randn(1,N);
l32 = 0.23e-3 + 0.08e-3*randn(1,N);
l12 = min( max(l12,0), ADC0 );
l22 = min( max(l22,0), ADC0 );
l32 = min( max(l32,0), ADC0 );
l22 = min(l12,l22);
l32 = min(l22,l32);
% Random diffusion tensors:
DT1 = zeros(6,N); % 6xN
DT1(1,:) = l11.*(xx1(1,:).*xx1(1,:)) + l21.*(yy1(1,:).*yy1(1,:)) + l31.*(zz1(1,:).*zz1(1,:));
DT1(2,:) = l11.*(xx1(1,:).*xx1(2,:)) + l21.*(yy1(1,:).*yy1(2,:)) + l31.*(zz1(1,:).*zz1(2,:));
DT1(3,:) = l11.*(xx1(1,:).*xx1(3,:)) + l21.*(yy1(1,:).*yy1(3,:)) + l31.*(zz1(1,:).*zz1(3,:));
DT1(4,:) = l11.*(xx1(2,:).*xx1(2,:)) + l21.*(yy1(2,:).*yy1(2,:)) + l31.*(zz1(2,:).*zz1(2,:));
DT1(5,:) = l11.*(xx1(2,:).*xx1(3,:)) + l21.*(yy1(2,:).*yy1(3,:)) + l31.*(zz1(2,:).*zz1(3,:));
DT1(6,:) = l11.*(xx1(3,:).*xx1(3,:)) + l21.*(yy1(3,:).*yy1(3,:)) + l31.*(zz1(3,:).*zz1(3,:));
% ---
DT2 = zeros(6,N); % 6xN
DT2(1,:) = l12.*(xx2(1,:).*xx2(1,:)) + l22.*(yy2(1,:).*yy2(1,:)) + l32.*(zz2(1,:).*zz2(1,:));
DT2(2,:) = l12.*(xx2(1,:).*xx2(2,:)) + l22.*(yy2(1,:).*yy2(2,:)) + l32.*(zz2(1,:).*zz2(2,:));
DT2(3,:) = l12.*(xx2(1,:).*xx2(3,:)) + l22.*(yy2(1,:).*yy2(3,:)) + l32.*(zz2(1,:).*zz2(3,:));
DT2(4,:) = l12.*(xx2(2,:).*xx2(2,:)) + l22.*(yy2(2,:).*yy2(2,:)) + l32.*(zz2(2,:).*zz2(2,:));
DT2(5,:) = l12.*(xx2(2,:).*xx2(3,:)) + l22.*(yy2(2,:).*yy2(3,:)) + l32.*(zz2(2,:).*zz2(3,:));
DT2(6,:) = l12.*(xx2(3,:).*xx2(3,:)) + l22.*(yy2(3,:).*yy2(3,:)) + l32.*(zz2(3,:).*zz2(3,:));
% Random quadratic forms:
QF1 = ...
    (gi(:,1).*gi(:,1).*bi(:))*DT1(1,:) + ...
    2*(gi(:,1).*gi(:,2).*bi(:))*DT1(2,:) + ...
    2*(gi(:,1).*gi(:,3).*bi(:))*DT1(3,:) + ...
    (gi(:,2).*gi(:,2).*bi(:))*DT1(4,:) + ...
    2*(gi(:,2).*gi(:,3).*bi(:))*DT1(5,:) + ...
    (gi(:,3).*gi(:,3).*bi(:))*DT1(6,:); % GxN
QF2 = ...
    (gi(:,1).*gi(:,1).*bi(:))*DT2(1,:) + ...
    2*(gi(:,1).*gi(:,2).*bi(:))*DT2(2,:) + ...
    2*(gi(:,1).*gi(:,3).*bi(:))*DT2(3,:) + ...
    (gi(:,2).*gi(:,2).*bi(:))*DT2(4,:) + ...
    2*(gi(:,2).*gi(:,3).*bi(:))*DT2(5,:) + ...
    (gi(:,3).*gi(:,3).*bi(:))*DT2(6,:); % GxN
% Random synthetic signal:
Si1 = exp(-QF1).*fc1; % GxN
Si2 = exp(-QF2).*fc2; % GxN
Si  = Si1 + Si2;
% Rearrange into an attenuation signal:
f    = reshape(f,[1,1,length(f)]); % 1x1xP
atti = f.*Si + (1-f).*exp(-ADC0*bi(:)); % GxNxP + Gx1xP -> GxNxP
atti = permute(atti,[2,3,4,1]); % NxPx1xG
% Add Rician noise:
[N,P,O,G] = size(atti);
ni   = ( randn([N,P,O,G]) + 1i*randn([N,P,O,G]) )/PSNR;
n0   = ( randn([N,P,O,1]) + 1i*randn([N,P,O,1]) )/PSNR;
atti = abs(atti+ni)./abs(1+n0);
end
