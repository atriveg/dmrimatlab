function SHodf = atti2shodf( atti, gi, varargin )
% function SHodf = atti2shodf( atti, gi, 'opt1', val1, 'opt2', val2, ... )
%
%   Takes an attenuation signal atti (DWI signal over the baseline) sampled
%   at a set of gradients gi and computes corresponding Orientation
%   Distribution Functions (ODF), either probabilistic (Jacobian weighted)
%   or non probabilistic (non-weighted) in the basis of spherical harmonics
%   by using (Tikhonov regularized) least squares over a mono-exponential 
%   model. This function basically implementes the ODF estimators described
%   in the following papers:
%
%     - Tristan-Vega, A., C-F. Westin, and S. Aja-Fernandez, "A new 
%     methodology for the estimation of fiber populations in the white 
%     matter of the brain with the Funk?Radon transform", NeuroImage, 
%     vol. 49, no. 2: Elsevier, pp. 1301?1315, 2010.
%
%     - Tristan-Vega, A., S. Aja-Fernandez, and C-F. Westin, "On the 
%     Blurring of the Funk-Radon Transform in Q?Ball Imaging", Medical 
%     Image Computing and Computer-Assisted Intervention?MICCAI 2009: 
%     Springer Berlin Heidelberg, pp. 415?422, 2009.
%
%     - Tristan-Vega, A., C-F. Westin, and S. Aja-Fernandez, "Estimation of
%     fiber orientation probability density functions in high angular 
%     resolution diffusion imaging", NeuroImage, vol. 47, no. 2: Elsevier, 
%     pp. 638?650, 2009.
%
%   which we ask you to cite in case you use this software for your
%   research.
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradient over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%
%      SHodf: a MxNxPx((L+1)(L+2)/2) array with the coefficients of the
%         ODFs computed at each voxel. You can recover the ODF at desired
%         directions by using sh2signal.m.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      type: a string that can be either: 'opdt', 'copdt', 'popdt'
%         (probabilistic, jacobian-weighted estimators), 'qballs',
%         'cqballs', or 'pqballs' (non probabilistic estimators) (default:
%         'opdt').
%      q0: this parmeter is used only for 'cqballs' and 'pqballs'. This is
%         the module of the variable q for the measured shell, the dual of 
%         R in the Fourier domain, defined as q0 = gamma*delta*G/(2*pi),
%         with G the module of the sensitizing gradients, gamma the 
%         gyromagnetic ratio, and delta the gradient duration. Then 
%         b = 4*pi^2*tau*q0^2, with tau the effective diffusion time
%         (default: 35 mm^(-1)).
%      L: an even integer with the maximum order of the SH to be used
%         (default: 6).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem (default 0.006).
%      chunksz: the LLS problem reduces to the product of the dwi signal
%         by an inverse matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the LLS
%         inverse at each step, hence taking advantage of matlab's
%         capabilities (default: 100).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-5, 1-1.0e-5).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% Check the mandatory input argments:
if(nargin<2)
    error('At least the atti volume and the gradients table must be supplied');
end
[M,N,P,G] = size(atti);
NV = M*N*P; % Total number of voxels to be processed
if(~ismatrix(gi))
    error('gi must be a 2-d matlab matrix');
end
if(size(gi,1)~=G)
    error('The number of rows in gi must match the 4-th dimension of dwi');
end
if(size(gi,2)~=3)
    error('The gradients table gi must have size Gx3');
end

% Parse the optional input arguments:
opt.type = 'opdt';      optchk.type = [true,false];   % string with variable size
opt.q0 = 35;            optchk.q0 = [true,true];      % always 1x1 double
opt.L = 6;              optchk.L = [true,true];       % always 1x1 double
opt.lambda = 0.006;     optchk.lambda = [true,true];  % always 1x1 double
opt.chunksz = 100;      optchk.chunksz = [true,true]; % always 1x1 double
opt.tl = 1.0e-5;        optchk.tl = [true,true];      % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];      % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Make sure the dwi lay within the proper range:
atti(atti>opt.tu) = opt.tu;
atti(atti<opt.tl) = opt.tl;


% Compute SH-related matrixes:
B   = GenerateSHMatrix( opt.L, gi );    % GxK, where K=(L+1)(L+2)/2
LR  = GenerateSHEigMatrix( opt.L );     % KxK
FRT = GenerateFRTMatrix( opt.L );       % KxK
WLS = (B'*B+(opt.lambda).*LR.*LR)\(B'); % (KxK)^(-1) * (KxG) -> KxG

% Prepare the data to be processed chunk-by-chunk:
atti = reshape(atti,[NV,G]);  % NVxG
mask = opt.mask(:);           % NVx1
% Mask...
atti   = atti(mask,:);        % PVxG

switch(opt.type)
    case 'opdt'
        a    = -1/(4*pi*pi);
        b    = 0;
        WLS1 = FRT*LR*WLS;
        % ---------------------
        sig2 = log(atti);
        sig2 = 2*sig2.*(3+2*sig2).*atti;
        WLS2 = FRT*WLS;
        % ---------------------
        WLS  = WLS1'; % GxK, for convenience, see loop below
        WLS2 = WLS2'; % GxK, for convenience, see loop below
    case 'copdt'
        a    = -1/(16*pi*pi);
        b    = sqrt(1/(4*pi));
        atti = Ein(-log(atti));
        WLS  = FRT*LR*WLS;
        WLS  = WLS'; % GxK, for convenience, see loop below
    case 'popdt'
        a    = 1/(16*pi*pi);
        b    = sqrt(1/(4*pi));
        atti = log(-log(atti));
        WLS  = FRT*LR*WLS;
        WLS  = WLS'; % GxK, for convenience, see loop below
    case 'qballs'
        a    = 1;
        b    = 0;
        WLS  = FRT*WLS;
        WLS  = WLS'; % GxK, for convenience, see loop below
    case 'cqballs'
        a    = -(1/4)*(opt.q0)*(opt.q0);
        b    = 0;
        atti = (1-atti)./log(atti);
        WLS  = FRT*WLS;
        WLS  = WLS'; % GxK, for convenience, see loop below
    case 'pqballs'
        a    = -(1/4)*(opt.q0)*(opt.q0);
        b    = 0;
        atti = 1./log(atti);
        WLS  = FRT*WLS;
        WLS  = WLS'; % GxK, for convenience, see loop below
    otherwise
        error(['Unrecognized ODF estimator: ',opt.type]);
end

% Now, fit the data chunk-by-chunk via least squares where the mask is true
PV     = size(atti,1);
SHmask = zeros(PV,size(WLS,2)); % NVxK
for ck=1:ceil(PV/opt.chunksz)
    idi = (ck-1)*opt.chunksz+1;
    idf = min(ck*opt.chunksz,PV);
    SHmask(idi:idf,:) = atti(idi:idf,:)*WLS; % (chunksz x G) * (G x K) -> (chunksz x K)
    if(strcmp(opt.type,'opdt'))
        SHmask(idi:idf,:) = SHmask(idi:idf,:) + sig2(idi:idf,:)*WLS2;
    end
end
SHmask      = a*SHmask;
SHmask(:,1) = SHmask(:,1) + b;
% Cast the result to the proper size:
SHodf = zeros(NV,size(WLS,2)); % NVx(L+1)(L+2)/2
SHodf(mask,:) = SHmask;
SHodf = reshape(SHodf,[M,N,P,size(WLS,2)]);
