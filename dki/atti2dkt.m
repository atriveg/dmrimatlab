function [dkt,dti,S0] = atti2dkt( atti, gi, bi, dti, S0, varargin )
% function [dkt,dti,S0] = atti2dkt( atti, gi, bi, dti, S0, ...
%                                   'opt1', value1, 'opt2', value2, ... )
%
%   Given the attenuation signal atti, computed as S_i/S_0, the gradients
%   table gi, and the b-values bi it corresponds to, computes the diffusion
%   tensor (dti) and the kurtosis tensor (dkt) by fitting the attenuation
%   signal in the logarithmic domain with linear least squares, according
%   to the model:
%
%     log(atti) = log(S) - log(S_0)
%       = - b * sum_{i=1,j=1}^3 ni·nj·D_{ij}
%         + b^2/6*MD^2 * sum_{i=1,j=1,k=1,l=1}^3 ni·nj·nk·nl·W_{ijkl}
% 
%   where MD is the mean diffusivity at each voxel computed from D_{ij}.
%   Alternatively, S_0 can also be included as a variable to optimize (note
%   the attenuation, not the diffusion, signal is fed to the function, so
%   thay S_0 is indeed a correction factor close to 1, in the same way as
%   in atti2dti).
%   To fit the model, either ordinary LS, weigthed LS or quadratic
%   programming can be used. In the latter case, we use a method inspired
%   in:
%
%     Tabesh, Jensen, Ardekani, and Helpern (2011): "Estimation of Tensors
%     and Tensor-Derived Measures in Diffusional Kurtosis Imaging".
%     Magnetic Resonance in Medicine 65:823-836,
% 
%   in which Nc evenly distributed orientations are designed and, for each
%   of them, three linear constraints are imposed: 1) non-negativity of the
%   diffusion tensor signal, 2) lower bound for the Kurtosis value, and 3)
%   monotonically decreasing behavior of the attenuation signal.
%
%   NOTE: kurtosis imaging requires the acqusition of at least two 
%   different shells to properly work, so that {gi,bi} must comprise a
%   multi-shell scheme. Besides, note the kurtosis model is usually assumed
%   to hold for b-values up to 3,000 s/mm^2 for in vivo tissues [Lazar et 
%   al., MRM 60(4): 774-781 (2008)].
%
%   ALTERNATIVELY, if the user passes non-empty values of dti and S0, the
%   function will return them untouched, and the kurtosis tensor will be
%   computed by simply fitting the residual signal that the diffusion 
%   tensor model cannot explain (with least squares, without any
%   constraints).
%
%   MANDATORY INPUTS:
%
%      atti: a MxNxPxG double array containing the signal sampled at G
%         directions at each voxel within the MxNxP image frame.
%      gi: a Gx3 matrix with the directions sampled table, each row
%         corresponding to a unit vector.
%      bi: a Gx1 vector with the b-values corresponding to each gradient
%         direction.
%      dti: either MxNxPx6 or empty ([]). In the former case, you can use
%         atti2dti to compute this variable from atti, gi and bi. In the
%         latter (RECOMMENDED), it will be internally computed.
%      S0: either MxNxP or empty ([]). In the former case, you can
%         use atti2dti to compute this variable from atti, gi and bi, and 
%         the attenuation signal will be first corrected with this value. 
%         In the latter, it will be assumed to be 1 unless the function is
%         explicitly asked to estimate it.
%
%   OUTPUTS
%
%      dkt: MxNxPx15, where the last dimension contains the 15 unique
%         components of the 4th order Kurtosis tensor, W_{ijkl}. The
%         ordering of these components is consistent with that used in
%         HOT-related functions. Accordingly, you can evaluate the kurtosis
%         tensor using hot2signal.
%      dti: Can be either:
%         MxNxPx6 if the 'unroll' option is switchwd off
%         MxNxPx3x3 if the 'unroll' option is switched on.
%         In the former case, no duplicates of tensor entries are returned,
%         so that squeeze(tensor(x,y,x,:)) = [D11,D12,D13,D22,D23,D33]'.
%      S0: MxNxP, the normalized baseline (a value near 1), so that you
%         can recover your ESTIMATED baseline dividing by S0.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      estimS0: whether (true) or not (false) include the normalization
%         term for the baseline in the estimation loop (default: false).
%      wls: whether (true) or not (false) using weighted least squares
%         besides of ordinary least squares (default: true).
%      wlsit: in case wls is switched on, the maximum number of iterations
%         for the problem (default: 5).
%      wsc: in case wls is switched on, the minimum weight to be applied in
%         the WLS problem with respect to the maximum one (default: 0.01).
%      rcondth: in case wls is switched on, minimum allowed reciprocal 
%         condition number for matrix inversions (default: 1.0e-6).
%      qp: whether (true) or not (false) using quadratic programming with
%         linear constraints to fit the signal (default: true).
%      nconst: in case qp is switched on, the number of orientations for
%         which linear constraints are imposed. Note 3 constraints are
%         designed for each orientation, so that the actual number of
%         constraints in the QP problem becomes 3*nconst (default: 30).
%      negkrt: in case qp is switched on, whether (true) or not (false)
%         allow negative (i.e in the range [-2,0)) values of the kurtosis.
%         Otherwise, non-negative values are enforced (default: false).
%
%   Other optional arguments:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-5, 1-1.0e-5).
%      bth: this is used as a threshold below which two b-values are
%         considered virtually identical. It is internally used to
%         determine if {gi,bi} comprise an actual multi-shell sampling
%         (default: 50).
%      unroll: wether (true) or not (false) output the dti as a 3x3
%         matrix at each voxel instead of a 6x1 vector (in the latter case
%         duplicates of the entries of the diffusion tensor are removed so
%         that only [D11,D12,D13,D22,D23,D33] are returned) (default:
%         false).
%      maxthreads: the algorithm is run with multiple threads. This is the
%         maximum allowed number of threads, which can indeed be reduced if
%         it exceeds the number of logical cores (default: the number of
%         logical cores in the machine).

%%% -----------------------------------------------------------------------
% Check the mandatory input argments:
if(nargin<5)
    error('At lest the atti, the gradients table, the b-values, dti and S0 must be supplied');
end
[M,N,P,G] = size(atti);
assert( isequal(size(gi),[G,3] ), 'The b-vectors gi must have size Gx3, with G the 4-th dimension of atti' );
assert( isequal(size(bi),[G,1] ), 'The b-values bi must have size Gx1, with G the 4-th dimension of atti' );
if(~isempty(dti))
    nd = ndims(dti);
    switch(nd)
        case 4
            assert( size(dti,4)==6, 'The dti volume provided does not look like a diffusion tensor volume' );
        case 5
            assert( (size(dti,4)==3) && (size(dti,5)==3), 'The dti volume provided does not look like a diffusion tensor volume' );
            dti = cat( 4, dti(:,:,:,1,1), dti(:,:,:,1,2), dti(:,:,:,1,3), ...
                dti(:,:,:,2,2), dti(:,:,:,2,3), dti(:,:,:,3,3) );
        otherwise
            error('The dti volume provided does not look like a diffusion tensor volume');
    end
    [M2,N2,P2,~,~] = size(dti);
    assert( isequal([M,N,P],[M2,N2,P2]), 'The FoV of dti does not match the FoV of atti' );
end
if(~isempty(S0))
    assert( isequal([M,N,P],size(S0)), 'The FoV of S0 does not match the FoV of atti' );
end
%%% -----------------------------------------------------------------------
% Parse the optional input arguments:
opt.estimS0 = false;    optchk.estimS0 = [true,true];    % always 1x1 boolean
opt.wls = true;         optchk.wls = [true,true];        % always 1x1 boolean
opt.wlsit = 5;          optchk.wlsit = [true,true];      % always 1x1 double
opt.wsc = 0.01;         optchk.wsc = [true,true];        % always 1x1 double
opt.rcondth = 1.0e-6;   optchk.rcondth = [true,true];    % always 1x1 double
opt.qp = true;          optchk.qp = [true,true];         % always 1x1 boolean
opt.nconst = 30;        optchk.nconst = [true,true];     % always 1x1 double
opt.negkrt = false;     optchk.negkrt = [true,true];     % always 1x1 boolean
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean with the size of the image field
opt.tl = 1.0e-5;        optchk.tl = [true,true];         % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];         % always 1x1 double
opt.bth = 50;           optchk.bth = [true,true];        % always 1x1 double
opt.unroll = false;     optchk.unroll = [true,true];     % always 1x1 boolean
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true]; % always 1x1 char
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%% -----------------------------------------------------------------------
[~,~,Ns] = auto_detect_shells(bi,opt.bth);
assert( Ns>=2, 'Kurtosis imaging requires multi-shell data, but the bi you provided suggests you have single-shell');
%%% -----------------------------------------------------------------------
% Make sure the atti lay within the proper range:
atti(atti>opt.tu) = opt.tu;
atti(atti<opt.tl) = opt.tl;
% This is just an estimate of the range for which the Kurtosis model
% holds, see Lazar et al., MRM 60(4): 774-781 (2008):
bmax = 3000;
% Proceed with computations depending on the input information provided
if( isempty(dti) )
    % ---------------------------------------------------------------------
    % Regardles of the value of S0, use the true joint optimization
    % procedure
    % ---------------------------------------------------------------------
    % Set-up the problem
    atti = reshape(atti,[M*N*P,G]);
    atti = atti( opt.mask(:), : );
    if( (~isempty(S0)) && (~opt.estim0) )
        S0   = reshape(S0,[M*N*P,1]);
        S0   = S0(opt.mask,:);
        atti = atti./S0;
    end
    if(opt.estimS0)
        options.estimS0 = 'y';
    else
        options.estimS0 = 'n';
    end
    options.wlsit = opt.wlsit;
    options.wsc = opt.wsc;
    options.rcondth = opt.rcondth;
    if(opt.qp)
        options.mode = 'q';
    elseif(opt.wls)
        options.mode = 'w';
    else
        options.mode = 'o';
    end
    if(opt.negkrt)
        options.negkrt = 'y';
    else
        options.negkrt = 'n';
    end
    % ---------------------------------------------------------------------
    % Build the inequality constraints:
    gic = designGradients( opt.nconst, 'plot', false, 'verbose', false );
    % ------------------
    % Kurtosis-related:
    u   = eye(15);
    u   = reshape(u,[1,1,15,15]); % 1 x 1 x 15 x 15
    AKc = hot2signal( u, gic );   % 1 x 1 x 15 x G
    AKc = permute(AKc,[4,3,1,2]); % G x 15;
    % ------------------
    % Diffusion tensor-related:
    ADc = [ gic(:,1).*gic(:,1), 2*gic(:,1).*gic(:,2), 2*gic(:,1).*gic(:,3), ...
        gic(:,2).*gic(:,2), 2*gic(:,2).*gic(:,3), gic(:,3).*gic(:,3) ];
    ADc = 1000*ADc; % The 1000* factor is necessary for inner consistency with the mex
    options.Drec = ADc;
    % ------------------
    Ain = [ zeros(opt.nconst,15),                -ADc;    % positive diffusion
                            -AKc, zeros(opt.nconst,6);    % positive kurtosis
                             AKc,        -(3/bmax)*ADc ]; % decaying attenuation signal
    % Add the column to account for the log(S0) unknown, if necessary:
    if(opt.estimS0)
        Ain = [ Ain, zeros(3*opt.nconst,1) ];
    end
    % ------------------
    options.At = Ain';
    % ---------------------------------------------------------------------
    % Solve the problem:
    [dkt_,dti_,S0_,~] = atti2dkt_( double(atti'), double(gi), double(bi), options, opt.maxthreads );
    % ---------------------------------------------------------------------
    % Reshape output values: 
    dkt = zeros(M*N*P,15);
    dti = zeros(M*N*P,6);
    S0  = zeros(M*N*P,1);
    dkt(opt.mask,:) = dkt_';
    dti(opt.mask,:) = dti_';
    S0(opt.mask,:)  = S0_';
    dkt = reshape(dkt,[M,N,P,15]);
    dti = reshape(dti,[M,N,P,6]);
    S0  = reshape(S0,[M,N,P]);
else
    % The user has explicitly provided dti, so that we will fit just the
    % residual of the signal to the 4-th order tensor representing the
    % kurtosis.
    if( isempty(S0) )
        % Use just ones:
        S0 = ones(M,N,P);
    end
    % From the DTI volume, we can just recover the signal:
    sig = dti2signal( dti, gi, 'mask', opt.mask ); % M x N x P x G
    % So that the residual of the tensor model within the log domain is:
    res = log(atti) - log(S0) + sig.*reshape(bi,[1,1,1,G]); % M x N x P x G
    % Normalize by the (squared) b-values (times 6), so that we end up
    % with the kurtosis signal times the squared mean diffusivity:
    res = 6*res./reshape(bi.*bi,[1,1,1,G]);
    % Make sure this signal is in range:
    res(res<0) = 0; % positive kurtosis, from experimental knowledge
    rth = max( 0, sig.*3/bmax );  % M x N x P x G
    res(res>rth) = rth(res>rth); % decreasing atti, mathematical constraint
    % This residual can be fitted to the basis of SH up to order 4:
    dkt = signal2sh( res, gi, 'L', 4, 'mask', opt.mask, 'lambda', 0.0 );
    % And, from SH, we can retrieve the rank-4 kurtosis tensor:
    dkt = sh2hot( dkt, 'mask', opt.mask, 'maxthreads', opt.maxthreads );
    % It only remains to normalize by the mean diffusivity
    MD  = sum( dti(:,:,:,[1,4,6]), 4 )/3;
    dkt = dkt./(MD.*MD);
end

if( (nargout>1) && opt.unroll )
    dti = dti(:,:,:,[1,2,3,2,4,5,3,5,6]);
    dti = reshape(dti,[M,N,P,3,3]);
end

end
