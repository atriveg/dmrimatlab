function pa = atti2apa( atti, gi, bi, varargin )
% function pa = atti2apa( atti, gi, bi, 'option1', value1, ... )
%
%   Computes the (apparent) Propagator Anisotropy (PA) according to a 
%   mono-exponential model for single-shell acquisitions as described in
%   AMURA:
%
%      atti(u,b) = exp(-b*D_0(u)),
%
%   where D_0(u) is the Apparent Diffusion Coefficient measured at b=b0 for
%   each direction 'u' within the unit sphere. The PA is defined in terms
%   of the 'angle' between the true attenuation signal and its isotropic
%   counterpart, which in AMURA is modeled as that with a homogeneous
%   Diffusion Coefficente: D_iso = mean{D_0(u)}.
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
%      pa: MxNxP double array with the computed PA.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).
%      epsilon: to improve the contrast of the raw PA, a gamma correction
%         is usually performed with the form:
%              pa = pa0^(3*epsilon)/(1-3*pa0^epsilon+3*pa0^(2*epsilon));
%         Use empty brackets, [], to avoid this correction and work with
%         the raw PA (default: 0.4).
%
%   Parameters related to SH computations:
%
%      L: the maximim order for the SH expansions used to represent
%         spherical functions (default: 8).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem of fitting SH coefficients (note the order
%         L used for SH is internally computed) (default: 0.001).
%
%   Sanity checks on the attenuation signal:
%
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%

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
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
opt.epsilon = 0.4;      optchk.epsilon = [true,false]; % might be an empty array
% -------------------------------------------------------------------------
opt.L = 8;              optchk.L =[true,true];         % always 1x1 double
opt.lambda = 0.001;     optchk.lambda = [true,true];   % always 1x1 double
% -------------------------------------------------------------------------
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
assert(abs(opt.L-floor(opt.L))<100*eps,'Optional argument ''L'' must be integer');
assert(abs(opt.L/2-floor(opt.L/2))<100*eps,'Optional argument ''L'' must be even');
% -------------------------------------------------------------------------
% Check optional input argument epsilon:
if( ~isempty(opt.epsilon) )
    assert(numel(opt.epsilon)==1,'Optional argument ''epsilon'' must be either an empty array or have size 1x1');
end
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably:
mask  = reshape( opt.mask, [M*N*P,1] ); % M*N*P x 1
atti  = reshape( atti, [M*N*P,G] );     % M*N*P x G
atti  = atti(mask,:);                   % Q x G
% -------------------------------------------------------------------------
% Compute the ADC from atti (keep the same variable to save memory)
is_broadcast_available = is_broadcast_available_test;
atti(atti<opt.tl) = opt.tl;
atti(atti>opt.tu) = opt.tu;
atti = -log(atti);
if(is_broadcast_available)
    atti = atti./(bi');                       % Q x G
else
    atti = bsxfun( @(x,y)(x./y), atti, bi' ); % Q x G
end
% -------------------------------------------------------------------------
% Compute SH-related matrices:
B  = GenerateSHMatrix( opt.L, gi );    % GxK, where K=(L+1)(L+2)/2
LR = GenerateSHEigMatrix( opt.L );     % KxK
B  = (B'*B+(opt.lambda).*LR.*LR)\(B'); % (KxK)^(-1) * (KxG) -> KxG
B  = B(1,:);                           % 1 x G
% -------------------------------------------------------------------------
% Compute the spherical averages
%   D average
if(is_broadcast_available)
    Dav = sum( atti.*B, 2 );                 % Q x 1
else
    Dav = sum( ...
        bsxfun( @(x,y)(x.*y), atti, B ), ... % Q x G
        2 );                                 % Q x 1
end
Dav = Dav/sqrt(4*pi);
%   INT (D+DAV)^(-3/2)
if(is_broadcast_available)
    CP = sum( B.*((Dav+atti).^(-3/2)), 2 );         % Q x 1
else
    CP = sum( ...
        bsxfun( @(x,y)(x.*y), B, ...                % Q x G
        bsxfun(@(x,y)(x+y),Dav,atti).^(-3/2) ), ... % Q x G
        2 );                                        % Q x 1
end
CP = CP.*CP;
%   INT (D)^(-3/2)
if(is_broadcast_available)
    Dnr = sum( B.*(atti.^(-3/2)), 2 );               % Q x 1
else
    Dnr = sum( ...
        bsxfun( @(x,y)(x.*y), B, atti.^(-3/2) ), ... % Q x G
        2 );                                         % Q x 1
end
% -------------------------------------------------------------------------
% Compute the squared cosine:
CS2 = (4/sqrt(pi))*CP./(Dnr.*Dav.^(-3/2));
SN  = sqrt( max(min(1-CS2,1),0) );
% -------------------------------------------------------------------------
if(~isempty(opt.epsilon))
    SN = SN.^(3*opt.epsilon)./ ...
        ( 1 - 3*SN.^(opt.epsilon) + 3*SN.^(2*opt.epsilon));
end
% -------------------------------------------------------------------------
% Assign the output
pa = zeros(M*N*P,1);
pa(mask,:) = SN;
pa = reshape(pa,[M,N,P]);
% -------------------------------------------------------------------------
