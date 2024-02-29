function sh = hydidsi2shodf( eap, dti, Qx, Qy, Qz, lattice, varargin )
% function sh = hydidsi2shodf( eap, dti, Qx, Qy, Qz, lattice, ...
%                                'option1', value1, 'option2, value2, ... )
%
%   Computes the SH coefficients of the ODF from a Cartesian sampling of
%   the Ensemble Average Propagator (EAP) as computed with atti2hydidsi. To
%   that end, the function:
%
%     1- Generates a set of evenly spaced directions in (half) the unit
%        sphere, twice as the number of SH coefficients to fit.
%     2- ANALYTICALLY computes the value of the following integral at each
%        unit direction u:
%
%          int_{q=0}^{infinity} q * (-1/4/pi^2) * Laplacian{E}(q*u) dq
%
%        according to the HYDI-DSI representation (NOTE: this is done by
%        the mex function hydidsiQIntegrals_).
%     3- NUMERICALLY fits the first SH coefficients up to order L to the
%        previously sampled signal.
%     4- ANALYTICALLY computes the Funk-Radon transform of the signal by
%        multiplying its SH coefficients by the FRT eigenvalues.
%
%   Mandatory inputs:
%
%     (As returned by atti2hydidsi):
%
%      eap: a MxNxPxK double array containing the evenly spaced samples of
%         the EAP within the MxNxP field of view.
%      dti: a MxNxPx6 double arrary containing the diffusion tensor
%         estimation at each voxel, which is used to define the actual
%         lattice where the EAP is sampled.
%      Qx, Qy, Qz: all MxNxP, the bandwidths of the EAP along each
%         (transformed) axis.
%
%     (As passed to atti2hydidsi):
%
%      lattice: a 3x1 vector describing the shape (radii) of the lattice
%         where the EAP is estimated.
%
%   Outputs:
%
%      sh: A MxNxPx(L+1)(L+2)/2 double array with the SH coefficients of
%         the ODF at each imaged voxel.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      L: an even integer with the maximum order of the SH to be used
%         (default: 8).
%      ADC0: estimated diffusivity of free water at body temperature 
%         (Diso). Should use the same as in atti2micro (default: 3.0e-3).
%      rescale: whether (true) or not (false) re-normalize the final SH
%         coefficients so that the DC component exactly equals 1/sqrt(4*pi)
%         (so that the ODF exactly integrates to 1) (default: true).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with isotropic ODFs
%         (default: all trues).
%      maxthreads: (only for POSIX systems) the maximum number of threads
%         used by mex functions (default: automatically determined).

% Check the mandatory input arguments:
if(nargin<6)
    error('At lest the eap volume, the dti volume, the bandwidths, and the lattice must be supplied');
end
[M,N,P,K] = size(eap);
assert(isequal(size(dti),[M,N,P,6]),sprintf('The size of dti should be [%i,%i,%i,6]',M,N,P));
assert(isequal(size(Qx),[M,N,P]),sprintf('The size of Qx should be [%i,%i,%i]',M,N,P));
assert(isequal(size(Qy),[M,N,P]),sprintf('The size of Qy should be [%i,%i,%i]',M,N,P));
assert(isequal(size(Qz),[M,N,P]),sprintf('The size of Qz should be [%i,%i,%i]',M,N,P));
assert(isequal(size(lattice),[3,1]),'The size of lattice should be [3,1]');
assert(K==(prod(2*lattice+1)+1)/2,'The last dimension of eap must match the lattice size');
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.L = 8;              optchk.L = [true,true];          % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.rescale = true;     optchk.rescale = [true,true];    % always 1x1 boolean
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
opt.maxthreads = 1e6;   optchk.maxthreads = [true,true]; % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Avoid repeated calls to these functions, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
% Unroll to work comfortably:
eap = reshape( eap, [M*N*P,K] );
dti = reshape( dti, [M*N*P,6] );
Qx  = Qx(:);
Qy  = Qy(:);
Qz  = Qz(:);
eap = eap( opt.mask, : );
dti = dti( opt.mask, : );
Qx  = Qx( opt.mask );
Qy  = Qy( opt.mask );
Qz  = Qz( opt.mask );
% -------------------------------------------------------------------------
% Generate an appropriate set of gradient directions:
Gi = designGradients( (opt.L+1)*(opt.L+2), 'plot', false, 'verbose', false );
% Analytically compute the integrals:
integrals = hydidsiQIntegrals_( eap', dti', Qx', Qy', Qz', ...
    lattice, opt.ADC0, Gi, opt.maxthreads ); % G x N
% -------------------------------------------------------------------------
% Fit these integrals to the SH basis. NOTE the DC component has not to be
% computed, since we know the ODF sum up to 1 -> C_0^0 = 1/sqrt(4*pi)
B  = GenerateSHMatrix(opt.L,Gi); % G x (L+1)(L+2)/2
integrals = (B')*integrals;   % (L+1)(L+2)/2 x N
integrals = (B'*B)\integrals; % (L+1)(L+2)/2 x N
% -------------------------------------------------------------------------
% Analytically compute the FRT:
F = diag(GenerateFRTMatrix(opt.L));
if( is_broadcast_available )
    integrals = F.*integrals; % (L+1)(L+2)/2 x N 
else
    integrals = bsxfun( @(x,y)(x.*y), F, integrals ); % (L+1)(L+2)/2 x N
end
% NOTE: a factor x 0.5 must be applied because we are using the central
% section theorem: computing the integral of E(q) in the transverse plane 
% is the same as computing the integral of P(R) from -inf to -inf along the
% longitudinal direction, but we only need the integral from 0 to inf:
integrals = integrals/2;
if(opt.rescale)
    masses = integrals(1,:);
    if(is_broadcast_available)
        integrals = integrals./masses;
    else
        integrals = bsxfun( @(x,y)(x./y), integrals, masses );
    end
end
% -------------------------------------------------------------------------
% Put things back in place:
LL = (opt.L+1)*(opt.L+2)/2;
sh = zeros( M*N*P, LL );
sh(opt.mask,:) = integrals';
sh = reshape( sh, [M,N,P,LL] );
