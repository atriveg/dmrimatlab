function [eap,dti,lattice,Qx,Qy,Qz] = atti2wuhydidsi( atti, gi, bi, varargin )
% function [eap,dti,lattice,Qx,Qy,Qz] = atti2wuhydidsi( atti, gi, bi, 
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Computes the Ensemble Average Propagator (eap) at a regular Cartesian
%   lattice from a set of scattered q-space measurements (atti) using
%   Hybrid Diffusion Imaging using q-space re-gridding. To that end, the
%   R-space is first aligned with the eigenvectors of a tensor model (the z
%   axis aligns with the maximum diffusion direction) and an appropriate
%   lattice is designed so that the whole extent of E(q) is covered. From 
%   this estimate, the q-space is interpolated using Matlab's built-in 
%   'griddatan', so that computing the EAP reduces to calculate the 3-D
%   FFT.
%
%   Mandatory inputs:
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradients over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%         These are the scattered q-space samples.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that a multi-shell acquisition is arranged. NOTE: at least 6
%         evenly spaced gradients with bi < 2000 should be present so that
%         the tensor model can be reliably estimated.
%
%   The outputs of the function area as follows:
%
%      eap: a MxNxPxK array with (half) the unrolled Cartesian samples of
%         the eap. If the lattice has radii Nx x Ny x Nz, a grand total of
%         Nl=(2*Nx+1)(2*Ny+1)*(2*Nz+1) samples would be needed. However,
%         due to the antipodal symmetry of the eap we need just K=(Nl+1)/2.
%         The samples are arranged as follows:
%            >> [x,y,z] = meshgrid(-Nx:Nx,-Ny:Ny,-Nz:Nz);
%            >> x = x(:); y = y(:); z = z(:);
%            >> Nl = size(x,1);
%            >> x = x(floor(Nl/2):end);
%            >> y = y(floor(Nl/2):end);
%            >> z = z(floor(Nl/2):end);
%      dti: MxNxPx6, a double array with the estimated tensor model  at
%         each voxel. Note this is necessary to interpret the data in eap,
%         since the lattice is aligned at each voxel to the eigenvectors of
%         dti (and perhaps scaled according to its eigenvalues).
%      lattice: 3x1, a double array with the radii of the lattice at each
%         voxel, i.e. [Nx;Ny;Nz]. This is useful in case this parameter is
%         automatically determined instead of fed as an input argument.
%      Qx, Qy, Qz: MxNxP each, double arrays with the estimated bandwidths
%         of the q-space signal along each dimension. They are necessary to
%         interpret the computed eap.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   General optional parameters:
%
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      fixedq: wether (true) or not (false) using a fixed bandwidth for the
%         EAP at all imaged voxels. If set true, then the bandwidth is
%         fixed to the maximum b-value acquired, and the next three
%         optional parameters are not used (default: false).
%      ADC0: estimated diffusivity of free water at body temperature. It is
%         used to determine the lower and upper bounds of the eigenvalues
%         of the dti and perform sanity checks (default: 3.0e-3).
%      tau: the effective diffusion time of the acquisition, necessary to
%         determine the actual support of the q-space  (default: 35e-3).
%      Rth: the threshold used to estimate the support of the EAP. The last
%         lattice node along each `x'-`y'-`z' axis will be placed in the
%         position where the DTI-EAP falls to Rth times its value at the 
%         origin, which is trivially computed from the eigenvalues l1, l2,
%         and l3 (default: 0.01).
%      lattice: the radii of the Cartesian grid where the eap will be
%         sampled, whose size will be typically in the range 3 x 3 x 3. You
%         may pass an empty array [] to let the function determine the
%         proper size. In any case, the third component should be greater
%         or equal than the first and the second (default: []).
%      const: wether (true) or not (false) re-normalize the EAP computed
%         via 3-D FFT so that it is positive and sum up to 1. Hint: set
%         this flag false if fixedq is true (default: true).
%
%   Other general options:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).
%      verbose: wether (true) or not (false) show a progress bar to check
%         the progress of the algorithm (default: false).

% Check the mandatory input argments:
if(nargin<3)
    error('At lest the atti volume, the gradient table, and the b-values must be supplied');
end
[M,N,P,G] = size(atti);
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
    error('The number of b-values bi must match the 4-th dimension of atti');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
opt.fixedq = true;      optchk.fixedq = [true,true];   % always 1x1 boolean
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
opt.tau = 35.0e-3;      optchk.tau = [true,true];      % always 1x1 double
opt.Rth = 0.01;         optchk.Rth = [true,true];      % always 1x1 double
opt.lattice = [];       optchk.lattice = [true,false]; % variable size
opt.const = false;      optchk.const = [true,true];    % always 1x1 boolean
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
opt.verbose = false;    optchk.verbose = [true,true];  % always 1x1 boolean
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
if(isempty(opt.lattice))
    % Auto determine the lattice based on the number of q-space samples
    lt = floor((nthroot(2*G-1,3)-1)/2);
    lt = [lt;lt;lt];
    while( (prod(2*lt+1)+1)/2 <= G )
        opt.lattice = lt;
        lt(3) = lt(3) + 1;
    end
else
    assert(isequal(size(opt.lattice),[3,1]),'lattice should be either 3x1 or empty');
end
assert( (opt.Rth<1)&&(opt.Rth>0), 'Rth must be in the range (0,1), and should be close to 0' );
lattice = opt.lattice; % It will be returned by the function
Nl      = prod(2*lattice+1);
K       = (Nl+1)/2;
% -------------------------------------------------------------------------
% Sanity checks and helper variables:
atti(atti<opt.tl) = opt.tl;
atti(atti>opt.tu) = opt.tu;
pdti = (bi<=2000);
Ndti = length(find(pdti));
assert(Ndti>=6,'At least 6 b-values must be less or equal than 2000 to proceed with DTI estimation');
is_broadcast_available = is_broadcast_available_test;
use_parallel           = use_parallel_test;
% -------------------------------------------------------------------------
% Estimate the tensor model:
dti = atti(:,:,:,pdti); % M x N x P x Ndti
dti = -log(dti);        % M x N x P x Ndti
if(is_broadcast_available)
    dti = dti./reshape(bi(pdti),[1,1,1,Ndti]);
else
    dti = bsxfun(@(x,y)(x./y),dti,reshape(bi(pdti),[1,1,1,Ndti]));
end
dti = signal2dti( dti, gi(pdti,:), bi(pdti), 'mask', opt.mask, ...
    'wls', false, 'unroll', false, 'chunksz', 256 );
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask = reshape(opt.mask,[M*N*P,1]);
atti = reshape(atti,[M*N*P,G]);
atti = atti(mask,:);
dti  = reshape(dti,[M*N*P,6]);
dti  = dti(mask,:);
Q    = size(atti,1);
eapM = zeros(Q,K);
QxM  = zeros(Q,1);
QyM  = zeros(Q,1);
QzM  = zeros(Q,1);
% -------------------------------------------------------------------------
% No other way than looping voxel-to-voxel. For the parallel
% implementation, extract parameters from the opt structure to reduce
% memory exchange:
verbose = opt.verbose;
fixedq  = opt.fixedq;
ADC0    = opt.ADC0;
tau     = opt.tau;
lRth    = -4*tau*log(opt.Rth);
const   = opt.const;
qi      = sqrt(bi/(4*pi*pi*tau));

NVs   = 50;
Qchnk = round(Q/NVs);
if(verbose)
    str   = repmat('-',[1,NVs]);
    fprintf(1,'Completed: [%s]',str);
end
if(use_parallel) % Parallel implementation
    parfor q=1:Q
        [eap,Qx,Qy,Qz] = one_voxel_eap( atti(q,:), dti(q,:), gi, qi, ...
            lattice, ADC0, lRth, const, fixedq );
        eapM(q,:) = eap;
        QxM(q)    = Qx;
        QyM(q)    = Qy;
        QzM(q)    = Qz;
        if(verbose)
            if(rem(q,Qchnk)==0)
                fprintf(1,repmat('\b',[1,NVs+2]));
                nvs = floor(q/Qchnk);
                str = [ repmat('*',[1,nvs]), repmat('-',[1,NVs-nvs]) ];
                fprintf(1,'[%s]',str);
            end
        end
    end
else             % Non-parallel implementation
    for q=1:Q
        [eap,Qx,Qy,Qz] = one_voxel_eap( atti(q,:), dti(q,:), gi, qi, ...
            lattice, ADC0, lRth, const, fixedq );
        eapM(q,:) = eap;
        QxM(q)    = Qx;
        QyM(q)    = Qy;
        QzM(q)    = Qz;
        if(verbose)
            if(rem(q,Qchnk)==0)
                fprintf(1,repmat('\b',[1,NVs+2]));
                nvs = floor(q/Qchnk);
                str = [ repmat('*',[1,nvs]), repmat('-',[1,NVs-nvs]) ];
                fprintf(1,'[%s]',str);
            end
        end
    end
end
if(verbose)
    fprintf(1,'\n');
end
% -------------------------------------------------------------------------
% Reshape back to output the result:
eap  = zeros(M*N*P,K);
Qx   = zeros(M*N*P,1);
Qy   = zeros(M*N*P,1);
Qz   = zeros(M*N*P,1);
eap(mask,:) = eapM;
Qx(mask)    = QxM;
Qy(mask)    = QyM;
Qz(mask)    = QzM;
eap  = reshape(eap,[M,N,P,K]);
Qx   = reshape(Qx,[M,N,P]);
Qy   = reshape(Qy,[M,N,P]);
Qz   = reshape(Qz,[M,N,P]);
% -------------------------------------------------------------------------

function [eap,Qx,Qy,Qz] = one_voxel_eap( atti, dti, gi, qi, lattice, ...
    ADC0, lRth, const, fixedq )
% The function that does it all:
% atti:    1 x G
% dti:     1 x 6
% gi:      G x 3
% bi:      G x 1
% lattice: 3 x 1 -> K = (prod(2*lattice+1)+1)/2
% qcrop:   true/false
% ADC0:    1 x 1
% eap:     1 x K
% -------------------------------------------------------------------------
% Check the tensor model
D     = dti([1,2,3;2,4,5;3,5,6]);
[U,L] = eig(D);
L     = abs(diag(L));
L(L>ADC0)    = ADC0;
L(L<ADC0/60) = ADC0/60;
[L,id]       = sort(L);
y     = U(:,id(2));
z     = U(:,id(3));
x     = cross(y,z);
% -------------------------------------------------------------------------
% Determine the R- and q-domains support. This depends on wether or not the
% option qcrop is set.
gi = gi*[x,y,z];
if(fixedq)
    Qx = 2*max(qi(:));
    Qy = 2*max(qi(:));
    Qz = 2*max(qi(:));
else
    Qx = lattice(1)/sqrt(L(1)*lRth);
    Qy = lattice(2)/sqrt(L(2)*lRth);
    Qz = lattice(3)/sqrt(L(3)*lRth);
end
Q  = Qx*Qy*Qz;
% -------------------------------------------------------------------------
% Create the q-space grid where data will be interpolated:
qox = [(0:1:lattice(1)-1),(-lattice(1):1:-1)]*Qx/lattice(1)/2;
qoy = [(0:1:lattice(2)-1),(-lattice(2):1:-1)]*Qy/lattice(2)/2;
qoz = [(0:1:lattice(3)-1),(-lattice(3):1:-1)]*Qz/lattice(3)/2;
[qox,qoy,qoz] = meshgrid(qox,qoy,qoz);
% -------------------------------------------------------------------------
% Prepare input data for interpolation. NOTE: Besides the actually measured
% data, we need to include:
%   - Their antipodes, associated to the same atti values
%   - The origin, to ensure the EAP has integral 1
%   - The bounds of the q-space, to ensure the convex hull covers all the
%     grid to interpolate.
if(fixedq)
    [qixp,qiyp,qizp] = meshgrid(0,0,0);
    Eip = 1;
else
    [qixp,qiyp,qizp] = meshgrid([-Qx,0,Qx]/2,[-Qy,0,Qy]/2,[-Qz,0,Qz]/2);
    Eip     = zeros(27,1);
    Eip(14) = 1;
end
qixyz = [ ...
    [ qixp(:); gi(:,1).*qi; -gi(:,1).*qi ], ...
    [ qiyp(:); gi(:,2).*qi; -gi(:,2).*qi ], ...
    [ qizp(:); gi(:,3).*qi; -gi(:,3).*qi ]    ];
Ei  = [ Eip; atti'; atti' ];
% -------------------------------------------------------------------------
% Actually interpolate the data in the grid.
E   = griddatan( ...
    qixyz, ...
    Ei, ...
    [qox(:),qoy(:),qoz(:)] ); % E is Nl x 1
% Values outside the convex hull computed with griddatan map to NaNs. We
% just place zeros assuming they are out-of-bandwidth:
E(isnan(E)) = 0;
E   = reshape(E,size(qox));
% -------------------------------------------------------------------------
% Compute the 3-D FFT
eap = (Q/numel(E))*fftn(E);
eap = fftshift(eap); % center the spectrum
eap(:,:,2*lattice(3)+1) = eap(:,:,1);
eap(:,2*lattice(1)+1,:) = eap(:,1,:);
eap(2*lattice(2)+1,:,:) = eap(1,:,:);
eap = ( eap + conj(eap(end:-1:1,end:-1:1,end:-1:1)) )/2; % Force antipodal symmetry.
if(const)
    norm = eap(1:end-1,1:end-1,1:end-1);
    norm = sum(abs(norm(:)))/Q;
end
eap = eap( ceil(numel(eap)/2):end );          % keep only unique values; row vector
% -------------------------------------------------------------------------
% Re-normalize if necessary:
if(const)
    eap = abs(eap)./norm;
else
    eap = real(eap);
end
% -------------------------------------------------------------------------
