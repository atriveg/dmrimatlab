function [eap,dti,lattice,Qx,Qy,Qz,res,lapl,lopt] = atti2hydidsi( atti, gi, bi, varargin )
% function [eap,dti,lattice,Qx,Qy,Qz,res,lapl,lopt] = atti2hydidsi( atti, gi, bi, 
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Computes the Ensemble Average Propagator (eap) at a regular Cartesian
%   lattice from a set of scattered q-space measurements (atti) using
%   Hybrid Diffusion Imaging without q-space re-gridding. To that end, the
%   R-space is first aligned with the eigenvectors of a tensor model (the z
%   axis aligns with the maximum diffusion direction) and an appropriate
%   lattice is designed so that the whole extent of the EAP is covered
%   according to the computed eigenvalues. From this estimate, the q-space
%   is cropped to avoid aliasing, meaning that some q-space samples might
%   be discarded for some voxels. Accordingly, the tensor model estimated
%   is necessary to interpret the computed eap (and it is returned by this
%   function). The eap is estimated by solving a quadratic programming
%   problem that involves the encoding matrix that relates q-space
%   measurments to R-space Fourier coefficients (i.e. EAP samples).
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
%         dti (and scaled according to its eigenvalues).
%      lattice: 3x1, a double array with the radii of the lattice at each
%         voxel, i.e. [Nx;Ny;Nz]. This is useful in case this parameter is
%         automatically determined instead of fed as an input argument.
%      Qx, Qy, Qz: MxNxP each, double arrays with the estimated bandwidths
%         of the q-space signal along each dimension. They are necessary to
%         interpret the computed eap.
%      res: MxNxP, the residual of the QP problem at each voxel.
%      lapl: MxNxP, the energy of the Laplacian at each voxel.
%      lopt: MxNxP, the optimal value of lambda computed at each voxel
%         (makes sense only if Generalized Cross Validation is used to
%         optimally set the value of lambda at each voxel, see below).
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   General optional parameters:
%
%      lambda: this is a regularization parameter that penalizes highly
%         irregular solutions of the EAP. In precise terms, the quadratic
%         deviation from the linear model is added lambda times the energy
%         of the Laplacian of the EAP to compute the final cost function.
%         NOTE: you can pass Inf, NaN or any negative value of lambda to
%         the fucntion. In this case, its value will be adaptively computed
%         at each voxel using Generalized Cross Validation (GCV), BUT: this
%         will blow up the overall computation time, since GCV implies
%         repeated matrix inversions (default 1.0e-3).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
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
%      const: wether (true) or not (false) constrain the estimated EAP to
%         be positive and have unit mass. With these constraints, the
%         estimation has a rigorous physical meaning but a quadratic
%         problem must be solved at each voxel. Without it, we just need to
%         solve an unconstrained least squares problem, which will be
%         pretty faster (default: true, i.e. use constraints).
%
%   Options when const is true:
%
%         These options are directly passed to dmriQuadprog (check the 
%         documentation therein):
%
%         miters (default: 100)
%         otol (default: 1.0e-6)
%         stol (default: 1.0e-6)
%         ctol (default: 1.0e-8)
%
%   Other general options:
%
%      usemex: whether (true) or not (false) use the C/MEX implementation
%         of the method. This implementation can be notably faster (or not)
%         depending on the CPU used. It makes extensive use of BLAS
%         routines, so the performance is computer dependent. In POSIX
%         systems the implementation is multi-threaded, so that you don't
%         need a working license for the parallel computing toolbox to run
%         the code in parallel (default: false). NOTE: if 'usemex' is set
%         true, the 'const' option is ignored and the constrained problem
%         is always solved.
%      maxthreads: if 'usemex' is set true and this is running in a POSIX
%         system, the maximum number of threads used (default:
%         automatically determined).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).
%      verbose: whether (true) or not (false) show a progress bar to check
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
opt.lambda = 1.0e-3;    optchk.lambda = [true,true];     % always 1x1 double
opt.tl = 1.0e-7;        optchk.tl = [true,true];         % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];         % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.tau = 35.0e-3;      optchk.tau = [true,true];        % always 1x1 double
opt.Rth = 0.01;         optchk.Rth = [true,true];        % always 1x1 double
opt.lattice = [];       optchk.lattice = [true,false];   % variable size
opt.const = true;       optchk.const = [true,true];      % always 1x1 boolean
opt.usemex = false;     optchk.usemex = [true,true];     % always 1x boolean
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true];
% -------------------------------------------------------------------------
opt.miters = 200;       optchk.miters = [true,true];
opt.otol = 1.0e-6;      optchk.otol = [true,true];
opt.stol = 1.0e-6;      optchk.stol = [true,true];
opt.ctol = 1.0e-8;      optchk.ctol = [true,true];
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
opt.verbose = false;    optchk.verbose = [true,true];    % always 1x1 boolean
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
use_parallel = use_parallel_test;
% -------------------------------------------------------------------------
% Estimate the tensor model:
dti = atti(:,:,:,pdti);  % M x N x P x Ndti
dti = atti2dti( dti, gi(pdti,:), bi(pdti), 'mask', opt.mask, ...
    'wls', true, 'nonlinear', true, 'wlsit', 3, 'wsc', 0.01, ...
    'tol', 1.0e-12, 'rcondth', 1.0e-6, 'fixmode', 'a', ...
    'maxthreads', opt.maxthreads, 'unroll', false );
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask  = reshape(opt.mask,[M*N*P,1]);
atti  = reshape(atti,[M*N*P,G]);
atti  = atti(mask,:);
dti2  = reshape(dti,[M*N*P,6]);
dti2  = dti2(mask,:);
Q     = size(atti,1);
eapM  = zeros(Q,K);
QxM   = zeros(Q,1);
QyM   = zeros(Q,1);
QzM   = zeros(Q,1);
resm  = zeros(Q,1);
laplm = zeros(Q,1);
loptm = zeros(Q,1);
% -------------------------------------------------------------------------
% No other way than looping voxel-to-voxel. For the parallel
% implementation, extract parameters from the opt structure to reduce
% memory exchange:
verbose = opt.verbose;
lambda  = opt.lambda;
ADC0    = opt.ADC0;
tau     = opt.tau;
lRth    = -4*tau*log(opt.Rth);
const   = opt.const;
qi      = sqrt(bi/(4*pi*pi*tau));
optim.miters = opt.miters;
optim.otol   = opt.otol;
optim.stol   = opt.stol;
optim.ctol   = opt.ctol;

% Compute the DFT matrix for the Laplacian penalty:
[DFT,u,v,w] = dft_matrix(lattice);
% -----
if(opt.usemex)
    [eapM,QxM,QyM,QzM,resm,laplm,loptm] = atti2hydidsi_( double(atti'), double(dti2'), gi, qi, ...
        lattice, lambda, DFT, [u(:),v(:),w(:)], ADC0, lRth, optim, opt.maxthreads );
    eapM = eapM';
else
    NVs   = 50;
    Qchnk = round(Q/NVs);
    if(verbose)
        str   = repmat('-',[1,NVs]);
        fprintf(1,'Completed: [%s]',str);
    end
    if(use_parallel) % Parallel implementation
        parfor q=1:Q
            [eap,Qx,Qy,Qz,resn,lapln,lopt] = one_voxel_eap( atti(q,:), dti2(q,:), gi, qi, ...
                lattice, lambda, DFT, u, v, w, ADC0, lRth, const, optim );
            eapM(q,:) = eap;
            QxM(q)    = Qx;
            QyM(q)    = Qy;
            QzM(q)    = Qz;
            resm(q)   = resn;
            laplm(q)  = lapln;
            loptm(q)  = lopt;
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
            [eap,Qx,Qy,Qz,resn,lapln,lopt] = one_voxel_eap( atti(q,:), dti2(q,:), gi, qi, ...
                lattice, lambda, DFT, u, v, w, ADC0, lRth, const, optim, q );
            eapM(q,:) = eap;
            QxM(q)    = Qx;
            QyM(q)    = Qy;
            QzM(q)    = Qz;
            resm(q)   = resn;
            laplm(q)  = lapln;
            loptm(q)  = lopt;
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
end
% -------------------------------------------------------------------------
% Reshape back to output the result:
eap  = zeros(M*N*P,K);
Qx   = zeros(M*N*P,1);
Qy   = zeros(M*N*P,1);
Qz   = zeros(M*N*P,1);
res  = zeros(M*N*P,1);
lapl = zeros(M*N*P,1);
lopt = zeros(M*N*P,1);
eap(mask,:) = eapM;
Qx(mask)    = QxM;
Qy(mask)    = QyM;
Qz(mask)    = QzM;
res(mask)   = resm;
lapl(mask)  = laplm;
lopt(mask)  = loptm;
eap  = reshape(eap,[M,N,P,K]);
Qx   = reshape(Qx,[M,N,P]);
Qy   = reshape(Qy,[M,N,P]);
Qz   = reshape(Qz,[M,N,P]);
res  = reshape(res,[M,N,P]);
lapl = reshape(lapl,[M,N,P]);
lopt = reshape(lopt,[M,N,P]);
% -------------------------------------------------------------------------

function [eap,Qx,Qy,Qz,resn,lapln,lopt] = one_voxel_eap( atti, dti, gi, qi, lattice, ...
    lambda, DFT, u, v, w, ADC0, lRth, const, optim, q )
% The function that does it all:
% atti:    1 x G
% dti:     1 x 6
% gi:      G x 3
% bi:      G x 1
% lattice: 3 x 1 -> K = (prod(2*lattice+1)+1)/2
% lambda:  1 x 1
% DFT:     K' x K
% u,v,w:   K' x 1
% ADC0:    1 x 1
% eap:     1 x K
% -------------------------------------------------------------------------
% Check the tensor model
D     = dti([1,2,3;2,4,5;3,5,6]);
if(any(isnan(D)))
    U = eye(3);
    L = ADC0*ones(3,1);
else
    [U,L] = eig(D);
end
L     = abs(diag(L));
L(L>ADC0)    = ADC0;
L(L<ADC0/60) = ADC0/60;
[L,id]       = sort(L);
y     = U(:,id(2));
z     = U(:,id(3));
x     = cross(y,z);
% -------------------------------------------------------------------------
% Determine the R- and q-domains support
gi = gi*[x,y,z];
Qx = lattice(1)/sqrt(L(1)*lRth);
Qy = lattice(2)/sqrt(L(2)*lRth);
Qz = lattice(3)/sqrt(L(3)*lRth);
Q  = Qx*Qy*Qz;
% -------------------------------------------------------------------------
% Determine which q-samples will be actually used:
qx = gi(:,1).*qi;
qy = gi(:,2).*qi;
qz = gi(:,3).*qi;
bw = ( abs(qx)<Qx/2 & abs(qy)<Qy/2 & abs(qz)<Qz/2 );
qx = qx(bw);
qy = qy(bw);
qz = qz(bw);
% -------------------------------------------------------------------------
% Sample the R-domain:
rx = (-lattice(1):lattice(1))/Qx;
ry = (-lattice(2):lattice(2))/Qy;
rz = (-lattice(3):lattice(3))/Qz;
[rx,ry,rz] = meshgrid(rx,ry,rz);
Nl = ceil(numel(rx)/2);
rx = rx(Nl:end);
ry = ry(Nl:end);
rz = rz(Nl:end);
% -------------------------------------------------------------------------
% Create the encoding matrix:
Fe = 2*pi*(qx*rx+qy*ry+qz*rz); % Ni x K
Fe = cos(Fe)/Q;                % Ni x K
Fe(:,2:end) = 2*Fe(:,2:end);
% -------------------------------------------------------------------------
% Create the Laplacian penalty:
u   = Qx*Qx*u/Q^(2/3);
v   = Qy*Qy*v/Q^(2/3);
w   = Qz*Qz*w/Q^(2/3);
DFT = DFT/Q;
DFT = diag(u+v+w)*DFT;
% -------------------------------------------------------------------------
%
%
%
% %%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTAL CODE (GCV)
if( isinf(lambda) || isnan(lambda) || (lambda<0) )
    lambda = dmri_gcv( Fe, DFT, double(atti(1,bw)'), ...
        1, ...        % initial (maximum) value probed
        0.000001, ... % final (minimum) value probed
        100 );        % maximum number of steps
    if(isnan(lambda))
        lambda = 0.00001;
    end
end
lopt = lambda;
% %%%%%%%%%%%%%%%%%%%%%%%
%
%
%
% -------------------------------------------------------------------------
% Solve the (regularized, constrained) quadratic programming problem
Fe0 = 2*ones(size(rx));
Fe0(1) = 1;
Fe0 = Fe0/Q;
H   = (Fe')*Fe+lambda*((DFT')*DFT);
b   = -(Fe')*(double(atti(1,bw)'));
% Initialize the estimate as an unconstrained problem in all cases:
eap = -H\b;
% Ensure the fulfillment of the constraints:
eap(eap<0) = 0;
eap = eap/(Fe0*eap);
if(const) % Constrained problem -> QP
    % Instead of using quadprog:
    %
    %    eap = quadprog(H,b,[],[],Fe0,1,zeros(size(eap)),[],eap,optim);
    %
    % We use a home-made function:
    [eap,~,~] = dmriQuadprog(H,b,Fe0',eap, ...
        optim.miters, ...
        optim.otol, ...
        optim.stol, ...
        optim.ctol );
end
resn  = ( double(atti(1,bw)') - Fe*eap );
lapln = DFT*eap;
resn  = (resn')*resn;
lapln = (lapln')*lapln;
eap   = eap';
% -------------------------------------------------------------------------

function [DFT,u,v,w] = dft_matrix(lattice)
% Computes a matrix that represents an ad-hoc Discrete Fourier Transform
% for this model, so that the Laplacian weighting will result from the
% product of q^2 with this one.
%
%   lattice: 3 x 1, [Nx;Ny;Nz]
% -------------------------------------------------------------------------
% Extended x-y-z domain; we need an extra border at +-(Nx+1), +-(Ny+1),
% +-(Nz+1) padded with zeros so that P(R) can be periodic without
% discontinuities. However, the samples at -(Nx+1), -(Ny+1), -(Nz+1) are
% not included since they must equal those at Nx+1, Ny+1, Nz+1
k = -lattice(1):lattice(1)+1;
l = -lattice(2):lattice(2)+1;
m = -lattice(3):lattice(3)+1;
[k,l,m] = meshgrid(k,l,m);
k = k(:);
l = l(:);
m = m(:);
% -------------------------------------------------------------------------
% Regular (complex) DFT matrix, so that:
%     P = Q/(8(Nx+1)(Ny+1)(Nz+1)) * DFT*E
DFT = ...
    (k*(k'))/(lattice(1)+1) + ...
    (l*(l'))/(lattice(2)+1) + ...
    (m*(m'))/(lattice(3)+1);    % 8(Nx+1)(Ny+1)(Nz+1) x 8(Nx+1)(Ny+1)(Nz+1) 
DFT = exp(1i*pi*DFT);           % 8(Nx+1)(Ny+1)(Nz+1) x 8(Nx+1)(Ny+1)(Nz+1)
% -------------------------------------------------------------------------
% Since DFT is Hermitic, its inversion is trivial:
%     iDFT = (DFT')/(8(Nx+1)(Ny+1)(Nz+1));
% however, the term 8(Nx+1)(Ny+1)(Nz+1) will cancel the other one:
%     E = 8(Nx+1)(Ny+1)(Nz+1)/Q * iDFT*P
%       = 8(Nx+1)(Ny+1)(Nz+1)/Q * 1/(8(Nx+1)(Ny+1)(Nz+1)) * (DFT')*P
%       = 1/Q (DFT')*P
% so that, we just take the complex transpose:
DFT = DFT';
% -------------------------------------------------------------------------
% Since both E and P are antipodal symmetric, we can remove those columns
% (and rows) corresponding to the hemisphere z<0 (and qz<0):
pp1 = (   (m<-1/2)   ...
    |   ( abs(m)<1/2 & k<-1/2 )   ...
    |   ( abs(m)<1/2 & abs(k)<1/2 & l<-1/2)   );
% -------------------------------------------------------------------------
% Since the border of P is zero-padded, we can remove the corresponding
% columns:
pp2 = ( (k>lattice(1)+1/2) | (l>lattice(2)+1/2) | (m>lattice(3)+1/2) );
% -------------------------------------------------------------------------
% Finally, the origin of the q-space can be discarded since it will be
% 0-weighted and will not contribute to the Laplacian:
pp3 = ( (abs(k)<1/2) & (abs(l)<1/2) & (abs(m)<1/2) );
% -------------------------------------------------------------------------
% Arrange the final matrix (factor 2 stands for the antipodal symmetry of
% P):
DFT = real(   DFT( ~pp1 & ~pp3, ~pp1 & ~ pp2 )   );
DFT(:,2:end) = 2*DFT(:,2:end);
u   = k( ~pp1 & ~pp3 ); % Return q-locations as well
v   = l( ~pp1 & ~pp3 ); % Return q-locations as well
w   = m( ~pp1 & ~pp3 ); % Return q-locations as well
u   = -pi*pi*u.*u./(lattice(1)+1)/(lattice(1)+1);
v   = -pi*pi*v.*v./(lattice(2)+1)/(lattice(2)+1);
w   = -pi*pi*w.*w./(lattice(3)+1)/(lattice(3)+1);
% -------------------------------------------------------------------------
