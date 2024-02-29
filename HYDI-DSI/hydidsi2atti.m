function atti = hydidsi2atti( eap, dti, lattice, Qx, Qy, Qz, gi, bi, varargin )
% function atti = hydidsi2atti( eap, dti, lattice, Qx, Qy, Qz, gi, bi, 
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Predicts the value of the attenuation signal, according to the HYDI-DSI
%   representation retrieved with the function atti2hydidsi, for a set of
%   desired q-space locations specified by directions gi and b-values bi
%
%   Mandatory inputs:
%
%      eap: a MxNxPxK double array containing a regular sampling of the
%         positive, unit mass Ensemble Average Propagator estimated with
%         atti2hydidsi. To interpret this sampling, the following five
%         input arguments (as returned by atti2hydidsi) are mandatory:
%      dti: MxNxPx6, a double array with the estimated tensor model  at
%         each voxel. The EAP lattice is aligned at each voxel to the 
%         eigenvectors of dti (and scaled according to its eigenvalues).
%      lattice: 3x1, a double array with the radii of the lattice at each
%         voxel, i.e. [Nx;Ny;Nz].
%      Qx, Qy, Qz: MxNxP each, double arrays with the estimated bandwidths
%         of the q-space signal along each dimension.
%
%     Besides, you must provide the q-space points where the signal
%     representation will be evaluated. THESE POINTS ARE NOT NECESSARILY
%     THE SAME (BUT MIGHT BE) AS THOSE USED TO CALL atti2hydidsi:
%
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that a multi-shell acquisition is arranged.
%
%   The outputs of the function area as follows:
%
%      atti: a MxNxPxG array with the attenaution signal predicted at each
%         voxel
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      tau: the effective diffusion time of the acquisition, necessary to
%         determine the actual support of the q-space. THE SAME USED IN
%         atti2hydidsi!!! (default: 35e-3).
%      ADC0: estimated diffusivity of free water at body temperature. It is
%         used to determine the lower and upper bounds of the eigenvalues
%         of the dti and perform sanity checks (default: 3.0e-3).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).

% Check the mandatory input argments:
if(nargin<8)
    error('A minimum of 8 parameters are required, check the help on this functon');
end
[M,N,P,K] = size(eap);
assert(isequal(size(lattice),[3,1]),'lattice should have size 3x1');
assert(K==(prod(2*lattice+1)+1)/2,'The fourth dimension of EAP does not match the lattice size');
if(~ismatrix(gi)||~ismatrix(bi))
    error('gi and bi must be 2-d matlab matrixes');
end
if(size(gi,1)~=size(bi,1))
    error('The number of rows in gi must match the number of entries of bi');
end
G = size(gi,1);
if(size(gi,2)~=3)
    error('The gradients table gi must have size Gx3');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end
assert(isequal(size(dti),[M,N,P,6]),'dti should have size MxNxPx6');
assert(isequal(size(Qx),[M,N,P]),'The size of Qx does not match the first three dimensions of eap');
assert(isequal(size(Qy),[M,N,P]),'The size of Qy does not match the first three dimensions of eap');
assert(isequal(size(Qz),[M,N,P]),'The size of Qz does not match the first three dimensions of eap');
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.tau  = 35.0e-3;     optchk.tau = [true,true];        % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask  = reshape(opt.mask,[M*N*P,1]);
eap   = reshape(eap,[M*N*P,K]);
eap   = eap(mask,:);
dti   = reshape(dti,[M*N*P,6]);
dti   = dti(mask,:);
Qx    = reshape(Qx,[M*N*P,1]);
Qx    = Qx(mask,:);
Qy    = reshape(Qy,[M*N*P,1]);
Qy    = Qy(mask,:);
Qz    = reshape(Qz,[M*N*P,1]);
Qz    = Qz(mask,:);
Q     = size(eap,1);
attim = zeros(Q,G);
% -------------------------------------------------------------------------
% Proceed:
tau     = opt.tau;
ADC0    = opt.ADC0;
qi      = sqrt(bi/(4*pi*pi*tau));
if(use_parallel_test) % Parallel implementation
    parfor q=1:Q
        PHI = one_voxel_phi( dti(q,:), lattice, Qx(q), Qy(q), Qz(q), gi, qi, ADC0 );
        attim(q,:) = eap(q,:)*(PHI');
    end
else             % Non-parallel implementation
    for q=1:Q
        PHI = one_voxel_phi( dti(q,:), lattice, Qx(q), Qy(q), Qz(q), gi, qi, ADC0 );
        attim(q,:) = eap(q,:)*(PHI');
    end
end
%--------------------------------------------------------------------
% Reshape back to output the result:
atti = zeros(M*N*P,G);
atti(mask,:) = attim;
atti = reshape(atti,[M,N,P,G]);
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PHI = one_voxel_phi( dti, lattice, Qx, Qy, Qz, gi, qi, ADC0 )
% The function that does it all:
% dti:        1 x 6
% lattice:    3 x 1 -> K = (prod(2*lattice+1)+1)/2
% Qx, Qy, Qz: 1 x 1
% gi:         G x 3
% qi:         G x 1
% ADC0:       1 x 1
% PHI:        G x K
% -------------------------------------------------------------------------
% Check the tensor model as in atti2hydidsi
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
[~,id]       = sort(L);
y     = U(:,id(2));
z     = U(:,id(3));
x     = cross(y,z);
% -------------------------------------------------------------------------
% Determine the R- and q-domains support
gi = gi*[x,y,z];
Q  = Qx*Qy*Qz;
% -------------------------------------------------------------------------
% Determine which q-samples will be actually used:
qx = gi(:,1).*qi;
qy = gi(:,2).*qi;
qz = gi(:,3).*qi;
bw = ( abs(qx)<Qx/2 & abs(qy)<Qy/2 & abs(qz)<Qz/2 ); % Other rows will be just cropped to zero
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
PHI = 2*pi*(qx*rx+qy*ry+qz*rz); % Ni x K
PHI = cos(PHI)/Q;                % Ni x K
PHI(:,2:end) = 2*PHI(:,2:end);
PHI(~bw,:) = 0; % Out of bandwidth values
