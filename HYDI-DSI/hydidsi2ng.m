function ng = hydidsi2ng(eap, lattice, Qx, Qy, Qz, l1, l2, l3, tau, varargin)
% function ng = hydidsi2index(eap, lattice, Qx, Qy, Qz, l1, l2, l3, tau,
%                      'opt1', value1, 'opt2', value2, ..., 'optN', valueN)
%
%   Computes the Non-gaussianity of the EAP computed with HYDI-DSI, i.e.
%   the fraction of the energy of the EAP that cannot be explained by the
%   Gaussian approximation
%
%   Mandatory inputs (all of them are returned by atti2hydi):
%
%      eap: a MxNxPxK double array containing the values of the EAP at a
%         regular Cartessian lattice. It is computed with atti2hydidsi.
%      lattice: 3x1, a double array with the radii of the lattice at each
%         voxel, i.e. [Nx;Ny;Nz]. MUST fulfill:
%            ((2*Nx+1)*(2*Ny+1)*(2*Nz+1)+1)/2 = K.
%      Qx, Qy, Qz: MxNxP each, double arrays with the estimated bandwidths
%         of the q-space signal along each dimension.
%      l1, l2, l3: MxNxP each, double arrays with the eigenvalues of the
%         tensor model (acending order). You can retrieve these inputs,
%         from the "dti" output of atti2hydidsi, calling dti2spectrum.
%         BEWARE: dti2spectrum will return the eigenvalues in descending
%         order, but this function asks for them in ascending order.
%      tau: 1x1, the effective diffusion time of the acquisition, necessary
%         to determine the actual support of the q-space. MUST match the
%         value used in the call to atti2hydidsi.
%
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).

% Check the mandatory input argments:
if(nargin<9)
    error('At lest the eap volume, the lattice size, the 3 bandwidths, the 3 eigenvalues, and tau must be supplied');
end
[M,N,P,K] = size(eap);
assert(isequal(size(lattice),[3,1]),'lattice must be 3x1');
assert((prod(2*lattice+1)+1)/2==K,'lattice should match the fourth dimension of eap, type ''help hydi2index'' for details');
assert(isequal(size(Qx),[M,N,P]),'the size of Qx must match the first three dimensiosn of eap');
assert(isequal(size(Qy),[M,N,P]),'the size of Qy must match the first three dimensiosn of eap');
assert(isequal(size(Qz),[M,N,P]),'the size of Qz must match the first three dimensiosn of eap');
assert(isequal(size(l1),[M,N,P]),'the size of l1 must match the first three dimensiosn of eap');
assert(isequal(size(l2),[M,N,P]),'the size of l2 must match the first three dimensiosn of eap');
assert(isequal(size(l3),[M,N,P]),'the size of l3 must match the first three dimensiosn of eap');
assert(isscalar(tau),'tau must be a scalar');
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask = reshape(opt.mask,[M*N*P,1]);
eap  = reshape(eap,[M*N*P,K]);
Qx   = reshape(Qx,[M*N*P,1]);
Qy   = reshape(Qy,[M*N*P,1]);
Qz   = reshape(Qz,[M*N*P,1]);
l1   = reshape(l1,[M*N*P,1]);
l2   = reshape(l2,[M*N*P,1]);
l3   = reshape(l3,[M*N*P,1]);
eap  = eap(mask,:); % S x K
Qx   = Qx(mask,:);  % S x 1
Qy   = Qy(mask,:);  % S x 1
Qz   = Qz(mask,:);  % S x 1
l1   = l1(mask,:);  % S x 1
l2   = l2(mask,:);  % S x 1
l3   = l3(mask,:);  % S x 1
% -------------------------------------------------------------------------
% Prepare the lattice where computations are done:
x = -lattice(1):lattice(1);
y = -lattice(2):lattice(2);
z = -lattice(3):lattice(3);
[x,y,z] = meshgrid(x,y,z);
x = x(K:end);
y = y(K:end);
z = z(K:end);
x = x(:)'; % 1 x K
y = y(:)'; % 1 x K
z = z(:)'; % 1 x K
% -------------------------------------------------------------------------
% Create the tensor approximation for this lattice:
if(is_broadcast_available)
    x = x./Qx; % S x K
    y = y./Qy; % S x K
    z = z./Qz; % S x K
else
    x = bsxfun( @(a,b)(a./b), x, Qx ); % S x K
    y = bsxfun( @(a,b)(a./b), y, Qy ); % S x K
    z = bsxfun( @(a,b)(a./b), z, Qz ); % S x K
end
% Compute the argument of the exponential function:
x = x.*x; % S x K
y = y.*y; % S x K
z = z.*z; % S x K
if(is_broadcast_available)
    eargx = x./l1; % S x K
    eargy = y./l2; % S x K
    eargz = z./l3; % S x K
else
    eargx = bsxfun( @(a,b)(a./b), x, l1 ); % S x K
    eargy = bsxfun( @(a,b)(a./b), y, l2 ); % S x K
    eargz = bsxfun( @(a,b)(a./b), z, l3 ); % S x K
end
earg = -(eargx+eargy+eargz)/4/tau;
% Compute the Gaussian EAP:
fptau = 4*pi*tau;
den   = sqrt(abs(l1.*l2.*l3)*fptau*fptau*fptau); % S x 1
geap  = exp(earg);
if(is_broadcast_available)
    geap = geap./den;
else
    geap = bsxfun( @(a,b)(a./b), geap, den );
end
% -------------------------------------------------------------------------
% Numerically approximate the NG:
ngM = eap-geap; % S x M
ngM = sqrt(sum(ngM.*ngM,2)./sum(eap.*eap,2));
% -------------------------------------------------------------------------
% Reshape back to output the result:
ng = zeros(M*N*P,1);
ng(mask) = ngM;
ng = reshape(ng,[M,N,P]);
% -------------------------------------------------------------------------
