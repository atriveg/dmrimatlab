function index = hydidsi2index(eap, lattice, Qx, Qy, Qz, varargin)
% function index = hydidsi2index(eap, lattice, Qx, Qy, Qz
%                      'opt1', value1, 'opt2', value2, ..., 'optN', valueN)
%
%   Computes scalar diffusion indices from the Ensemble Average Propagator
%   computed by atti2hydidsi, namely RTOP, RTPP, RTAP, and MSD.
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
%
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      kind: string, one of 'rtop', 'rtpp', 'rtap', or 'msd', describing
%         the actual index to coompute (default: 'rtop').
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).

% Check the mandatory input argments:
if(nargin<5)
    error('At lest the eap volume, the lattice size, and the 3 bandwidths must be supplied');
end
[M,N,P,K] = size(eap);
assert(isequal(size(lattice),[3,1]),'lattice must be 3x1');
assert((prod(2*lattice+1)+1)/2==K,'lattice should match the fourth dimension of eap, type ''help hydi2index'' for details');
assert(isequal(size(Qx),[M,N,P]),'the size of Qx must match the first three dimensiosn of eap');
assert(isequal(size(Qy),[M,N,P]),'the size of Qy must match the first three dimensiosn of eap');
assert(isequal(size(Qz),[M,N,P]),'the size of Qz must match the first three dimensiosn of eap');
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.kind = 'rtop';      optchk.kind = [true,false];    % variable length string
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask = reshape(opt.mask,[M*N*P,1]);
eap  = reshape(eap,[M*N*P,K]);
Qx   = reshape(Qx,[M*N*P,1]);
Qy   = reshape(Qy,[M*N*P,1]);
Qz   = reshape(Qz,[M*N*P,1]);
eap  = eap(mask,:);
Qx   = Qx(mask,:);
Qy   = Qy(mask,:);
Qz   = Qz(mask,:);
% -------------------------------------------------------------------------
% Proceed:
x = -lattice(1):lattice(1);
y = -lattice(2):lattice(2);
z = -lattice(3):lattice(3);
[x,y,z] = meshgrid(x,y,z);
x = x(K:end);
y = y(K:end);
z = z(K:end);
switch(lower(opt.kind))
    case 'rtop'
        indM = eap(:,1);
    case 'rtpp'
        pp   = abs(z)<0.5; % In-plane samples
        pp   = pp(2:end);  % Remove the origin
        indM = 2*sum(eap(:,pp),2)+eap(:,1);
        indM = indM./(Qx.*Qy);
    case 'rtap'
        pp   = ( abs(x)<0.5 & abs(y)<0.5 ); % In-axis samples
        pp   = pp(2:end);                   % Remove the origin
        indM = 2*sum(eap(:,pp),2)+eap(:,1);
        indM = indM./(Qz);
    case 'msd'
        x    = (1./Qx)*x;
        y    = (1./Qy)*y;
        z    = (1./Qz)*z;
        rr   = x.*x+y.*y+z.*z;
        indM = 2*sum(eap.*rr,2);
        indM = indM./(Qx.*Qy.*Qz);
    case 'mass' % Dumb sanity check
        indM = eap(:,1)+2*sum(eap(:,2:end),2);
        indM = indM./(Qx.*Qy.*Qz);
    otherwise
        error(['Unknown index kind (rtop/rtpp/rtap/msd): ',opt.kind]);
end
% -------------------------------------------------------------------------
% Reshape back to output the result:
index = zeros(M*N*P,1);
index(mask) = indM;
index = reshape(index,[M,N,P]);
% -------------------------------------------------------------------------
