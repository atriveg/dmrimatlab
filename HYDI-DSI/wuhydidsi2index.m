function index = wuhydidsi2index(eap, lattice, Qx, Qy, Qz, varargin)
% function index = wuhydidsi2index(eap, lattice, Qx, Qy, Qz
%                      'opt1', value1, 'opt2', value2, ..., 'optN', valueN)
%
%   Computes scalar diffusion indices from the Ensemble Average Propagator
%   computed by atti2wuhydidsi, namely RTOP, RTPP, RTAP, and MSD.
%
%   Mandatory inputs (all of them are returned by atti2hydi):
%
%      eap: a MxNxPxK double array containing the values of the EAP at a
%         regular Cartessian lattice. It is computed with atti2wuhydidsi.
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
%      kind: string, one of 'rtop', 'rtpp', 'rtap', 'msd', or 'negativity',
%         describing the actual index to coompute (default: 'rtop').
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
x  = -lattice(1):lattice(1);
y  = -lattice(2):lattice(2);
z  = -lattice(3):lattice(3);
[x,y,z] = meshgrid(x,y,z);
% ---------------------------------------------
ptr    = zeros(size(x));
Nl     = numel(ptr);
ptr(:) = 1:Nl;
ptr(1:floor(Nl/2)) = (Nl:-1:ceil(Nl/2)+1);
ptr    = ptr(1:end-1,1:end-1,1:end-1);
ptr    = ptr(:)' - floor(Nl/2);
% ---------------------------------------------
x = x(1:end-1,1:end-1,1:end-1); x = x(:)';
y = y(1:end-1,1:end-1,1:end-1); y = y(:)';
z = z(1:end-1,1:end-1,1:end-1); z = z(:)';
% ---------------------------------------------
switch(lower(opt.kind))
    case 'rtop'
        indM = eap(:,1);
    case 'rtpp'
        pp   = abs(z)<0.5; % In-plane samples
        indM = sum( eap(:,ptr(pp)), 2 )./(Qx.*Qy);
    case 'rtap'
        pp   = ( (abs(x)<0.5) & (abs(y)<0.5) ); % In-axis samples
        % Count only once the origin and boundary points (avoid duplication
        % due to the periodical conditions):
        indM = sum( eap(:,ptr(pp)), 2 )./(Qz);
    case 'msd'
        x    = (1./Qx)*x;
        y    = (1./Qy)*y;
        z    = (1./Qz)*z;
        rr   = x.*x+y.*y+z.*z;
        indM = sum( eap(:,ptr).*rr, 2 )./(Qx.*Qy.*Qz);
    case 'mass' % Dumb sanity check
        indM = sum( eap(:,ptr), 2 )./(Qx.*Qy.*Qz);
    case 'negativity' % Check negative valeus of the EAP, since it is not enforced to be posiive
        bads = ( eap(:,ptr) < 0 );
        negE = eap(:,ptr).*eap(:,ptr).*double(bads);
        indM = sum(negE,2)./sum(eap(:,ptr).*eap(:,ptr),2);
    otherwise
        error(['Unknown index kind (rtop/rtpp/rtap/msd/negativity): ',opt.kind]);
end
% -------------------------------------------------------------------------
% Reshape back to output the result:
index = zeros(M*N*P,1);
index(mask) = indM;
index = reshape(index,[M,N,P]);
% -------------------------------------------------------------------------
