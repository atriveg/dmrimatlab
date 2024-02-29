function DTI = shadc2dti(SH,varargin)
% function DTI = shadc2dti(SH,'opt1,val1,...)
%    This function implements the precomputed transition matrices between 
%    the coefficients of a real, symmetric function defined over the unit 
%    sphere in the basis of even Spherical Harmonics up to order L=2 and 
%    the coefficients of this same function in the basis of order-2 
%    tensors. These two sets are known to be basis of the same functional 
%    space. The details on this implementation may be found in:
%
%          Maxime Descoteaux, Elaine Angelino, Shaun Fitzgibbons, and 
%          Rachid Deriche, "Apparent Diffusion Coefficients from High 
%          Angular Resolution Diffusion Images: Estimation and 
%          Applications." Tech. Report 5681, Institut National de Recherche 
%          en Informatique et Automatique, September 2005.
%
%    The input:
%
%       SH: a MxNxPxK, for K=((L+1)(L+2)/2), L=2,4,6... is a double array
%           with the Spherical Harmononics expansion of a given ADC with
%           order at above L=2.
%
%    The output:
%
%       DTI: the corresponding diffusion tensor that best matches th SH 
%            expansion (this is trivial since SH functionscare an 
%            orthonormal basis.). Its size may be (see below):
%               > MxNxPx3x3, it the 'unroll' option is true (default);
%               > MxNxPx6, if the 'unroll' option is false.
%
%    Optional arguments may be passed as name/value pairs in the regular
%    matlab style:
%
%       unroll: wether (true) or not (false) output the result as a 3x3
%         matrix at each voxel instead of a 6x1 vector (in the latter case
%         duplicates of the entries of the diffusion tensor are removed so
%         that only [D11,D12,D13,D22,D23,D33] are returned) (default:
%         true).
%
%       chunksz: the computation reduces to the product of the SH coeffs
%         by a matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the 
%         corresponding matrix, hence taking advantage of matlab's
%         capabilities (default: 1000).
%
%       mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% Check the mandatory input argments:
if(nargin<1)
    error('At lest the coefficients volume must be supplied');
end
[M,N,P,K] = size(SH);
NV = M*N*P;                % Total number of voxels to be processed
L = (sqrt(8*(K-1)+9)-3)/2; % SH order
if( abs(L-round(L))>1.0e-9 || L<2 )
    error('Weird size of the SH volume. Its fourth dimension should have size 6, 15, 28, 45, ..., (L+1)(L+2)/2, with L=2,4,6,...');
end
% Keep SH of orders 0 and 2 and remove the rest. Since SH are orthonormal,
% this is actually the projection of the signal onto the space of order-2
% symmetric tensors:
SH = SH(:,:,:,1:6); % MxNxPx6

% Parse the optional input arguments:
opt.unroll = true;      optchk.unroll = [true,true];  % always 1x1 boolean
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Transition matrix analitically computed with Maple:
B = [ ...
    2/3,         0,           2/3,        0,          0,           2/3;
    0,           0,          -2/sqrt(15), 0,          0,           2/sqrt(15);
    0,           0,           0,          4/sqrt(15), 0,           0;
    4/sqrt(45),  0,          -2/sqrt(45), 0,          0,          -2/sqrt(45);
    0,          -4/sqrt(15),  0,          0,          0,           0;
    0,           0,           0,          0,          4/sqrt(15),  0
    ];
B = B*sqrt(pi); % 6x6
B = B\eye(6);   % B^(-1), 6x6
B = B'; % For convenience, see loop below

% Now, process the data chunk-by chunk where the mask is true:
SH   = reshape(SH,[NV,6]);   % NVx6
mask = opt.mask(:);          % NVx1
% Mask...
SH      = SH(mask,:);  % PVx6
PV      = size(SH,1);
DTImask = zeros(PV,6); % PVx6
for ck=1:ceil(PV/opt.chunksz)
    idi = (ck-1)*opt.chunksz+1;
    idf = min(ck*opt.chunksz,PV);
    DTImask(idi:idf,:) = SH(idi:idf,:)*B; % (chunksz x 6) * (6 x 6) -> (chunksz x 6)
end
% Cast the result to the proper size:
DTIshrink = zeros(NV,6); % NVx6
DTIshrink(mask,:) = DTImask;
DTIshrink = reshape(DTIshrink,[M,N,P,6]);
% Cast into the full tensor:
if(opt.unroll)
    DTI = zeros(M,N,P,3,3);
    DTI(:,:,:,1,1) = DTIshrink(:,:,:,6);
    DTI(:,:,:,1,2) = DTIshrink(:,:,:,5);
    DTI(:,:,:,1,3) = DTIshrink(:,:,:,4);
    DTI(:,:,:,2,1) = DTIshrink(:,:,:,5);
    DTI(:,:,:,2,2) = DTIshrink(:,:,:,3);
    DTI(:,:,:,2,3) = DTIshrink(:,:,:,2);
    DTI(:,:,:,3,1) = DTIshrink(:,:,:,4);
    DTI(:,:,:,3,2) = DTIshrink(:,:,:,2);
    DTI(:,:,:,3,3) = DTIshrink(:,:,:,1);
else
    DTI = DTIshrink(:,:,:,end:-1:1);
end
