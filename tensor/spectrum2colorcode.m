function varargout = spectrum2colorcode(u1,l1,l2,l3,varargin)
% [R,G,B] = spectrum2colorcode(u1,l1,l2,l3,varargin)
%     RGB = spectrum2colorcode(u1,l1,l2,l3,varargin)
%
%   Computes the standard color-code representation for diffusion tensor
%   data based on their spectra as computed by the dti2spectrum function.
%
%   All mandatory inputs u1, l1, l2, and l3 are returnned by dti2spectrum.
%
%   INPUTS:
%
%      u1: a MxNxPx3 double array with the unit-norm eigenvector associated
%         to the largest eigenvalue, l1.
%      l1: a MxNxP volume with the first (largest) eigenvalue as provided
%         by dti2spectrum.
%      l2: a MxNxP volume with the second eigenvalue as provided by
%         dti2spectrum.
%      l3: a MxNxP volume with the third (samlles) eigenvalue as provided
%         by dti2spectrum.
%
%   OUPUTS (when called with 1 < output arguments <= 3 ):
%
%      R: a MxNxP volume with the red (x-alignment) channel of the
%         color-code map weighted by the FA, in the range [0,1].
%      G: a MxNxP volume with the red (y-alignment) channel of the
%         color-code map weighted by the FA, in the range [0,1].
%      B: a MxNxP volume with the red (z-alignment) channel of the
%         color-code map weighted by the FA, in the range [0,1].
%
%   OUPUTS (when called with 1 output argument):
%
%      RGB: a MxNxPx3 volume with the stacked RGB (xyz) channels weigthed
%         by the FA, in the range [0,1].
%
%   For example, you can plot the 17-th slice of a "tensor" volume as:
%
%      >> [u1,~,~,l1,l2,l3] = dti2spectrum( tensor, 'mask', mask );
%      >> RGB = spectrum2colorcode( u1, l1, l2, l3, 'mask', mask);
%      >> slice = 17;
%      >> IMG = squeeze(RGB(:,:,slice,:));
%      >> imshow(IMG);
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% -------------------------------------------------------------------------
assert(ndims(u1)==4,'Argument u1 must be MxNxPx3');
[M,N,P,Q] = size(u1);
assert(Q==3,'Argument u1 must be MxNxPx3');
assert(isequal(size(l1),[M,N,P]),'l1 must be 3-D, and its 3 dimensions must match those of u1');
assert(isequal(size(l2),[M,N,P]),'l2 must be 3-D, and its 3 dimensions must match those of u1');
assert(isequal(size(l3),[M,N,P]),'l3 must be 3-D, and its 3 dimensions must match those of u1');
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];  % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Compute the FA:
FA = spectrum2scalar(l1,l2,l3,'scalar','fa','mask',opt.mask);
% -------------------------------------------------------------------------
% Unroll to work comfortably:
u1 = reshape(u1,[M*N*P,3]);
u1 = u1(opt.mask,:);
FA = FA(opt.mask);
% -------------------------------------------------------------------------
% Compute channels:
Rm = abs(u1(:,1)).*FA;
Gm = abs(u1(:,2)).*FA;
Bm = abs(u1(:,3)).*FA;
% -------------------------------------------------------------------------
% Reshape:
R = zeros(M,N,P);
R(opt.mask) = Rm;
G = zeros(M,N,P);
G(opt.mask) = Gm;
B = zeros(M,N,P);
B(opt.mask) = Bm;
% -------------------------------------------------------------------------
% Assign outputs:
if(nargout==1)
    RGB = zeros(M,N,P,3);
    RGB(:,:,:,1) = R;
    RGB(:,:,:,2) = G;
    RGB(:,:,:,3) = B;
    varargout{1} = RGB;
elseif(nargout==2)
    varargout{1} = R;
    varargout{2} = G;
elseif(nargout==3)
    varargout{1} = R;
    varargout{2} = G;
    varargout{3} = B;
else
    error('Too many output arguments');
end
% -------------------------------------------------------------------------
end
