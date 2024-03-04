function vol = spectrum2scalar(l1,l2,l3,varargin)
% vol = spectrum2scalar( l1, l2, l3, 'option1', value1, ... )
%
%   Computes diffusion tensor-related scalar measures from the spectrum
%   of the diffusion tensor, provided by the three MxNxP volumes of
%   eigenvalues l1, l2, and l3.
%
%   These eigenvalues can be computed from the diffusion tensor using the
%   dti2spectrum function.
%
%   INPUTS:
%
%      l1: a MxNxP volume with the first (largest) eigenvalue as provided
%         by dti2spectrum.
%      l2: a MxNxP volume with the second eigenvalue as provided by
%         dti2spectrum.
%      l3: a MxNxP volume with the third (samlles) eigenvalue as provided
%         by dti2spectrum.
%
%   OUPUTS:
%
%      vol: a MxNxP scalar map with the scalar measure computed by the
%      function.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      scalar: A string, one of:
%
%         + 'fa' [DEFAULT]: the fractional anisotropy is computed.
%         + 'md': the mean diffusivity is computed.
%         + 'ra': the relative anisotropy is computed.
%         + 'cl': Westin's linear coefficient is computed.
%         + 'cp': Westin's planar coefficient is computed.
%         + 'cs': Westin's spherical coefficient is computed.
%
%        Westin's coeffcient can also be computed with a RMS normalizing
%        factor as:
%
%         + 'clsq': Westin's linear coefficient is computed.
%         + 'cpsq': Westin's planar coefficient is computed.
%         + 'cssq': Westin's spherical coefficient is computed.
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% -------------------------------------------------------------------------
assert(ndims(l1)<=3,'The three input arguments must be 3-D volumes');
[M,N,P] = size(l1);
assert(isequal(size(l2),[M,N,P]),'Input 2 must be the same size as input 1');
assert(isequal(size(l3),[M,N,P]),'Input 3 must be the same size as input 1');
% -------------------------------------------------------------------------
opt.scalar = 'fa';      optchk.scalar = [true,false]; % variable length string
opt.mask = true(M,N,P); optchk.mask   = [true,true];  % boolean with the size of the image field
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Unroll to work comfortably:
l1 = l1(opt.mask);
l2 = l2(opt.mask);
l3 = l3(opt.mask);
% -------------------------------------------------------------------------
% Actually compute:
switch(lower(opt.scalar))
    case 'fa'
        num  = (l1-l2).*(l1-l2) + (l2-l3).*(l2-l3) + (l3-l1).*(l3-l1);
        den  = l1.*l1 + l2.*l2 + l3.*l3;
        vol2 = sqrt(num./den)/sqrt(2);
    case 'md'
        vol2 = (l1+l2+l3)/3;
    case 'ra'
        num  = (l1-l2).*(l1-l2) + (l2-l3).*(l2-l3) + (l3-l1).*(l3-l1);
        den  = l1+l2+l3;
        vol2 = (sqrt(num)./den)/sqrt(2);
    case 'cl'
        vol2 = (l1-l2)./l1;
    case 'cp'
        vol2 = (l2-l3)./l1;
    case 'cs'
        vol2 = l3./l1;
    case 'clsq'
        den  = sqrt(l1.*l1 + l2.*l2 + l3.*l3);
        vol2 = (l1-l2)./den;
    case 'cpsq'
        den  = sqrt(l1.*l1 + l2.*l2 + l3.*l3);
        vol2 = sqrt(2)*(l2-l3)./den;
    case 'cssq'
        den  = sqrt(l1.*l1 + l2.*l2 + l3.*l3);
        vol2 = sqrt(3)*l3./den;
    otherwise
        error(['Unknown scalar measure type: ',opt.scalar]);
end
% -------------------------------------------------------------------------
% Re-arrange the output:
vol = zeros(M,N,P);
vol(opt.mask) = vol2;
end
