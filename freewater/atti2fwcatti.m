function [atti,gi,bi] = atti2fwcatti( atti, gi, bi, f, varargin )
% function [attic,gi,bi] = atti2fwcatti( atti, gi, bi, f,
%                                    'opt1', value1, 'opt2', value2, ... )
%
%   Corrects the attenuation signal to remove its free-water compartment
%   as computed by atti2freewater, according to the signal model:
%
%     Si(g,b) = f*\int\int_{u}\Phi(u)*k(g,b;u) + (1-f)*exp(-b*D0)
%             = f*S_c(g,b) + (1-f)*S_f(b)
%
%   i.e. finds the partial volume of water with constrained diffusion, S_c,
%   from the raw attenuation signal Si and the partial volume fraction f 
%   (as returned by atti2freewater). For further details, see:
%
%       Antonio Tristan-Vega; Guillem Paris; Rodrigo de Luis-Garcia;
%       Santiago Aja-Fernandez. "Accurate free-water estimation in white
%       matter from fast diffusion MRI acquisitions using the spherical
%       means technique". Magnetic Resonance in Medicine 87(2),
%       pp. 1028–1035. Wiley, 2022.
%
%   MANDATORY INPUTS:
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradient over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      bi: a Gx1 vector with the b-values used at each gradient direction,
%         so that multi-shell acquisitions are allowed.
%      f: a MxNxP aray in the range (0,1) with the partial volume fraction
%         of non-free, i.e. confined, water (the first output returned by
%         atti2freewater).
%
%   OUPUTS:
%
%      attic: a MxNxPxG double array containing the corrected signal S_c
%         for the very same values of gi and bi as the original.
%      gi: same as the input, returned for convenience.
%      bi: same as the input, returned for convenience.
%
%   OPTIONAL arguments may be passed as name/value pairs in the regular
%   Matlab style:
%
%      nonan: whether (true) or not (false) remove nans from the corrected
%         atti signal at those voxels where f becomes 0 (default: true).
%      clip: whether (true) or not (false) clip out of bounds values of
%         the corrected atti signal to the allowed range [tl,tu]
%         (default: false).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the atti will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-7, 1-1.0e-7).
%      ADC0: estimated diffusivity of free water at body temperature. Do
%         not change this value unless you have a good reason to do so
%         (default: 3.0e-3).
%
%   Other general options:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.

% Check the mandatory input argments:
if(nargin<4)
    error('At least the atti volume, the gradients table, and the b-values must be supplied');
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
    error('The number of b-values bi must match the number of entries of gi');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end
assert( isequal(size(f),[M,N,P]), 'The FoV of f does not match that of atti' );
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.nonan = true;       optchk.nu = [true,true];
opt.clip = false;       optchk.nu = [true,true];
opt.tl = 1.0e-7;        optchk.tl = [true,true];       % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];       % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
atti = reshape(atti,[M*N*P,G]);
mask = opt.mask(:);
raw  = atti(mask,:);
raw(raw<opt.tl) = opt.tl;
raw(raw>opt.tu) = opt.tu;
f    = f(:);
f    = f(mask);
f(f<0) = 0;
f(f>1) = 1;
% -------------------------------------------------------------------------
raw  = ( raw - (1-f).*exp(-(bi')*opt.ADC0) );
raw(f~=0,:) = raw(f~=0,:)./f(f~=0);
raw(f==0,:) = nan;
if(opt.nonan)
    raw(isnan(raw)) = 0;
end
if(opt.clip)
    raw(raw<opt.tl) = opt.tl;
    raw(raw>opt.tu) = opt.tu;
end
% -------------------------------------------------------------------------
atti(mask,:) = raw;
atti = reshape(atti,[M,N,P,G]);
% -------------------------------------------------------------------------
end
