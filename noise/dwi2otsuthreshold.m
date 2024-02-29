function [th,mask] = dwi2otsuthreshold(dwi,bi,varargin)
% function [th,mask] = dwi2otsuthreshold(dwi,bi,'opt1',val1,'opt2',val2...)
%
%    Computes a threshold that divides the baseline images of a dwi into
%    background/ foreground.
%
%       dwi: MxNxPx(B+G), the DWI image (NOT the attenuation signal, it has 
%            to include the baselines themselves).
%       bi: (B+G)x1, the b-values, including the zeros corresponding to the
%            baseline images.
%       th: 1x1, the threshold.
%       mask: MxNxP, a volume of logicals with the values where the
%            baseline is above the threshold.
%
%    Optional arguments provided as option/value pairs:
%
%       bins: histogram bins to compute statistics (default: 2048)
%       allchnl: wether (true) or not (false) average all DWI channels in
%          the volume instead of just the baselines (default: true)
%       kern: use a separable kernel to filter the input pixels and
%          smooth the target image (default: [1;1;1]/3, meaning a 3-D
%          kernel ones(3,3,3)/27 is actually used)
%       nobck: if set "true", the just-zero-voxels are ignored from
%          computations as long as they are assumed to belong to the
%          artificially removed background (default: true)

if(nargin<2)
    error('At lest the dwi volume and the b-values must be supplied');
end
[~,~,~,G] = size(dwi);
if(size(bi,1)~=G)
    error('The number of b-values bi must match the 4-th dimension of dwi');
end
if(size(bi,2)~=1)
    error('The b-values vector must be a column vector');
end

% Parse the optional input arguments:
opt.bins = 2048;       optchk.bins = [true,true];    % always 1x1 double
opt.allchnl = true;    optchk.allchnl = [true,true]; % always 1x1 logical
opt.kern = [1;1;1]/3;  optchk.kern = [true,false];   % size is variable
opt.nobck = true;      optchk.nobck = [true,true];   % always 1x1 logical
opt = custom_parse_inputs(opt,optchk,varargin{:});

if(opt.allchnl)
    dwi = mean(dwi,4);
else
    idx = (abs(bi)<1); % Gx1
    dwi = mean(dwi(:,:,:,idx),4); % MxNxP
end

if(numel(opt.kern)>1)
    krn = opt.kern(:)./sum(opt.kern(:));
    dwi = convn(dwi,krn,'same');
    dwi = convn(dwi,krn','same');
    dwi = convn(dwi,permute(krn,[2,3,1]),'same');
end

if(opt.nobck)
    [counts,pos] = hist(dwi(dwi>1.0e-6),opt.bins);
else
    [counts,pos] = hist(dwi(:),opt.bins);
end

freqs = counts/sum(counts);

frqleft  = cumsum(freqs);
muleft   = cumsum(freqs.*pos)./frqleft;
frqleft  = [0,frqleft];
muleft   = [0,muleft];

frqright = cumsum(freqs(end:-1:1));
muright  = cumsum(freqs(end:-1:1).*pos(end:-1:1))./frqright;
frqright = [frqright(end:-1:1),0];
muright  = [muright(end:-1:1),0];

varbetween = frqleft.*frqright.*(muleft-muright).*(muleft-muright);

[~,Mpos] = max(varbetween(2:end-1));
th   = (pos(Mpos)+pos(Mpos+1))/2;
mask = (dwi>th);




