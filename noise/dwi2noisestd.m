function [noiseStd,mask] = dwi2noisestd(dwi,bi,varargin)
% function [noise,mask] = dwi2noisestd(dwi,bi,'opt1',val1,'opt2',val2,...)
%
%    Computes the standard deviation of noise, noiseStd, in a stationary
%    Rician distributed DWI image, where noiseStd is the standard deviation
%    of Gaussian noise in the original complex domain:
%
%        A_noisy = sqrt( (A_c+sigmaStd*randn) + 1i*(A_s+sigmaStd*randn) )
%
%       dwi: MxNxPx(B+G), the DWI image (NOT the attenuation signal, it has 
%            to include the baselines themselves).
%       bi: (B+G)x1, the b-values, including the zeros corresponding to the
%            baseline images.
%       noise: 1x1, the standard decviarion of noise in the complex
%            domain.
%       mask: MxNxP, the mask (logicals) computed with Otsu's method used
%            to skip computations in the background of the image.
%
%    The noise is estimated from the mode (histogram-based) of the local
%    variance of the DWI within the image backgound, which is oughly
%    segmented based on Otsu's method
%
%    Optional arguments provided as option/value pairs, related to Otsu's
%    method to segment teh background:
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
%
%    Optional arguments provided as option/value pairs, related to the
%    estimation of noise as the mode of the local variance:
%
%       nograd: if true, the noise is estimated based on statistics
%          computed just over the baseline (or baselines) of the DWI. If
%          false, the statistics are computed over the entire DWI volume,
%          including the gadient images (default: false)
%       estrad: estimation radii to compute local variances (default:
%          [1;1;0])
%       ebins: number of histogram bins to compute the mode of the local
%          variance in the image backgound (defaul: 1024)
%

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
opt.nograd = false;    optchk.nograd = [true,true];  % always 1x1 logical
opt.estrad = [1;1;0];  optchk.estrad = [true,true];  % always 3x1 double
opt.ebins = 1024;      optchk.ebins = [true,true];   % always 1x1 double
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Compute the background mask:
[~,mask] = dwi2otsuthreshold(dwi,bi,...
    'bins',opt.bins,'allchnl',opt.allchnl,...
    'kern',opt.kern,'nobck',opt.nobck);

% If the statistics will be computed just over the baseline images, there
% is no need to keep the whole volume here:
if(opt.nograd)
    bsidx = find(bi<1);
    if(~isempty(bsidx))
        dwi = dwi(:,:,:,bsidx);
    end
end

% Replicate the mask to 4-D
mask  = repmat(mask,[1,1,1,size(dwi,4)]);

% Compute the local mean:
lmean = compute_local_directional_mean(dwi,  1,opt.estrad(1));
lmean = compute_local_directional_mean(lmean,2,opt.estrad(2));
lmean = compute_local_directional_mean(lmean,3,opt.estrad(3));

% Compute the local variance:
lvar  = dwi.*dwi;
lvar  = compute_local_directional_mean(lvar,1,opt.estrad(1));
lvar  = compute_local_directional_mean(lvar,2,opt.estrad(2));
lvar  = compute_local_directional_mean(lvar,3,opt.estrad(3));

% Update the mask to skip artificial background voxels:
mask  = mask & (lmean>1.0e-2);

% Compute actual statistics:
N = prod(2*opt.estrad+1);
lvar(mask)   = lvar(mask)-(lmean(mask).*lmean(mask))/N;
lvar(mask)   = lvar(mask)/(N-1);
lmean(mask)  = lmean(mask)/N;
lvar(lvar>0) = sqrt(lvar(lvar>0));
lvar(lvar<0) = sqrt(-lvar(lvar<0));

% Correct the values to account for Gaussian/Rice/Rayleigh behavior:
lvar  = correct_local_noiseStd(lmean,lvar,mask);

% The overall mean of the local variance within the background:
stdMean = mean(lvar(mask));

% Compute the histogram within a proper range:
mask2 = ( mask & (lvar>stdMean/100) & (lvar<stdMean*3) );
[cc,pp]  = hist(lvar(mask2),opt.ebins);
[~,mPos] = max(cc);
noiseStd = pp(mPos);
mask     = mask(:,:,:,1);

function vol2 = compute_local_directional_mean(vol,dim,radius)

if(radius>0)
    N    = size(vol,dim);
    ptra = 1:N;
    ptrb = 1:N;
    vol2 = vol;
    for k=1:radius
        ptra = [ptra(N),ptra(1:N-1)];
        ptrb = [ptrb(2:N),ptrb(1)];
        switch(dim)
            case 1
                vol2 = vol2 + vol(ptra,:,:,:);
                vol2 = vol2 + vol(ptrb,:,:,:);
            case 2
                vol2 = vol2 + vol(:,ptra,:,:);
                vol2 = vol2 + vol(:,ptrb,:,:);
            case 3
                vol2 = vol2 + vol(:,:,ptra,:);
                vol2 = vol2 + vol(:,:,ptrb,:);
        end
    end
else
    vol2 = vol;
end


function lvar2  = correct_local_noiseStd(lmean,lvar,mask)


m1 = (lmean>2*lvar); % Almost pure Gaussian. Nothing to do...
m2 = (lvar>2*lmean); % Almost pure Rayleigh
m3 = ~(m1|m2);       % Strictly Rician
m2 = m2 & mask;
m3 = m3 & mask;

lvar2     = lvar;
lvar2(m2) = lvar(m2)/sqrt((4-pi)/2);
sigma     = 1/2 + 2.5528*2*(lvar(m3)./lmean(m3)-1/2)/3;
lvar2(m3) = sigma.*lmean(m3);
