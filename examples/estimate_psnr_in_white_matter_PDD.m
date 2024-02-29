% estimate_psnr_in_white_matter
clear;

% This script illustrates a simple pipeline to estimate the peak SNR (i.e.
% the value of the baseline over sigma, where sigma is the standard
% deviation of noise in the complex domin) in the white matter of a brain
% volume, starting from reading a nrrd file, masking the background,
% estimating the GFA to segment the white matter, estiamting sigma and
% plotting the hisotgram of the PSNR.

% Choose a valid 3D-Slicer nrrd file:
%filename = '/Users/atriveg/Downloads/P6_b800_r2_dwi.nrrd';
PATH = '/media/atriveg/DATOS/Documentos/DMRIDATA/PDD/CO_p07090/';
filename = 'CO_p07090_dwi.nii.gz';
bvalname = 'CO_p07090.bval';
bvecname = 'CO_p07090.bvec';

% Use 3rd party function from matlab bridge to load the volume in 3D-Slicer
% format:
vol = load_untouch_nii([PATH,filename]);

% Convert 3D-Slicer-like volume to a data set to be processed with the
% toolbox:
dwi    = double(vol.img);
biall  = load([PATH,bvalname]);
giall  = load([PATH,bvecname])';

% Only dwis with b=1000:
%pt  = [1,122,2:121];
% Only dwis with b=2500;
%pt    = [1,122:242];
% All dwis:
pt = [1,122,2:121,123:242];
dwi   = dwi(:,:,:,pt);
biall = biall(pt,1);
giall = giall(pt,1:3);

% Estimate the noise (sigma) together with a background mask computed with
% Otsu's method:
[sigma,mask] = dwi2noisestd(dwi,biall,'nograd',true,'allchnl',true);

% Obtain an attenuation signal from the DWI images and the baselines, from
% the 3D-Slicer-like volume to Matlab arrays:
baseline = mean(dwi(:,:,:,1:2),4);
atti     = bsxfun( @(x,y)(x./y), dwi(:,:,:,3:end), baseline );
gi       = giall(3:end,1:3);
bi       = biall(3:end,1);

% Compute the Spherical Hamonics expansion of the Apparent Diffusion
% Coefficient from the attenuation signal:
SH = atti2shadc( atti, gi, bi, 'mask', mask, 'L', 8, 'lambda', 0.001 );

% Compute the GFA from the SH coefficients:
GFA = sum(SH(:,:,:,2:end-1).*SH(:,:,:,2:end-1),4)./sum(SH.*SH,4);
GFA = sqrt(GFA);

% Roughly segment the white matter by thresholding the GFA:
th    = 0.15;
mask2 = ( mask & (GFA>=th) );

% Find the baselines within the dwi volume and average them:
bidx     = find(biall<1);
baseline = mean(dwi(:,:,:,bidx),4);

% Compute the PSNR in the selected voxels:
psnr = baseline(mask2)/sigma;

% Compute and plot the statistics:
[cc,pp] = hist(psnr,100);
cc  = cc./sum(cc);
ccs = cumsum(cc);
hf = figure(1001);
close(hf);
hf = figure(1001);
hold('on');
plot(pp,cc,'k','LineWidth',2);
p10 = find(ccs>0.1,1,'first');
p90 = find(ccs>0.9,1,'first');
plot([pp(p10),pp(p10)],[0,max(cc)],'b--');
plot([pp(p90),pp(p90)],[0,max(cc)],'b--');
grid('on');
title(sprintf('sigma=%1.1f, mean PSNR=%1.1f',sigma,mean(psnr)),'FontSize',18);
xlabel('PSNR','FontSize',16);
ylabel('frequency','FontSize',16);
set(get(hf,'CurrentAxes'),'FontSize',14);

hf=figure(2002);
close(hf);
figure(2002);
subplot(2,2,1); imshow(baseline(:,:,27),[]);
subplot(2,2,2); imshow(mask(:,:,27));
subplot(2,2,3); imshow(GFA(:,:,27),[]);
subplot(2,2,4); imshow(mask2(:,:,27));
























