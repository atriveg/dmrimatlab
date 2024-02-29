% estimate_nonfree_matter
clear;
close('all');

% This script illustrates a simple pipeline to estimate the percentage of
% confined water (i.e. intracellular) in the white matter starting from
% two nrrd files (3D-Slicer format) of the same volume acquired with two
% different b-values:

% Choose a valid 3D-Slicer nrrd file:
filename1 = '/home/atriveg/Documentos/DMRIDATA/GregorioMaranon/P6_b800_r2_dwi.nrrd';
filename2 = '/home/atriveg/Documentos/DMRIDATA/GregorioMaranon/P6_b1000_r2_dwi.nrrd';
filename3 = '/home/atriveg/Documentos/DMRIDATA/GregorioMaranon/P6_b1300_r2_dwi.nrrd';

% Use 3rd party function from matlab bridge to load the volume in 3D-Slicer
% format:
vol1 = nrrdread(filename1);
vol2 = nrrdread(filename2);
vol3 = nrrdread(filename3);

% Convert 3D-Slicer-like volume to a data set to be processed with the
% toolbox:
[dwi1,giall1,biall1] = slicervol2dwi(vol1); % giall and biall include baselines' gi=[0,0,0] and bi=0
[dwi2,giall2,biall2] = slicervol2dwi(vol2); % giall and biall include baselines' gi=[0,0,0] and bi=0
[dwi3,giall3,biall3] = slicervol2dwi(vol3); % giall and biall include baselines' gi=[0,0,0] and bi=0

x = 23:103;
y = 22:113;
z = [21,26,29,31,32,33,34,35,37,40,45];
dwi1 = dwi1(x,y,z,:);
dwi2 = dwi2(x,y,z,:);
dwi3 = dwi3(x,y,z,:);

% Estimate the noise (sigma) together with a background mask computed with
% Otsu's method. Since the two volumes are acquired in one session, there
% is no need to do it twice:
[sigma,mask] = dwi2noisestd(dwi1,biall1,'estrad',[1;1;1]);
fprintf(1,'Estimated sigma: %1.2f\n',sigma);

% Select the first 16 gradients from volume 1, the remining 45 from the
% second. All them from the third
idx1   = 1:62;%1:17;
dwi1   = dwi1(:,:,:,idx1);
giall1 = giall1(idx1,:);
biall1 = biall1(idx1,:);
idx2   = 1:62;%[1,18:62];
dwi2   = dwi2(:,:,:,idx2);
giall2 = giall2(idx2,:);
biall2 = biall2(idx2,:);
idx3   = 1:62;%1:62;
dwi3   = dwi3(:,:,:,idx3);
giall3 = giall3(idx3,:);
biall3 = biall3(idx3,:);

% Clean the volumes using the denoising filters. Use the fastest approach:
% fprintf(1,'Cleaning the DWI volumes to improve the PSNR. May take a while... ');
% dwi1 = dwi2cleandwi( dwi1, giall1, sigma, 'beta', 2.0, 'mask', mask, 'onlyUNLM', true );
% dwi2 = dwi2cleandwi( dwi2, giall2, sigma, 'beta', 2.0, 'mask', mask, 'onlyUNLM', true );
% dwi3 = dwi2cleandwi( dwi3, giall3, sigma, 'beta', 2.0, 'mask', mask, 'onlyUNLM', true );
% fprintf(1,'done\n');

% Merge the two DWI volumes:
dwi1(:,:,:,idx1(end)+idx2)  = dwi2;
dwi1(:,:,:,idx1(end)+idx2(end)+idx3) = dwi3;
gi = [ giall1; giall2; giall3 ];
bi = [ biall1; biall2; biall3 ];

% Compute an attenuation signal from the DWI:
[atti,gi,bi] = dwi2atti(dwi1,gi,bi);

% Compute the percentage of non-free water
fprintf(1,'Estimating non-free water. May take a while...\n');
[f,SH] = atti2freewater( atti, gi, bi, 'L', 2, 'lambda', 0.001, 'mask', mask, 'verbose', true, 'fdepth', 4, 'fbins', 2 );
%[f,SH] = atti2freewaterTensor( atti, gi, bi, 'mask', mask, 'verbose', ...
%                            true, 'O2', false, 'O3', false ); % DOESN'T WORK!!!
f = permute(f,[2,1,3]);
f = 1-f(end:-1:1,:,:);
figure(1);
imagesc([f(:,:,3),f(:,:,4),f(:,:,5);f(:,:,6),f(:,:,7),f(:,:,8)]);
colormap(jet);
colorbar;

mu = 5.0e-6; % 0.00015; % 5.0e-5
[lpar,lperp,f2] = atti2micro( atti, gi, bi, 'tl', 1.0e-6, 'tu', 1-1.0e-6, ...
    'ADC0', 3.0e-3, 'usef', true, 'mask', mask, 'bth', 100, 'mlperp',0.01e-3, ...
    'mu', mu, 'verbose', true );
f2 = permute(f2,[2,1,3]);
f2 = 1-f2(end:-1:1,:,:);
figure(2);
imagesc([f2(:,:,3),f2(:,:,4),f2(:,:,5);f2(:,:,6),f2(:,:,7),f2(:,:,8)]);
colormap(jet);
colorbar;

lpa = permute(lpar,[2,1,3]);
lpa = lpa(end:-1:1,:,:);
figure(3);
imagesc([lpa(:,:,3),lpa(:,:,4),lpa(:,:,5);lpa(:,:,6),lpa(:,:,7),lpa(:,:,8)],[0,3.0e-3]);
colormap(jet);
colorbar;

lpp = permute(lperp,[2,1,3]);
lpp = lpp(end:-1:1,:,:);
figure(4);
imagesc([lpp(:,:,3),lpp(:,:,4),lpp(:,:,5);lpp(:,:,6),lpp(:,:,7),lpp(:,:,8)],[0,3.0e-4]);
colormap(jet);
colorbar;
