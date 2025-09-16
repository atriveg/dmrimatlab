%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_sh2hot2sh.m
rng(13061981);
clear;
close('all');
clc;
format('compact');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('test_data.mat');
tl = 100*eps;
tu = 1.0-tl;
atti(atti<tl) = tl;
atti(atti>tu) = tu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xc = 25:76;
yc = 24:78;
zc = 8:10;
atti = atti(xc,yc,zc,:);
mask = mask(xc,yc,zc);
lattice = [1;1;2];
maxiters = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH2 = atti2shadc( atti, gi, bi, 'L', 2, 'mask', mask, 'lambda', 0.001, ...
    'tl', tl, 'tu', tu );
GFA = sqrt( sum(SH2(:,:,:,2:end).*SH2(:,:,:,2:end),4)./sum(SH2.*SH2,4) );
GFA = GFA.*double(mask);
th  = quantile(GFA(mask),0.9);
pp  = find(GFA>=th);
pp  = pp(randperm(length(pp),1));
mask2 = false(size(mask));
mask2(pp) = true;
% mask = mask2; %%% BEWARE!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[eap2,dti2,lattice2,Qx2,Qy2,Qz2,res2,lapl2] = atti2hydidsi( atti, gi, bi, ...
    'lambda', 1.0e-3, ...
    'tl', tl, ...
    'tu', tu, ...
    'ADC0', 3.0e-3, ...
    'tau', 35.0e-3, ...
    'Rth', 0.01, ...
    'lattice', lattice, ...
    'const', true, ...
    'miters', maxiters, ...
    'otol', 1.0e-6, ...
    'stol', 1.0e-6, ...
    'ctol', 1.0e-8, ...
    'mask', mask, ...
    'usemex', true, ...
    'maxthreads', 32 );
fprintf(1,'It took %1.5f seconds to compute with mex code\n',toc);
tic;
[eap,dti,lattice,Qx,Qy,Qz,res,lapl] = atti2hydidsi( atti, gi, bi, ...
    'lambda', 1.0e-3, ...
    'tl', tl, ...
    'tu', tu, ...
    'ADC0', 3.0e-3, ...
    'tau', 35.0e-3, ...
    'Rth', 0.01, ...
    'lattice', lattice, ...
    'const', true, ...
    'miters', maxiters, ...
    'otol', 1.0e-6, ...
    'stol', 1.0e-6, ...
    'ctol', 1.0e-8, ...
    'mask', mask, ...
    'usemex', false, ...
    'maxthreads', 32 );
fprintf(1,'It took %1.5f seconds to compute with matlab code\n',toc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(1));
hf1 = figure(1);
subplot(1,3,1);
imshow([Qx(:,:,1)',Qx(:,:,2)',Qx(:,:,3)';Qx2(:,:,1)',Qx2(:,:,2)',Qx2(:,:,3)'],[]);
colorbar;
colormap(jet);
subplot(1,3,2);
imshow([Qy(:,:,1)',Qy(:,:,2)',Qy(:,:,3)';Qy2(:,:,1)',Qy2(:,:,2)',Qy2(:,:,3)'],[]);
colorbar;
colormap(jet);
subplot(1,3,3);
imshow([Qz(:,:,1)',Qz(:,:,2)',Qz(:,:,3)';Qz2(:,:,1)',Qz2(:,:,2)',Qz2(:,:,3)'],[]);
colorbar;
colormap(jet);
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(2));
hf2 = figure(2);
subplot(1,2,1);
imshow([res(:,:,1)',res(:,:,2)',res(:,:,3)';res2(:,:,1)',res2(:,:,2)',res2(:,:,3)'],[]);
colorbar;
colormap(jet);
subplot(1,2,2);
imshow([lapl(:,:,1)',lapl(:,:,2)',lapl(:,:,3)';lapl2(:,:,1)',lapl2(:,:,2)',lapl2(:,:,3)'],[]);
colorbar;
colormap(jet);
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(3));
hf3 = figure(3);
imshow([eap(:,:,1,1)',eap(:,:,2,1)',eap(:,:,3,1)';eap2(:,:,1,1)',eap2(:,:,2,1)',eap2(:,:,3,1)'],[]);
colorbar;
colormap(jet);
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
