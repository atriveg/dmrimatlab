% test_jaLMMSE

disp('NOTE: This script is mainly intended to give you an idea on how');
disp('   long each filter configuration will take in your computer');
disp('   (be ready to wait for several tens of minutes),');
disp('   and also for you to have a quick glance of the results');
disp('   you may expect from each configuration');

load('test_dwi2cleandwi.mat');

DWInoisy = double(DWInoisy);

NG = 5;

tic;
DWIfiltered1 = dwi2cleandwi( DWInoisy, Grads, sigma, ...
    'rs', [2;2;2], ...
    'rc', [1;1;1], ...
    'beta', 1.5, ...
    'Ng', 0, ...
    'onlyUNLM', true, ...
    'filterOutliers', true, ...
    'mask', mask );
etime = toc;
disp(['Only UNLM version with mask and all gradients completed in ',num2str(etime),' seconds']);

tic;
DWIfiltered2 = dwi2cleandwi( DWInoisy, Grads, sigma, ...
    'rs', [2;2;2], ...
    'rc', [1;1;1], ...
    'beta', 1.5, ...
    'Ng', 0, ...
    'onlyUNLM', true, ...
    'filterOutliers', true, ...
    'mask', [] );
etime = toc;
disp(['Only UNLM version without mask and all gradients completed in ',num2str(etime),' seconds']);

disp('NOTE: With the UNLM version, the value of Ng becomes irrelevant');

tic;
DWIfiltered3 = dwi2cleandwi( DWInoisy, Grads, sigma, ...
    'rs', [2;2;2], ...
    'rc', [1;1;1], ...
    'beta', 1.5, ...
    'Ng', 0, ...
    'onlyUNLM', false, ...
    'filterOutliers', true, ...
    'mask', mask );
etime = toc;
disp(['Whole Wiener version with mask and all gradients completed in ',num2str(etime),' seconds']);

tic;
DWIfiltered4 = dwi2cleandwi( DWInoisy, Grads, sigma, ...
    'rs', [2;2;2], ...
    'rc', [1;1;1], ...
    'beta', 1.5, ...
    'Ng', NG, ...
    'onlyUNLM', false, ...
    'filterOutliers', true, ...
    'mask', mask );
etime = toc;
disp(['Whole Wiener version with mask and partial (',num2str(NG),') gradients completed in ',num2str(etime),' seconds']);

tic;
DWIfiltered5 = dwi2cleandwi( DWInoisy, Grads, sigma, ...
    'rs', [2;2;2], ...
    'rc', [1;1;1], ...
    'beta', 1.5, ...
    'Ng', 0, ...
    'onlyUNLM', false, ...
    'filterOutliers', true, ...
    'mask', [] );
etime = toc;
disp(['Whole Wiener version without mask and all gradients completed in ',num2str(etime),' seconds']);

slice = 33;
grad  = 3;

figure(1);

subplot(2,2,1);
imshow(DWInoisy(:,:,slice,grad),[]);
title('Original noisy slice');

subplot(2,2,2);
imshow(DWIfiltered1(:,:,slice,grad),[]);
title('UNLM filter without LMMSE correction');

subplot(2,2,3);
imshow(DWIfiltered3(:,:,slice,grad),[]);
title('Whole filter mixing all gradients');

subplot(2,2,4);
imshow(DWIfiltered4(:,:,slice,grad),[]);
title('Whole filter mixing only part of the gradients');