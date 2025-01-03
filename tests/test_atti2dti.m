% test_signal2dti
%%% -----------------------------------------------------------------------
close('all');
clear;
clc;
rng(13061981);
%%% -----------------------------------------------------------------------
load test_data.mat;
%%% -----------------------------------------------------------------------
bmin = 0;
bmax = 2000;
pp   = ( (bi<=bmax) & (bi>=bmin) );
atti = atti(:,:,:,pp);
bi   = bi(pp);
gi   = gi(pp,:);
%%% -----------------------------------------------------------------------
tl = 1.0e-9;
tu = 1-tl;
atti(atti<tl) = tl;
atti(atti>tu) = tu;
%%% -----------------------------------------------------------------------
signal = -log(atti); % M x N x P x G
signal = bsxfun(@(x,y)(x./y),signal,permute(bi,[2,3,4,1]));
%%% -----------------------------------------------------------------------
tic;
[tensor1,S01] = atti2dti( atti, gi, bi, ...
    'wls', false, ...
    'nonlinear', false, ...
    'wlsit', 1, ...
    'wsc', 0.01, ...
    'tol', 1.0e-6, ...
    'maxiters', 100, ...
    'rcondth', 1.0e-6, ...
    'fixmode', 'z', ...
    'maxthreads', 4, ...
    'unroll', false, ...
    'mask', mask );
T = toc;
fprintf(1,'It took %f seconds to compute the OLS scheme\n',T);
NI = 5;
% ----
tic;
[tensor2,S02] = atti2dti( atti, gi, bi, ...
    'wls', true, ...
    'nonlinear', false, ...
    'wlsit', 5, ...
    'wsc', 0.01, ...
    'tol', 1.0e-6, ...
    'maxiters', 100, ...
    'rcondth', 1.0e-6, ...
    'fixmode', 'a', ...
    'maxthreads', 4, ...
    'unroll', false, ...
    'mask', mask );
T = toc;
fprintf(1,'It took %f seconds to compute the WLS scheme\n',T);
% ----
tic;
[tensor3,S03] = atti2dti( atti, gi, bi, ...
    'wls', true, ...
    'nonlinear', true, ...
    'wlsit', 3, ...
    'wsc', 0.01, ...
    'tol', 1.0e-6, ...
    'maxiters', 100, ...
    'rcondth', 1.0e-6, ...
    'fixmode', 'a', ...
    'maxthreads', 4, ...
    'unroll', false, ...
    'mask', mask );
T = toc;
fprintf(1,'It took %f seconds to compute the NLLS scheme\n',T);
%%% -----------------------------------------------------------------------
[u11,u21,u31,l11,l21,l31] = dti2spectrum( tensor1, 'mask', mask );
[u12,u22,u32,l12,l22,l32] = dti2spectrum( tensor2, 'mask', mask );
[u13,u23,u33,l13,l23,l33] = dti2spectrum( tensor3, 'mask', mask );
%%% -----------------------------------------------------------------------
FA1 = sqrt(3/2)*sqrt( (l11-l21).^2 + (l11-l31).^2 + (l21-l31).^2 )./sqrt(l11.^2+l21.^2+l31.^2);
FA2 = sqrt(3/2)*sqrt( (l12-l22).^2 + (l12-l32).^2 + (l22-l32).^2 )./sqrt(l12.^2+l22.^2+l32.^2);
FA3 = sqrt(3/2)*sqrt( (l13-l23).^2 + (l13-l33).^2 + (l23-l33).^2 )./sqrt(l13.^2+l23.^2+l33.^2);
FA1(~mask) = 0;
FA2(~mask) = 0;
FA3(~mask) = 0;
%%% -----------------------------------------------------------------------
u11 = abs(bsxfun(@(x,y)(x.*y),u11,FA1));
u12 = abs(bsxfun(@(x,y)(x.*y),u12,FA2));
u13 = abs(bsxfun(@(x,y)(x.*y),u13,FA3));
%%% -----------------------------------------------------------------------
sl = [4,7,9,11,14];
RGB1 = [ ...
    permute(u11(:,:,sl(1),:),[2,1,4,3]), ...
    permute(u11(:,:,sl(2),:),[2,1,4,3]), ...
    permute(u11(:,:,sl(3),:),[2,1,4,3]), ...
    permute(u11(:,:,sl(4),:),[2,1,4,3]), ...
    permute(u11(:,:,sl(5),:),[2,1,4,3]) ];
RGB2 = [ ...
    permute(u12(:,:,sl(1),:),[2,1,4,3]), ...
    permute(u12(:,:,sl(2),:),[2,1,4,3]), ...
    permute(u12(:,:,sl(3),:),[2,1,4,3]), ...
    permute(u12(:,:,sl(4),:),[2,1,4,3]), ...
    permute(u12(:,:,sl(5),:),[2,1,4,3]) ];
RGB3 = [ ...
    permute(u13(:,:,sl(1),:),[2,1,4,3]), ...
    permute(u13(:,:,sl(2),:),[2,1,4,3]), ...
    permute(u13(:,:,sl(3),:),[2,1,4,3]), ...
    permute(u13(:,:,sl(4),:),[2,1,4,3]), ...
    permute(u13(:,:,sl(5),:),[2,1,4,3]) ];
%%% -----------------------------------------------------------------------
close(figure(1));
hf1 = figure(1);
subplot(3,1,1);
imshow(RGB1); title('Regular LS');
subplot(3,1,2);
imshow(RGB2); title('Weighted LS');
subplot(3,1,3);
imshow(RGB3); title('Nonlinear LS');
%%% -----------------------------------------------------------------------
win1 = quantile( S01(mask), [.01,.99] );
win2 = quantile( S02(mask), [.01,.99] );
win3 = quantile( S03(mask), [.01,.99] );
%%% -----------------------------------------------------------------------
S01 = [ ...
    S01(:,:,sl(1))', ...
    S01(:,:,sl(2))', ...
    S01(:,:,sl(3))', ...
    S01(:,:,sl(4))', ...
    S01(:,:,sl(5))' ...
    ];
S02 = [ ...
    S02(:,:,sl(1))', ...
    S02(:,:,sl(2))', ...
    S02(:,:,sl(3))', ...
    S02(:,:,sl(4))', ...
    S02(:,:,sl(5))' ...
    ];
S03 = [ ...
    S03(:,:,sl(1))', ...
    S03(:,:,sl(2))', ...
    S03(:,:,sl(3))', ...
    S03(:,:,sl(4))', ...
    S03(:,:,sl(5))' ...
    ];
close(figure(2));
hf2 = figure(2);
subplot(3,1,1);
imshow(S01,win1); colormap(jet); colorbar; title('Regular LS');
subplot(3,1,2);
imshow(S02,win2); colormap(jet); colorbar; title('Weighted LS');
subplot(3,1,3);
imshow(S03,win3); colormap(jet); colorbar; title('Nonlinear LS');
%%% -----------------------------------------------------------------------
