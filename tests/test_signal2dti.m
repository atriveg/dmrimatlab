% test_signal2dti
load test_data.mat;
%%% -----------------------------------------------------------------------
bmin = 0;
bmax = inf;
pp   = ( (bi<=bmax) & (bi>=bmin) );
atti = atti(:,:,:,pp);
bi   = bi(pp);
gi   = gi(pp,:);
%%% -----------------------------------------------------------------------
tl = 1.0e-5;
tu = 1-tl;
atti(atti<tl) = tl;
atti(atti>tu) = tu;
%%% -----------------------------------------------------------------------
signal = -log(atti); % M x N x P x G
signal = bsxfun(@(x,y)(x./y),signal,permute(bi,[2,3,4,1]));
%%% -----------------------------------------------------------------------
tic;
tensor1 = signal2dti( signal, gi, [], 'mask', mask, 'wls', false );
T = toc;
fprintf(1,'It took %f seconds to compute the LS scheme\n',T);
NI = 5;
tic;
tensor2 = signal2dti( signal, gi, bi, 'mask', mask, 'wls', true, 'wlsit', NI, 'wsc', 1.0e-6 );
T = toc;
fprintf(1,'It took %f seconds to compute the WLS scheme\n',T);
tic;
sh      = atti2shadc( atti, gi, bi, 'mask', mask, 'L', 2, 'lambda', 0.001, 'tl', tl, 'tu', tu );
tensor3 = shadc2dti( sh,'mask', mask, 'unroll', false );
T = toc;
fprintf(1,'It took %f seconds to compute the SH-based scheme\n',T);
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

close(figure(1));
hf = figure(1);
subplot(3,1,1);
imshow(RGB1); title('Regular LS');
subplot(3,1,2);
imshow(RGB2); title('Weighted LS');
subplot(3,1,3);
imshow(RGB3); title('SH-based');
%%% -----------------------------------------------------------------------
