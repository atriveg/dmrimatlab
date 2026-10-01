function test_dki

sf = check_software_platform;
if(sf==2)
    pkg load optim;
end

% -------------------------------------------------------------------------
close('all');
S    = load('test_data','atti','mask','gi','bi');
atti = S.atti;
mask = S.mask;
gi   = S.gi;
bi   = S.bi;
% -------------------------------------------------------------------------
wls    = true;
wlsit  = 5;
wsc    = 0.01;
qp     = true;
nconst = 60;
estS0  = false;
negkrt = false;
% -------------------------------------------------------------------------
pdti = (bi<=1500);
[dti_,S0_] = atti2dti( atti(:,:,:,pdti), gi(pdti,:), bi(pdti,1), ...
    'mask', mask, 'nonlinear', true );
[dkt,dti,S0] = atti2dkt( atti, gi, bi, [], [], ...
    'mask', mask, 'wls', wls, 'wlsit', wlsit, 'wsc', wsc, 'qp', qp, 'nconst', nconst, ...
    'estimS0', estS0, 'negkrt', negkrt );
[mk,ak,rk] = dkt2kurtosis( dkt, dti, 'mask', mask );
atti2 = dkt2atti( dkt, dti, S0, gi, bi, 'mask', mask );
atti3 = dkt2atti( zeros(size(dkt)), dti, S0, gi, bi, 'mask', mask );
% -------------------------------------------------------------------------
G   = size(atti,4);
% ----------------------
tv1 = [23,55,4];
tv2 = [74,23,7];
tv3 = [45,36,13];
% ----------------------
opts.mode    = 'q';
opts.wlsit   = wlsit;
opts.wsc     = wsc;
opts.nconst  = nconst;
opts.negkrt  = negkrt;
opts.estimS0 = estS0;
[Ain,ADc] = create_reconst_matrixes_dki(nconst,estS0);
% ----------------------
sol1M  = atti2dki_voxel( double(reshape(atti(tv1(1),tv1(2),tv1(3),:),[G,1])), double(gi), double(bi), Ain, ADc, opts );
sol1MK = sol1M(1:15);
sol1MD = sol1M(16:21);
if(opts.estimS0)
    sol1MS = sol1M(22);
else
    sol1MS = 1;
end
sol2M  = atti2dki_voxel( double(reshape(atti(tv2(1),tv2(2),tv2(3),:),[G,1])), double(gi), double(bi), Ain, ADc, opts );
sol2MK = sol2M(1:15);
sol2MD = sol2M(16:21);
if(opts.estimS0)
    sol2MS = sol2M(22);
else
    sol2MS = 1;
end
sol3M  = atti2dki_voxel( double(reshape(atti(tv3(1),tv3(2),tv3(3),:),[G,1])), double(gi), double(bi), Ain, ADc, opts );
sol3MK = sol3M(1:15);
sol3MD = sol3M(16:21);
if(opts.estimS0)
    sol3MS = sol3M(22);
else
    sol3MS = 1;
end
% -------
sol1XK = squeeze(dkt(tv1(1),tv1(2),tv1(3),:));
sol2XK = squeeze(dkt(tv2(1),tv2(2),tv2(3),:));
sol3XK = squeeze(dkt(tv3(1),tv3(2),tv3(3),:));
sol1XD = squeeze(dti(tv1(1),tv1(2),tv1(3),:));
sol2XD = squeeze(dti(tv2(1),tv2(2),tv2(3),:));
sol3XD = squeeze(dti(tv3(1),tv3(2),tv3(3),:));
sol1XS = squeeze(S0(tv1(1),tv1(2),tv1(3)));
sol2XS = squeeze(S0(tv2(1),tv2(2),tv2(3)));
sol3XS = squeeze(S0(tv3(1),tv3(2),tv3(3)));
% -------
[sol1MD,sol1XD], %#ok<NOPRT>
[sol1MK,sol1XK], %#ok<NOPRT>
[sol1MS,sol1XS], %#ok<NOPRT>
[sol2MD,sol2XD], %#ok<NOPRT>
[sol2MK,sol2XK], %#ok<NOPRT>
[sol2MS,sol2XS], %#ok<NOPRT>
[sol3MD,sol3XD], %#ok<NOPRT>
[sol3MK,sol3XK], %#ok<NOPRT>
[sol3MS,sol3XS], %#ok<NOPRT>
% -------------------------------------------------------------------------
close(figure(1));
figure(1);
ch = [23,34,67,81,134,167];
for c=1:length(ch)
    % ----
    subplot(3,length(ch),c+0*length(ch));
    vol = atti(:,:,:,ch(c)).*mask;
    imshow( vol(:,:,9)', [0,1] );
    vls = vol(mask);
    vls = vls( ~isnan(vls) & ~isinf(vls) );
    title( sprintf('b=%1.0f, mu=%1.3f', bi(ch(c)), mean(vls) ) );
    % ----
    subplot(3,length(ch),c+1*length(ch));
    vol = atti3(:,:,:,ch(c));
    imshow( vol(:,:,9)', [0,1] );
    vls = vol(mask);
    vls = vls( ~isnan(vls) & ~isinf(vls) );
    title( sprintf('b=%1.0f, mu=%1.3f', bi(ch(c)), mean(vls) ) );
    % ----
    subplot(3,length(ch),c+2*length(ch));
    vol = atti2(:,:,:,ch(c));
    imshow( vol(:,:,9)', [0,1] );
    vls = vol(mask);
    vls = vls( ~isnan(vls) & ~isinf(vls) );
    title( sprintf('b=%1.0f, mu=%1.3f', bi(ch(c)), mean(vls) ) );
    % ----
end
drawnow;
% -------------------------------------------------------------------------
[u1_,~,~,l1_,l2_,l3_] = dti2spectrum( dti_, 'mask', mask );
[u1,~,~,l1,l2,l3]     = dti2spectrum( dti,  'mask', mask );
md_ = spectrum2scalar(l1_,l2_,l3_,'mask',mask,'scalar','md');
md  = spectrum2scalar(l1,l2,l3,'mask',mask,'scalar','md');
rgb_ = spectrum2colorcode(u1_,l1_,l2_,l3_,'mask',mask);
rgb  = spectrum2colorcode(u1,l1,l2,l3,'mask',mask);
% -------------------------------------------------------------------------
close(figure(2));
figure(2);
sl = [2,5,8,11,14,16];
for s=1:length(sl)
    % ----
    subplot(2,length(sl),s+0*length(sl));
    IMG = S0_(:,:,sl(s))';
    imshow(IMG,[]);
    colormap(parula);
    colorbar;
    title('S_0');
    % ----
    subplot(2,length(sl),s+1*length(sl));
    IMG = S0(:,:,sl(s))';
    imshow(IMG,[]);
    colormap(parula);
    colorbar;
    title('S_0');
end
drawnow;
% -------------------------------------------------------------------------
close(figure(3));
figure(3);
sl = [2,5,8,11,14,16];
for s=1:length(sl)
    % ----
    subplot(2,length(sl),s+0*length(sl));
    IMG = md_(:,:,sl(s))';
    imshow(IMG,[0,3.0e-3]);
    colormap(parula);
    colorbar;
    title('MD');
    % ----
    subplot(2,length(sl),s+1*length(sl));
    IMG = md(:,:,sl(s))';
    imshow(IMG,[0,3.0e-3]);
    colormap(parula);
    colorbar;
    title('MD');
end
drawnow;
% -------------------------------------------------------------------------
close(figure(4));
figure(4);
sl = [2,5,8,11,14,16];
for s=1:length(sl)
    % ----
    subplot(2,length(sl),s+0*length(sl));
    IMG = squeeze(rgb_(:,:,sl(s),:));
    IMG = permute(IMG,[2,1,3]);
    imshow(IMG);
    title('color FA');
    % ----
    subplot(2,length(sl),s+1*length(sl));
    IMG = squeeze(rgb(:,:,sl(s),:));
    IMG = permute(IMG,[2,1,3]);
    imshow(IMG);
    title('color FA');
end
drawnow;
% -------------------------------------------------------------------------
close(figure(5));
figure(5);
sl = [2,5,8,11,14,16];
for s=1:length(sl)
    % ----
    subplot(3,length(sl),s+0*length(sl));
    IMG = mk(:,:,sl(s))';
    imshow(IMG,[0,1.5]);
    colormap(parula);
    colorbar;
    title('MK');
    % ----
    subplot(3,length(sl),s+1*length(sl));
    IMG = ak(:,:,sl(s))';
    imshow(IMG,[0,1.25]);
    colormap(parula);
    colorbar;
    title('K_{||}');
    % ----
    subplot(3,length(sl),s+2*length(sl));
    IMG = rk(:,:,sl(s))';
    imshow(IMG,[0,2.0]);
    colormap(parula);
    colorbar;
    title('K_{\perp}');
end
drawnow;
% -------------------------------------------------------------------------
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function x = atti2dki_voxel(Si,gi,bi,Ain,Dreconst,opts)
% -------------------------------------------------------------------------
[~,mu,pwrs] = hot2signal( zeros(1,1,1,15), [0,0,1] );
AK = (mu').*realpow(gi(:,1),pwrs(:,1)') ...
    .*realpow(gi(:,2),pwrs(:,2)').*realpow(gi(:,3),pwrs(:,3)');
AK = AK.*(bi.*bi/6);
% ------
AD = [ gi(:,1).*gi(:,1), 2*gi(:,1).*gi(:,2), 2*gi(:,1).*gi(:,3), ...
        gi(:,2).*gi(:,2), 2*gi(:,2).*gi(:,3), gi(:,3).*gi(:,3) ];
AD = -bi.*AD*1000;
% ------
A  = [AK,AD];
if(opts.estimS0)
    sc = sqrt(mean(A(:).*A(:)));
    A  = [ A, sc*ones(size(A,1),1) ];
end
% ------
Si(Si<eps) = eps;
lSi = log(Si);
% -------------------------------------------------------------------------
% OLS
x = ((A')*A)\((A')*lSi);
% -------------------------------------------------------------------------
% WLS
wls_success = true;
if( opts.mode=='w' || opts.mode=='q' )
    Aw   = A;
    lSiw = lSi;
    for n=1:opts.wlsit
        % --------------------
        wi = A*x;
        wi = exp(2*wi);
        wi(wi>1) = 1;
        wi(wi<0) = 0;
        if(any(isnan(wi)))
            break;
        end
        wi( wi < (opts.wsc)*max(wi) ) = (opts.wsc)*max(wi);
        % --------------------
        lSiw = wi.*lSi;
        Aw   = wi.*A;
        % --------------------
        x2 = ((A')*Aw)\((A')*lSiw);
        if(any(isnan(x2)))
            wls_success = false;
            break;
        else
            x = x2;
        end
    end
end
% -------------------------------------------------------------------------
% QP
if( opts.mode=='q' && wls_success )
    %------------------
    sz  = size(Ain,1)/3;
    if(opts.negkrt)
        bin = [ zeros(sz,1); Dreconst*x(16:21); zeros(sz,1) ];
        bin = 2*bin.*bin;
    else
        bin = zeros(3*sz,1);
    end
    %------------------
    Q = ((A')*Aw);
    Q = (Q+Q')/2;
    f = -((A')*lSiw);
    %------------------
    options = optimoptions('quadprog','Display','none');
    [x,~,~,~,~] = quadprog( Q, f, Ain, bin, [], [], [], [], [], options );
end
% -------------------------------------------------------------------------
% Unnormalize
x(16:21) = x(16:21)*1000;
md       = (x(16)+x(19)+x(21))/3;
x(1:15)  = x(1:15)/(md*md);
if(opts.estimS0)
    x(22) = exp(x(22)*sc);
else
    x(22) = 1;
end
% -------------------------------------------------------------------------

end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [Ain,ADc] = create_reconst_matrixes_dki(nconst,estimS0)
gic = designGradients( nconst, 'plot', false, 'verbose', false );
% ------------------
% Kurtosis-related:
u   = eye(15);
u   = reshape(u,[1,1,15,15]); % 1 x 1 x 15 x 15
AKc = hot2signal( u, gic );   % 1 x 1 x 15 x G
AKc = permute(AKc,[4,3,1,2]); % G x 15;
% ------------------
% Diffusion tensor-related:
ADc = [ gic(:,1).*gic(:,1), 2*gic(:,1).*gic(:,2), 2*gic(:,1).*gic(:,3), ...
    gic(:,2).*gic(:,2), 2*gic(:,2).*gic(:,3), gic(:,3).*gic(:,3) ];
ADc = 1000*ADc; % The 1000* factor is necessary for inner consistency with the mex
% ------------------
% This is just an estimate of the range for which the Kurtosis model
% holds, see Lazar et al., MRM 60(4): 774-781 (2008):
bmax = 3000;
Ain = [ zeros(nconst,15), -ADc;    % positive diffusion
    -AKc, zeros(nconst,6);    % positive kurtosis
    AKc,        -(3/bmax)*ADc ]; % decaying attenuation signal
% Add the column to account for the log(S0) unknown, if necessary:
if(estimS0)
    Ain = [ Ain, zeros(3*nconst,1) ];
end
end
