function test_micro2eap

sf = check_software_platform;
if(sf==2)
    pkg load statistics;
end

%%% test_micro2eap.m
clc;
close('all');
L   = 8;
tau = 40e-3;
[atti,gi,bi] = generate_atti_signal;
atti = reshape( atti, [1,1,1,size(atti,1)] );

tic;
[lpar,lperp,~] = atti2micro( atti, gi, bi, 'tl', 1.0e-6, 'tu', 1-1.0e-6, ...
    'ADC0', 3.0e-3, 'usef', false, 'bth', 100, 'mlperp',0.0, ...
    'mu', 0.0, 'verbose', true );
fprintf(1,'It took %1.5f seconds to fit the u-structure model\n',toc);

tic;
sh = micro2shodf( atti, gi, bi, lpar, lperp, [], 'L', L, ...
    'tl', 1.0e-6, 'tu', 1-1.0e-6, ...
    'optimal', true, 'lambda', 0.0005 );
fprintf(1,'It took %1.5f seconds to fit the ODFs\n',toc);

% -------------------------------------------------------------------------
% Define a proper vector for q:
NHT  = 2*8192+1; % Number of points for Hankel transforms
bmin = 0;
bmax = 100000;
qmin = sqrt(bmin/(4*pi*pi*tau)); % mm^{-1}
qmax = sqrt(bmax/(4*pi*pi*tau)); % mm^{-1}
q    = linspace(qmin,qmax,NHT);
b    = 4*pi*pi*tau*(q.*q);
% Compute the profiles of the convolution weights:
e    = dmri_compute_convolution_weights_ODF( ...
    b', lpar(1), lperp(1), L ); % 1 x (L/2+1) x NHT
% -------------------------------------------------------------------------
close(figure(1));
figure(1);
for l=0:2:L
    subplot(1,L/2+1,l/2+1);
    plot(q,squeeze(e(1,l/2+1,:))','LineWidth',2,'Color',[0,0,0]);
    grid('on');
    xlabel('q [mm^{-1}]','FontSize',14);
    ylabel(['\Lambda_{',num2str(l),'}'],'FontSize',14);
end
% -------------------------------------------------------------------------
% Define the vector of values of R where the EAP will be computed:
NVR  = 200;
Rmin = 0;
Rmax = 0.05;
R0   = linspace(Rmin,Rmax,NVR);
Pr0  = zeros(L/2+1,NVR);
% -------------------------------------------------------------------------
% For each l-value, i.e. for each of the 2nd dimension, compute the 
% modified Hankel transform, i.e.:
%
% 4*pi*(-1)^(l/2)*int( q^2*e_l(q)*besselj(l,2*pi*q*R), q=0..infty );
arg  = 2*pi*(R0')*q;  % NVR x NHT
for l=0:2:L
    y = spbesselj(l,arg);            % NVR x NHT
    y = (q.*q).*y;                   % NVR x NHT, broadcasted
    y = y.*(squeeze(e(1,l/2+1,:))'); % NVR x NHT, broadcasted
    % Simpsons rule:
    ck = [1,4,1];
    H  = conv2(y,ck,'same');         % NVR x NHT
    H  = H(:,2:2:end-1);             % NVR x (NHT-1)/2
    dq = q(1:2:end);                 % 1 x (NHT+1)/2
    dq = diff(dq)/6;                 % 1 x (NHT-1)/2
    H  = H.*dq;                      % NVR x (NHT-1)/2, broadcasted
    H  = sum(H,2);                   % NVR x 1
    % Normalization:
    H  = 4*pi*H*realpow(-1,l/2);     % NVR x 1
    % Assignment:
    Pr0(l/2+1,:) = H';
end
% -------------------------------------------------------------------------
close(figure(2));
figure(2);
for l=0:2:L
    subplot(1,L/2+1,l/2+1);
    plot(R0,Pr0(l/2+1,:),'LineWidth',2,'Color',[0,0,0]);
    grid('on');
    xlabel('R [mm]','FontSize',14);
    ylabel(['\Gamma_{',num2str(l),'}'],'FontSize',14);
end
% -------------------------------------------------------------------------


% Design a proper gradient sampling:
NG = 128;
gi = designGradients(NG);
% Design a SH analysis matrix:
B   = GenerateSHMatrix( L, gi );    % NGxK, where K=(L+1)(L+2)/2
WLS = (B'*B)\(B');                  % (KxK)^(-1) * (KxGN) -> KxNG
% The results will be stored in:
SH1 = zeros( (L+1)*(L+2)/2, length(R0) );
SH2 = zeros( (L+1)*(L+2)/2, length(R0) );
% For each value of R...
ptr = dmri_sh_expand_coeffs(L);
for n=1:length(R0)
    % ----------------------------
    % Call the ad hoc function to compute the EAP:
    eap = micro2eap( sh, lpar, lperp, [], ...
        gi, ...
        ones(NG,1)*R0(n), ...
        'tau', tau, 'rith', 10 ); % 1x1x1xNG
    % Fit this signal in the basis of SH:
    SH1(:,n) = WLS*(eap(:));
    % ----------------------------
    % Multiply the SH coeffs times the Hankel transform at this R:
    SH2(:,n) = sh(:).*Pr0(ptr,n);
    % ----------------------------
end

close(figure(3));
figure(3);
p = 1;
for l=0:2:L
    for m=-l:l
        plot( R0, SH1(p,:), 'ro', R0, SH2(p,:), 'gx' );
        grid('on');
        title(sprintf('l=%d, m=%d',l,m));
        drawnow;
        pause;
        p =p+1;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [atti,gi,bi] = generate_atti_signal
[gi,bi] = generate_ms_sampling;
Diso = 3.0e-3;
f1 = 0.6;
% -------------------------------------
% Fiber 1
dir   = [1.1,0.2,0.1];
dir   = dir/norm(dir);
lpar  = 2.0e-3;
lperp = 0.3e-3;
ffw   = 0.1;
fin   = 0.4;
uTv   = sum(gi.*dir,2);
uTv   = uTv.*uTv;
% ---
cfw   = exp( -bi*Diso );
cin   = exp( -bi.*(lpar*uTv) );
cen   = exp( -bi.*((lpar-lperp)*uTv+lperp) );
% ---
S1    = ffw*cfw + fin*cin + (1-ffw-fin)*cen;
% -------------------------------------
% Fiber 2
dir   = [0.04,0.06,0.9];
dir   = dir/norm(dir);
lpar  = 2.1e-3;
lperp = 0.15e-3;
ffw   = 0.075;
fin   = 0.6;
uTv   = sum(gi.*dir,2);
uTv   = uTv.*uTv;
% ---
cfw   = exp( -bi*Diso );
cin   = exp( -bi.*(lpar*uTv) );
cen   = exp( -bi.*((lpar-lperp)*uTv+lperp) );
% ---
S2    = ffw*cfw + fin*cin + (1-ffw-fin)*cen;
% -------------------------------------
atti = f1*S1 + (1-f1)*S2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gi,bi] = generate_ms_sampling
bs = [500,1000,1500,2000,3000,5000];
Ns = [15,30,30,40,50,60];
bi = [];
gi = [];
for k=1:length(bs)
    bi = [ bi; bs(k)*ones(Ns(k),1) ];    %#ok<AGROW>
    gi = [ gi; designGradients(Ns(k)) ]; %#ok<AGROW>
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = spbesselj(nu,x)
% Compute the spherical Bessel function from the regular
% Bessel function:
y = sqrt(pi./(2*x)).*besselj(nu+0.5,x);
if(nu==0)
    y(x==0) = 1;
else
    y(x==0) = 0;
end
end

