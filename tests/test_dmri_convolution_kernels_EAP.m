% test_dmri_convolution_kernels_EAP.m
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf(1,'Design experiment... ');
ref   = 7;
G     = icosamplesSphere(ref,'O1',true);
Ri    = [0.00375;0.0075;0.015;0.03;0.06];
N     = 20;
lpar  = linspace(1.0e-3,3.0e-3,N);
lperp = linspace(0.2e-3,0.9e-3,N);
[lpar,lperp] = meshgrid(lpar,lperp);
lpar  = lpar(:);
lperp = lperp(:);
tau   = 35e-3;
L     = 8;
lmb   = 0;
fprintf(1,'Done\n');
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf(1,'Compute analytical weights... ');
elODF  = dmri_compute_convolution_weights_EAP(Ri,lpar,lperp,tau,L);
fprintf(1,'Done\n');
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cuv2   = G(:,3).*G(:,3);   % Gx1
cuv2   = cuv2';            % 1xG
delta  = 1./lpar-1./lperp; % N^2x1
B      = GenerateSHMatrix( L, G );  % GxK
T      = GenerateSHEigMatrix( L );  % KxK
SH     = (((B')*B + lmb*T*T)\eye(size(T,1)))*(B'); % KxG
ll     = 0:2:L;
id     = (ll+1).*(ll+2)/2;
id     = [0,id(1:end-1)] + ll + 1;
SH     = SH(id,:);             % (L/2+1)xG
crr    = sqrt(4*pi./(2*ll+1)); % 1x(L/2+1)
SH     = diag(crr)*SH;         % (L/2+1)xG       
SH     = SH';                  % Gx(L/2+1)
Km1    = (4*pi*tau)*(4*pi*tau)*(4*pi*tau);
Km1    = sqrt(Km1);
cols   = [ ...
    .5, .0, .0;
    .0, .5, .0;
    .0, .0, .5;
    .6, .6, .0;
    .0, .6, .6;
    .6, .0, .6
    ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf(1,'Numerically compute ODF EAP weights... ');
for nb=1:length(Ri)
    numis = zeros(N*N,L/2+1);
    ptr = 1:N;
    for n=1:N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        arg = delta(ptr,:)*cuv2; % NxG 
        arg = arg/(4*tau);       % NxG
        arg = Ri(nb)*Ri(nb)*arg; % NxG
        arg = exp(-arg);         % NxG
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        arg = arg*SH;            % Nx(L/2+1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        arg = arg./Km1;          % NxG
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numis(ptr,:) = arg;      % Nx(L/2+1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ptr = ptr+N;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    hf = figure(1000+nb);
    close(hf);
    figure(1000+nb);
    grid('on');
    hold('on');
    % ---
    nrm = exp( -Ri(nb)*Ri(nb) ./ (4*tau*lperp) ); % N^2x1
    elODF(:,:,nb) = elODF(:,:,nb)./nrm;           % N^2x(L+1)(L+2)/2xNb
    nrm = sqrt(lpar.*lperp.*lperp);
    elODF(:,:,nb) = elODF(:,:,nb).*nrm;
    % ---
    for nl=1:L/2+1
        plot(sqrt(-delta/4/tau),numis(:,nl),'o','Color',cols(nl,:));
        plot(sqrt(-delta/4/tau),elODF(:,nl,nb),'*','Color',cols(nl,:));
    end
    xlabel('(|\lambda_{||}^{-1}-\lambda_{\perp}^{-1}|/(4\tau))^{1/2} [mm^{-1}]','FontSize',20);
    ylabel('{(\lambda_{||}\lambda_{\perp}^2})^{1/2}e_{l,EAP}^{i}/exp({-R_i^2}/{4\tau\lambda_{\perp}})','FontSize',20);
    title(['EAP convolution coefficients at R_i = ',num2str(Ri(nb)),'   ---   *=analytical; o=numerical'],'FontSize',20);
end
fprintf(1,'Done\n');
toc;
