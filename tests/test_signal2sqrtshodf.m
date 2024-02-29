% test_signal2posshodf.m
% Change 'use_NR' to use Newton-Raphson's or Levenberg-Marquardt's
% Change 'wise_init' to choose the way the algorithm is initialized

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(13061981);
clear;
close('all');
clc;
wise_init = 1;
use_NR    = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L    = 2;
FOV  = 4;
psnr = 10;
tgt  = 3; % Voxel to plot within the FOV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts.nu = 0.005;
if(use_NR)
    opts.algorithm = 'N';
else
    opts.algorithm = 'L';
end
opts.T = 200;
opts.maxfails = 10;
opts.thCost = 0.000001;
opts.thGrad = 0.00001;
opts.thCond = 0.00001;
opts.rho0 = 1.0;
opts.minrcn = 1.0e-5;
opts.psi0eps = 0.001;
opts.maxthreads = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U1  = [ [1/2;sqrt(3)/2;0], [sqrt(3)/2;-1/2;0], [0;0;1] ];
l1  = diag([1.7e-3,0.3e-3,0.2e-3]);
DT1 = U1*l1*U1';
U2  = [ [sqrt(3)/2;-1/2;0], [1/2;sqrt(3)/2;0], [0;0;1] ];
l2  = diag([1.8e-3,0.1e-3,0.2e-3]);
DT2 = U2*l2*U2';
% ---
% The ODF for a tensor model reads:
%  Phi(r) = (r^T(D^{-1})r)^{-3/2}/(4*pi*sqrt(det(D)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gi   = icosamplesSphere(5);        % M x 3
odf1 = (DT1\(gi'))';               % M x 3
odf2 = (DT2\(gi'))';               % M x 3
odf1 = sum(gi.*odf1,2);            % M x 1
odf2 = sum(gi.*odf2,2);            % M x 1
odf1 = odf1.^(-3/2);               % M x 1
odf2 = odf2.^(-3/2);               % M x 1
odf1 = odf1/(4*pi)/sqrt(det(DT1)); % M x 1
odf2 = odf2/(4*pi)/sqrt(det(DT2)); % M x 1
odf  = (odf1+odf2)/2;
% ------
B  = GenerateSHMatrix(2*L,gi);
LM = GenerateSHEigMatrix(2*L);
LS = (B'*B+opts.nu*(LM'*LM))\(B');
SH = LS*odf;
SH = SH/SH(1)/sqrt(4*pi);          % Kp x 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lpar  = (linspace(2.4e-3,2.0e-3,FOV))'; % FOV x 1
lperp = (linspace(0.1e-3,0.4e-3,FOV))'; % FOV x 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bs  = [1000;2500;5000];
Ns  = length(bs);
gi  = icosamplesSphere(4);
Ni  = size(gi,1);
gi  = repmat(gi,[length(bs),1]);
N   = size(gi,1);
bi  = ones(N,1);
for s=1:Ns
    bi( (s-1)*Ni+1:s*Ni ) = bi( (s-1)*Ni+1:s*Ni ) * bs(s);
end
Ylm = GenerateSHMatrix(2*L,gi); % N x Kp
K   = (L+1)*(L+2)/2;
Kp  = (2*L+1)*(2*L+2)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lmb    = dmri_compute_convolution_weights_ODF(bs,lpar,lperp,2*L); % FOV x (L+1) x Ns
lmb    = permute( lmb, [2,3,1] ); % (L+1) x Ns x FOV
pshell = [zeros(Ni,1);ones(Ni,1);2*ones(Ni,1)];
lambda = ones(L+1,N,FOV);
for s=1:Ns
    lambda(:,(s-1)*Ni+1:s*Ni,:) = lambda(:,(s-1)*Ni+1:s*Ni,:).*lmb(:,s,:); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH2 = ones(Kp,N,FOV);
pos = 1;
for l=0:2:2*L
    li  = l/2+1;
    cnt = 2*l+1;
    SH2(pos:pos+cnt-1,:,:) = reshape(SH(pos:pos+cnt-1),[cnt,1,1]).*lambda(li,:,:);
    pos = pos+cnt;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = sum( SH2.*permute(Ylm,[2,1,3]), 1 ); % 1 x N x FOV
E = permute(E,[2,3,1]); % N x FOV
E = E + randn(size(E))/psnr + 1i*randn(size(E))/psnr;
E = abs(E);
E = min(E,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[bs,ps,Ns] = auto_detect_shells(bi,300);
p = cell(1,Ns);
for n=1:Ns
    p{n} = find( abs(ps-n)<0.5 );
end
SH0 = compute_shodf_from_micro( ...
    E', ...
    lpar, ...
    lperp, ...
    bs, ...
    gi, ...
    p, ...
    2*L, ...
    opts.nu, ...
    false, ...
    100 );
SH0 = [ones(1,FOV)/sqrt(4*pi);SH0']; % Kp x FOV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Es = Ylm*SH0;         % (N x Kp)*(Kp x FOV)->(N x FOV)
Es = sqrt(max(Es,0)); % (N x FOV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(wise_init)
    Bp   = GenerateSHMatrix(L,gi);                   % N x K
    LMp  = GenerateSHEigMatrix(L);                   % K x K
    LSp  = (Bp'*Bp+(opts.nu/1000)*(LMp'*LMp))\(Bp'); % K x N
    SHs  = LSp*Es;                                   % (K x N)*(N x FOV)->(K x FOV)
    psi0 = SHs./sqrt(sum(SHs.*SHs,1));               % K x FOV
else
    psi0 = rand(K,FOV);
    psi0(1,:) = 1;
    psi0 = psi0./sqrt(sum(psi0.*psi0,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(use_NR)
    NT        = 201;
    Q_        = zeros(NT,FOV);
    gradnorm_ = zeros(NT,FOV);
    cnorm_    = zeros(NT,FOV);
    nit_      = zeros(NT,FOV);
    psi_      = zeros(NT,FOV);
    mu_       = zeros(NT,FOV);
    msg = sprintf('Step %i of %i',0,NT);
    ml  = length(msg);
    fprintf(1,msg);
    for T=0:NT-1
        %%%
        fprintf(1,repmat('\b',[1,ml]));
        msg = sprintf('Step %i of %i',T,NT);
        ml  = length(msg);
        fprintf(1,msg);
        %%%
        opts.T = T;
        %%%
        [psi,nit,mu,Q,q0,q1] = ...
            signal2sqrtshodf( E, lmb, pshell, psi0, Ylm, opts );
        psi_(T+1,:)      = psi(1,:);
        nit_(T+1,:)      = nit;
        Q_(T+1,:)        = Q;
        gradnorm_(T+1,:) = q0;
        cnorm_(T+1,:)    = q1;
        mu_(T+1,:)       = mu;
    end
    fprintf(1,repmat('\b',[1,ml]));
    fprintf(1,'\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LW = 2;
    T  = 0:NT-1;
    % ---
    close(figure(2));
    hf2 = figure();
    hold('on');
    grid('on');
    T1 = T( nit_(:,tgt) < 0 );
    P1 = Q_( nit_(:,tgt) < 0, tgt );
    T2 = T( nit_(:,tgt) >= 0 );
    P2 = Q_( nit_(:,tgt) >= 0, tgt );
    plot(T1,P1,'LineWidth',LW,'Color',[.5,.0,.0]);
    plot(T2,P2,'LineWidth',LW,'Color',[.0,.5,.0]);
    xlabel('iteration');
    ylabel('Lagrangian');
    % ---
    close(figure(3));
    hf3 = figure(3);
    hold('on');
    grid('on');
    T1 = T( nit_(:,tgt) < 0 );
    P1 = gradnorm_( nit_(:,tgt) < 0, tgt );
    T2 = T( nit_(:,tgt) >= 0 );
    P2 = gradnorm_( nit_(:,tgt) >= 0, tgt );
    plot(T1,P1,'LineWidth',LW,'Color',[.5,.0,.0]);
    plot(T2,P2,'LineWidth',LW,'Color',[.0,.5,.0]);
    xlabel('iteration');
    ylabel('gradient norm');
    % ---
    close(figure(4));
    hf4 = figure(4);
    hold('on');
    grid('on');
    T1 = T( nit_(:,tgt) < 0 );
    P1 = cnorm_( nit_(:,tgt) < 0, tgt );
    T2 = T( nit_(:,tgt) >= 0 );
    P2 = cnorm_( nit_(:,tgt) >= 0, tgt );
    plot(T1,P1,'LineWidth',LW,'Color',[.5,.0,.0]);
    plot(T2,P2,'LineWidth',LW,'Color',[.0,.5,.0]);
    xlabel('iteration');
    ylabel('Constraint error');
    % ---
    close(figure(5));
    hf5 = figure(5);
    hold('on');
    grid('on');
    T1 = T( nit_(:,tgt) < 0 );
    P1 = mu_( nit_(:,tgt) < 0, tgt );
    T2 = T( nit_(:,tgt) >= 0 );
    P2 = mu_( nit_(:,tgt) >= 0, tgt );
    plot(T1,P1,'LineWidth',LW,'Color',[.5,.0,.0]);
    plot(T2,P2,'LineWidth',LW,'Color',[.0,.5,.0]);
    xlabel('iteration');
    ylabel('Lagrange multiplier');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    NT        = 201;
    Q_        = zeros(NT,FOV);
    gradnorm_ = zeros(NT,FOV);
    nit_      = zeros(NT,FOV);
    psi_      = zeros(NT,FOV);
    msg = sprintf('Step %i of %i',0,NT);
    ml  = length(msg);
    fprintf(1,msg);
    for T=0:NT-1
        %%%
        fprintf(1,repmat('\b',[1,ml]));
        msg = sprintf('Step %i of %i',T,NT);
        ml  = length(msg);
        fprintf(1,msg);
        %%%
        opts.T = T;
        %%%
        [psi,nit,mu,Q,q0,q1] = ...
            signal2sqrtshodf( E, lmb, pshell, psi0, Ylm, opts );
        psi_(T+1,:)      = psi(1,:);
        nit_(T+1,:)      = nit;
        Q_(T+1,:)        = Q;
        gradnorm_(T+1,:) = q0;
    end
    fprintf(1,repmat('\b',[1,ml]));
    fprintf(1,'\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LW = 2;
    T  = 0:NT-1;
    % ---
    close(figure(2));
    hf2 = figure();
    hold('on');
    grid('on');
    T1 = T( nit_(:,tgt) < 0 );
    P1 = Q_( nit_(:,tgt) < 0, tgt );
    T2 = T( nit_(:,tgt) >= 0 );
    P2 = Q_( nit_(:,tgt) >= 0, tgt );
    plot(T1,P1,'LineWidth',LW,'Color',[.5,.0,.0]);
    plot(T2,P2,'LineWidth',LW,'Color',[.0,.5,.0]);
    xlabel('iteration');
    ylabel('cost function');
    % ---
    close(figure(3));
    hf3 = figure(3);
    hold('on');
    grid('on');
    T1 = T( nit_(:,tgt) < 0 );
    P1 = gradnorm_( nit_(:,tgt) < 0, tgt );
    T2 = T( nit_(:,tgt) >= 0 );
    P2 = gradnorm_( nit_(:,tgt) >= 0, tgt );
    plot(T1,P1,'LineWidth',LW,'Color',[.5,.0,.0]);
    plot(T2,P2,'LineWidth',LW,'Color',[.0,.5,.0]);
    xlabel('iteration');
    ylabel('gradient norm');
    % ---
    close(figure(4));
    hf4 = figure(4);
    hold('on');
    grid('on');
    T1 = T( nit_(:,tgt) < 0 );
    P1 = psi_( nit_(:,tgt) < 0, tgt );
    T2 = T( nit_(:,tgt) >= 0 );
    P2 = psi_( nit_(:,tgt) >= 0, tgt );
    plot(T1,P1,'LineWidth',LW,'Color',[.5,.0,.0]);
    plot(T2,P2,'LineWidth',LW,'Color',[.0,.5,.0]);
    xlabel('iteration');
    ylabel('psi[0]');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
