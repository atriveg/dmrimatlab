function test_atti2freewater(opt)

sf = check_software_platform;
if(sf==2)
    pkg load statistics;
end

rng(15,'twister');

PLOTOTHERSCALARS = true; % C_l, C_p, C_s, FA...
SAMPLING  = 12;
if(SAMPLING==5)
    bextra = input('bextra: ');
end
EIGENVALS = 'gaussian';    % 'uniform'/'gaussian'
NOISE     = 'rice';        % 'rice'/'gauss'
METHOD    = 'micro-fixed'; % 'tensor'/'positive'/'micro'/'micro-fixed'
nu_       = 0.1;
lpar_     = 2.0e-3;
Nnb       = 4;
PSNR      = 18;
if(nargin>0)
    if(isfield(opt,'PLOTOTHERSCALARS'))
        PLOTOTHERSCALARS = opt.PLOTOTHERSCALARS;
    end
    if(isfield(opt,'SAMPLING'))
        SAMPLING = opt.SAMPLING;
    end
    if(isfield(opt,'EIGENVALS'))
        EIGENVALS = opt.EIGENVALS;
    end
    if(isfield(opt,'NOISE'))
        NOISE = opt.NOISE;
    end
    if(isfield(opt,'METHOD'))
        METHOD = opt.METHOD;
    end
    if(isfield(opt,'nu'))
        nu_ = opt.nu;
    end
    if(isfield(opt,'lpar'))
        lpar_ = opt.lpar;
    end
    if(isfield(opt,'gi'))
        gi = opt.gi;
        G  = size(gi,1);
    end
    if(isfield(opt,'bi'))
        bi = opt.bi;
    end
    if(isfield(opt,'Nnb'))
        Nnb = opt.Nnb;
    end
    if(isfield(opt,'PSNR'))
        PSNR  = opt.PSNR;
        sigma = 1/PSNR;
    end
end

switch(SAMPLING)
    case 0
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SYNTHETIC SAMPLING SCHEME:
        [gi,linit,lend] = icosamplesSphere(4);
        G  = size(gi,1);
        bi = zeros(G,1);
        bi(linit(1):lend(1),:) = 2000; %  6 grads
        bi(linit(2):lend(2),:) = 1667; % 10 grads
        bi(linit(3):lend(3),:) = 1333; % 15 grads
        bi(linit(4):lend(4),:) = 1000; % 30 grads
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        %%%%%%%%%%%
        PSNR  = 30;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAMPLING SCHEME TAKEN FORM HCP-MGH 1007
        load('test_atti2freewater.mat');
        idg = (bvals_reduced>1);
        bi  = bvals_reduced(idg,:);
        gi  = Gi_reduced(idg,:);
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        PSNR  = 20;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAMPLING SCHEME TAKEN FORM PDD case CO_p07090
        load('test_atti2freewater.mat');
        bi  = bi_PDD;
        gi  = gi_PDD;
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        % PSNR goes between 8 (10% quantile) and 22 (90% quantile) with a
        % typical value (mode) 12
        PSNR  = 20;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SYNTHETIC SAMPLING SCHEME:
        [gi,linit,lend] = icosamplesSphere(4);
        G  = size(gi,1);
        bi = zeros(G,1);
        bi(linit(1):lend(1),:) = 1000; %  6 grads
        bi(linit(2):lend(2),:) = 2000; % 10 grads
        bi(linit(3):lend(3),:) = 2000; % 15 grads
        bi(linit(4):lend(4),:) = 3000; % 30 grads
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        %%%%%%%%%%%
        PSNR  = 16;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 5
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ANOTHER SYNTHETIC SAMPLING SCHEME:
        gi1 = icosamplesSphere(1); %  6 grads
        gi2 = icosamplesSphere(4); % 61 grads
        gi  = [gi1;gi2];
        G   = size(gi,1);
        bi  = zeros(G,1);
        bi(1:6,1)   = bextra;
        bi(7:end,1) = 1000;
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        %%%%%%%%%%%
        PSNR  = 16;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 6
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % YET ANOTHER SYNTHETIC SAMPLING SCHEME:
        gi1 = icosamplesSphere(1); %  6 grads
        gi2 = icosamplesSphere(4); % 61 grads
        gi  = [ ...
            gi1;
            gi1;
            gi1;
            gi1;
            gi1;
            gi1;
            gi1;
            gi1;
            gi2 ];
        G   = size(gi,1);
        bi  = zeros(G,1);
        bi(1:6,1)    = 250;
        bi(7:12,1)   = 300;
        bi(13:18,1)  = 350;
        bi(19:24,1)  = 400;
        bi(25:30,1)  = 450;
        bi(31:36,1)  = 500;
        bi(37:42,1)  = 550;
        bi(43:48,1)  = 600;
        bi(49:end,1) = 1000;
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        %%%%%%%%%%%
        PSNR  = 16;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 7
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % YET ANOTHER SYNTHETIC SAMPLING SCHEME:
        gi1 = icosamplesSphere(1); %  6 grads
        gi2 = icosamplesSphere(4); % 61 grads
        gi  = [ ...
            gi1;
            gi2;
            gi1 ];
        G   = size(gi,1);
        bi  = zeros(G,1);
        bi(1:6,1)     = 500;
        bi(7:67,1)    = 1000;
        bi(68:73,1)   = 1500;
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        %%%%%%%%%%%
        PSNR  = 20;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 8
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAMPLING SCHEME TAKEN FROM HOY'S PAPER:
        gi  = designGradients(32,'verbose',true);
        G   = size(gi,1);
        bi  = [ 500*ones(G,1); 1500*ones(G,1) ];
        gi  = [gi;gi];
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        PSNR  = 40;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 9
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAMPLING SCHEME TAKEN FROM HOY'S PAPER PLUS CENTRAL SHELL:
        gi1 = designGradients(32,'verbose',true); % 32 gradients
        gi2 = gi1;                                % 32 gradients
        G1  = size(gi1,1);
        G2  = size(gi2,1);
        bi  = zeros(2*G1+G2,1);
        bi(1:G1)           = 500;
        bi(G1+1:2*G1)      = 1500;
        bi(2*G1+1:2*G1+G2) = 1000;
        gi  = [gi1;gi1;gi2];
        bpt = (abs(bi)>100);
        Nnb = 8;
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        PSNR  = 30;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 10
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAMPLING SCHEME TAKEN FROM DRYAD
        gi1 = designGradients(32,'verbose',true); % 32 gradients
        gi  = [gi1;gi1;gi1;gi1;gi1;gi1;gi1;gi1];
        G1  = size(gi1,1);
        bi  = zeros(6*G1,1);
        bi(1:G1)           = 200;
        bi(1*G1+1:2*G1)    = 400;
        bi(2*G1+1:3*G1)    = 600;
        bi(3*G1+1:4*G1)    = 800;
        bi(4*G1+1:5*G1)    = 1000;
        bi(5*G1+1:6*G1)    = 1200;
        bi(6*G1+1:7*G1)    = 1400;
        bi(7*G1+1:8*G1)    = 1600;
        bpt = (abs(bi)>100);
        Nnb = 8;
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        PSNR  = 40;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 11
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IMPROVED SAMPLING SCHEME TAKEN FROM HOY'S PAPER PLUS CENTRA SHELL:
        gi  = designGradients(64,'verbose',true);
        G   = size(gi,1);
        bi  = [ 500*ones(G,1); 1500*ones(G,1); 1000*ones(G,1) ];
        gi  = [gi;gi;gi];
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        PSNR  = 20;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 12
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % IMPROVED SAMPLING SCHEME TAKEN FROM HOY'S PAPER PLUS CENTRA SHELL:
        gi1 = designGradients(30,'verbose',true);
        gi2 = designGradients(20,'verbose',true);
        G1  = size(gi1,1);
        G2  = size(gi2,1);
        bi  = [ 1200*ones(G1,1); 500*ones(G2,1) ];
        gi  = [gi1;gi2];
        bpt = (abs(bi)>100);
        Nnb = 4;
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        PSNR  = 20;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

N = 200; % Number of synthetic voxels to be tested
DISO = 3.0e-3; % isotropic diffusion coefficient of free water
L = 4; % Order of the spherical harmonics
lambda = 0.001; % Tikhonov regularisation for Spherical Harmonics
tau = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precompute random rotations for each compartment
q1  = randn(N,4);
nq1 = sqrt(sum(q1.*q1,2));
q1  = q1./nq1;
%%%
ph2 = rand(N,1)*(pi/2);
q2  = [ cos(ph2/2), zeros(N,1), zeros(N,1), sin(ph2/2) ];
%%%
ph3 = (rand(N,1)-1/2)*pi/6;
ax3 = randn(N,3);
nax = sqrt(sum(ax3.*ax3,2));
ax3 = ax3./nax;
q3  = [ cos(ph3/2), ax3.*sin(ph3/2) ];
%%%
q   = { q1, multiply_quatrot(q1,q2), multiply_quatrot(q1,q3) };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f   = (1:1:20)/20;
hp1 = zeros(1,3);
for C=1:3
    f2cum  = zeros(N,length(f));
    SHcum  = zeros(N,(L+1)*(L+2)/2,length(f));
    DTcum  = zeros(N,6,length(f));
    labels = cell(1,length(f));
    txtl   = 0;
    for nf=1:length(f)
        % -- User feedback
        if(rem(nf,2)==1)
            labels{nf} = num2str(f(nf));
        else
            labels{nf} = '';
        end
        fprintf(1,repmat('\b',[1,txtl]));
        msg = sprintf('f=%1.2f',f(nf));
        fprintf(1,msg);
        txtl = length(msg);
        % -- Initialize the signal and the volume fractions of each
        %    compartment:
        Si = zeros(G,N); % GxN
        fi = ones(C,N)/C + (0.2/C)*(rand(C,N)-1/2); % CxN 
        fi = bsxfun( @(x,y)(x./y), fi, sum(fi,1) ); % CxN
        % -- For each compartment:
        for nc=1:C
            % - Compute a random rotation for the eigenvectors of the
            %   compartment:
            ii = [1,0,0];
            jj = [0,1,0];
            kk = [0,0,1];
            switch(nc)
                case 1
                    xx = apply_quatrot( q{nc}, ii )';
                    yy = apply_quatrot( q{nc}, jj )';
                    zz = apply_quatrot( q{nc}, kk )';
                case 2
                    xx = apply_quatrot( q{nc}, jj )';
                    yy = apply_quatrot( q{nc}, kk )';
                    zz = apply_quatrot( q{nc}, ii )';
                case 3
                    xx = apply_quatrot( q{nc}, kk )';
                    yy = apply_quatrot( q{nc}, ii )';
                    zz = apply_quatrot( q{nc}, jj )';
            end
            % - Generate random eigenvalues for this compartment:
            switch(EIGENVALS)
                case 'uniform'
                    l1 = 1.55 + 0.2*rand(1,N); % 1xN
                    l2 = 0.25 + 0.2*rand(1,N); % 1xN
                    l3 = 0.25 + 0.2*rand(1,N); % 1xN
                case 'gaussian'
                    l1 = 1.34 + 0.27*randn(1,N);
                    l2 = 0.39 + 0.10*randn(1,N);
                    l3 = 0.23 + 0.08*randn(1,N);
                    l1 = min( max(l1,0), DISO*1000 );
                    l2 = min( max(l2,0), DISO*1000 );
                    l3 = min( max(l3,0), DISO*1000 );
                    l2 = min(l1,l2);
                    l3 = min(l2,l3);
                otherwise
                    error(['Unknown eigenvalues distribution: ',EIGENVALS]);
            end
            % - Compute the 6 independent components of the diffusion tensor
            %   corresponding to this compartment:
            DT = zeros(6,N); % 6xN
            DT(1,:) = l1.*(xx(1,:).*xx(1,:)) + l2.*(yy(1,:).*yy(1,:)) + l3.*(zz(1,:).*zz(1,:));
            DT(2,:) = l1.*(xx(1,:).*xx(2,:)) + l2.*(yy(1,:).*yy(2,:)) + l3.*(zz(1,:).*zz(2,:));
            DT(3,:) = l1.*(xx(1,:).*xx(3,:)) + l2.*(yy(1,:).*yy(3,:)) + l3.*(zz(1,:).*zz(3,:));
            DT(4,:) = l1.*(xx(2,:).*xx(2,:)) + l2.*(yy(2,:).*yy(2,:)) + l3.*(zz(2,:).*zz(2,:));
            DT(5,:) = l1.*(xx(2,:).*xx(3,:)) + l2.*(yy(2,:).*yy(3,:)) + l3.*(zz(2,:).*zz(3,:));
            DT(6,:) = l1.*(xx(3,:).*xx(3,:)) + l2.*(yy(3,:).*yy(3,:)) + l3.*(zz(3,:).*zz(3,:));
            DT = 1.0e-3*DT; % 6xN
            % - Compute the quadratic form corresponing to the present
            %   compartment for all gradient directions:
            QF = ...
                (gi(:,1).*gi(:,1).*bi(:))*DT(1,:) + ...
                2*(gi(:,1).*gi(:,2).*bi(:))*DT(2,:) + ...
                2*(gi(:,1).*gi(:,3).*bi(:))*DT(3,:) + ...
                (gi(:,2).*gi(:,2).*bi(:))*DT(4,:) + ...
                2*(gi(:,2).*gi(:,3).*bi(:))*DT(5,:) + ...
                (gi(:,3).*gi(:,3).*bi(:))*DT(6,:); % GxN
            
            % - Compute the corresponding signal for this compartment:
            Si = Si + exp(-QF).*repmat(fi(nc,:),[G,1]); % GxN
        end
        % -- Divide by the total volume fraction:
        Si = Si./repmat(sum(fi,1),[G,1]); % GxN
        
        % -- Compute the overall signal accounting for the free water:
        Si2 = repmat((1-f(nf))*exp(-DISO*bi),[1,N]) + f(nf)*Si; % GxN
        
        % -- Add Rician/Gaussian noise to both the DWI and the baseline
        switch(NOISE)
            case 'rice'
                ni  = sigma*(randn(G,N)+1i*randn(G,N));           % GxN
                n0  = sigma*(randn(1,N)+1i*randn(1,N))/sqrt(Nnb); % GxN
                Si2 = abs(Si2+ni)./repmat(abs(1+n0),[G,1]);       % GxN
            case 'gauss'
                Si2 = Si2 + sigma*(randn(G,N)+randn(G,N));
            otherwise
                error('NOISE must be either ''gauss'' or ''rice''')
        end
        
        % -- Run the algorithm:
        switch(METHOD)
            case 'tensor'
                [f2,DT2] = atti2freewaterTensor( reshape(Si2',[N,1,1,G]), gi, bi, ...
                    'verbose', false, ...
                    'O1', true, ...
                    'O2', false, 'fiters', 100, 'fth', 0.0001, ...
                    'O3', true, 'fnriters', 100, 'fnrth', 0.000001, 'lmb', 1, ...
                    'ADC0', DISO, 'pos', false );
                DTcum(:,:,nf) = permute(DT2,[1,4,2,3]);
            case 'positive'
                [f2,SH2] = atti2freewater( reshape(Si2',[N,1,1,G]), gi, bi, ...
                    'L', L, 'lambda', lambda, 't0', 0.00001, ...
                    'verbose', false, ...
                    'O1', true, ...
                    'O2', true, 'tau', tau, 'fnriters', 100, 'fnrth', 0.000001, 'lmb', 1, ...
                    'ADC0', DISO, 'pos', false );
                SHcum(:,:,nf) = permute(SH2,[1,4,2,3]);
            case 'micro'
                [~,~,f2] = atti2micro( reshape(Si2',[N,1,1,G]), gi, bi, ...
                    'lambda', lambda, ...
                    'mu', 0.005, ...
                    'nu', 0.1, ...
                    'bth', 10, ...
                    'forcelpar', false, ...
                    'lpar', 2.0e-3, ...
                    'verbose', false, ...
                    'usef', true, ...
                    'regf', true, ...
                    'fmin', 0.001, 'fmax', 1, ...
                    'nmax', 100, 'dl', 1.0e-7, 'dC', 1.0e-6, ...
                    'nolpwarn', true, ...
                    'ADC0', DISO ...
                    );
            case 'micro-fixed'
                [~,~,f2] = atti2micro( reshape(Si2',[N,1,1,G]), gi, bi, ...
                    'lambda', lambda, ...
                    'mu', 0.00, ...
                    'nu', nu_, ...
                    'bth', 10, ...
                    'forcelpar', true, ...
                    'lpar', lpar_, ...
                    'verbose', false, ...
                    'usef', true, ...
                    'regf', true, ...
                    'fmin', 0.001, 'fmax', 1, ...
                    'nmax', 100, 'dl', 1.0e-7, 'dC', 1.0e-6, ...
                    'nolpwarn', true, ...
                    'ADC0', DISO ...
                    );
        end
        % -- Cast the output in place:
        f2cum(:,nf) = f2(:);
    end
    fprintf(1,repmat('\b',[1,txtl]));
    fprintf(1,'\n');
    
    FS1 = 14;
    FS2 = 15;
    FS3 = 14;
    FS4 = 12;
    CLS = eye(3)*0.5;
    % ---------------------------------------------------------------------
    hf = figure(1001);
    if(C==1)
        close(hf);
        hf = figure(1001);
        plot(f,f,'Color',[0,0,0],'LineWidth',2,'LineStyle','--');
    end
    xlabel('Actual value of fraction-f','FontSize',FS1);
    ylabel('Estimated value of fraction-f','FontSize',FS1);
    switch(METHOD)
        case 'tensor'
            titletxt = ['PSNR=',num2str(PSNR),'; Tensor model'];    
        case 'positive'
            titletxt = ['PSNR=',num2str(PSNR),'; Positive L=',num2str(L),'; \lambda=',num2str(lambda)];
        case 'micro'
            titletxt = ['PSNR=',num2str(PSNR),'; Micro L=',num2str(L),'; \lambda=',num2str(lambda)];
        case 'micro-fixed'
            titletxt = ['PSNR=',num2str(PSNR),'; Micro L=',num2str(L),'; \lambda=',num2str(lambda)];
    end
    title(titletxt,'Interpreter','tex','FontSize',FS2);
    hold('on');
    grid('on');
    hb1 = myboxplot(f2cum,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
    hp1(C) = plot(f,median(f2cum),'Color',CLS(C,:),'LineStyle','--');
    if(C==3)
        hl = legend(hp1,'1-Compartment','2-Compartments','3-Compartments');
        set(hl,'FontSize',FS1);
        set(hl,'Location','SouthEast');
        set(get(hf,'CurrentAxes'),'FontSize',FS4);
        axis([0,1,-0.05,1.05]);
    end
    % ---------------------------------------------------------------------
    if(PLOTOTHERSCALARS)
        switch(METHOD)
            case 'tensor'
                % -------------------------------------------------------------
                % DTcum: Nx6xnf
                DTcum = permute(DTcum,[2,4,1,3]); % 6x1xNxnf
                DT    = zeros(3,3,N,length(f));
                DT(1,1,:,:) = DTcum(1,1,:,:);
                DT(1,2,:,:) = DTcum(2,1,:,:);
                DT(2,1,:,:) = DTcum(2,1,:,:);
                DT(1,3,:,:) = DTcum(3,1,:,:);
                DT(3,1,:,:) = DTcum(3,1,:,:);
                DT(2,2,:,:) = DTcum(4,1,:,:);
                DT(2,3,:,:) = DTcum(5,1,:,:);
                DT(3,2,:,:) = DTcum(5,1,:,:);
                DT(3,3,:,:) = DTcum(6,1,:,:);
                % -------------------------------------------------------------
                FA = zeros(N,length(f));
                ml = zeros(N,length(f));
                cl = zeros(N,length(f));
                cp = zeros(N,length(f));
                cs = zeros(N,length(f));
                for n=1:N
                    for nf=1:length(f)
                        e  = sort(max(eig(DT(:,:,n,nf)),0));
                        e1 = e(3);
                        e2 = e(2);
                        e3 = e(1);
                        % ------
                        FA(n,nf) = sqrt(1/2) * ...
                            sqrt( (e(1)-e(2))*(e(1)-e(2)) + ...
                            (e(1)-e(3))*(e(1)-e(3)) + ...
                            (e(3)-e(2))*(e(3)-e(2)) ) / ...
                            sqrt(e(1)*e(1)+e(2)*e(2)+e(3)*e(3));
                        % ------
                        ml(n,nf) = e3;
                        % ------
                        cl(n,nf) = (e1-e2)/e1;%sqrt(e1*e1+e2*e2+e3*e3);
                        % ------
                        cp(n,nf) = (e2-e3)/e1;%/sqrt(e1*e1+e2*e2+e3*e3);
                        % ------
                        cs(n,nf) = e3/e1;%/sqrt(e1*e1+e2*e2+e3*e3);
                        % ------
                    end
                end
                % -------------------------------------------------------------
                hf2 = figure(2002);
                hf3 = figure(3003);
                hf4 = figure(4004);
                hf5 = figure(5005);
                hf6 = figure(6006);
                if(C==1)
                    close(hf2);
                    close(hf3);
                    close(hf4);
                    close(hf5);
                    close(hf6);
                    hf2 = figure(2002);
                    hf3 = figure(3003);
                    hf4 = figure(4004);
                    hf5 = figure(5005);
                    hf6 = figure(6006);
                end
                % -------------------------------------------------------------
                figure(hf2);
                xlabel('Actual value of fraction-f','FontSize',FS1);
                ylabel('Fractional anisotropy','FontSize',FS1);
                title(['PSNR=',num2str(PSNR),'; FA with tensor model'],'FontSize',FS2);
                hold('on');
                grid('on');
                hb2 = myboxplot(FA,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
                hp2(C) = plot(f,mean(FA),'Color',CLS(C,:),'LineStyle','--'); %#ok<SAGROW>
                % -------------------------------------------------------------
                figure(hf3);
                xlabel('Actual value of fraction-f','FontSize',FS1);
                ylabel('Minimum eigenvalue','FontSize',FS1);
                title(['PSNR=',num2str(PSNR),'; minimum eigenvalue with tensor model'],'FontSize',FS2);
                hold('on');
                grid('on');
                hb3 = myboxplot(ml,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
                hp3(C) = plot(f,mean(ml),'Color',CLS(C,:),'LineStyle','--'); %#ok<SAGROW>
                % -------------------------------------------------------------
                figure(hf4);
                xlabel('Actual value of fraction-f','FontSize',FS1);
                ylabel('Linear coefficient','FontSize',FS1);
                title(['PSNR=',num2str(PSNR),'; C_l with tensor model'],'FontSize',FS2);
                hold('on');
                grid('on');
                hb4 = myboxplot(cl,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
                hp4(C) = plot(f,mean(cl),'Color',CLS(C,:),'LineStyle','--'); %#ok<SAGROW>
                % -------------------------------------------------------------
                figure(hf5);
                xlabel('Actual value of fraction-f','FontSize',FS1);
                ylabel('Planar coefficient','FontSize',FS1);
                title(['PSNR=',num2str(PSNR),'; C_p with tensor model'],'FontSize',FS2);
                hold('on');
                grid('on');
                hb5 = myboxplot(cp,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
                hp5(C) = plot(f,mean(cp),'Color',CLS(C,:),'LineStyle','--'); %#ok<SAGROW>
                % -------------------------------------------------------------
                figure(hf6);
                xlabel('Actual value of fraction-f','FontSize',FS1);
                ylabel('Spherical coefficient','FontSize',FS1);
                title(['PSNR=',num2str(PSNR),'; C_s with tensor model'],'FontSize',FS2);
                hold('on');
                grid('on');
                hb6 = myboxplot(cs,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
                hp6(C) = plot(f,mean(cs),'Color',CLS(C,:),'LineStyle','--'); %#ok<SAGROW>
                % -------------------------------------------------------------
                if(C==3)
                    % ---------------------------------------------------------
                    figure(hf2);
                    hl2 = legend(hp2,'1-Compartment','2-Compartments','3-Compartments');
                    set(hl2,'FontSize',FS1);
                    set(hl2,'Location','SouthEast');
                    set(get(hf2,'CurrentAxes'),'FontSize',FS4);
                    % ---------------------------------------------------------
                    figure(hf3);
                    hl3 = legend(hp3,'1-Compartment','2-Compartments','3-Compartments');
                    set(hl3,'FontSize',FS1);
                    set(hl3,'Location','SouthEast');
                    set(get(hf3,'CurrentAxes'),'FontSize',FS4);
                    % ---------------------------------------------------------
                    figure(hf4);
                    hl4 = legend(hp4,'1-Compartment','2-Compartments','3-Compartments');
                    set(hl4,'FontSize',FS1);
                    set(hl4,'Location','SouthEast');
                    set(get(hf4,'CurrentAxes'),'FontSize',FS4);
                    % ---------------------------------------------------------
                    figure(hf5);
                    hl5 = legend(hp5,'1-Compartment','2-Compartments','3-Compartments');
                    set(hl5,'FontSize',FS1);
                    set(hl5,'Location','SouthEast');
                    set(get(hf5,'CurrentAxes'),'FontSize',FS4);
                    % ---------------------------------------------------------
                    figure(hf6);
                    hl6 = legend(hp6,'1-Compartment','2-Compartments','3-Compartments');
                    set(hl6,'FontSize',FS1);
                    set(hl6,'Location','SouthEast');
                    set(get(hf6,'CurrentAxes'),'FontSize',FS4);
                    % ---------------------------------------------------------
                end
                % -------------------------------------------------------------
            case 'positive'
                % -------------------------------------------------------------
                % SHcum: Nx(L+1)*(L+2)/2xnf
                SHcum = permute(SHcum,[2,1,3]); % (L+1)*(L+2)/2 x N x nf
                GPass = icosamplesSphere(5,'O1',true,'verbose',false);
                B1    = GenerateSHMatrix(L,GPass);
                B2    = GenerateSHMatrix(L+4,GPass);
                Msq   = (B2'*B2)\(B2');
                % -------------------------------------------------------------
                FA = zeros(N,length(f));
                ml = zeros(N,length(f));
                for n=1:N
                    for nf=1:length(f)
                        coeff_sqrt = SHcum(:,n,nf);
                        sig_sqrt   = B1*coeff_sqrt;
                        sig_def    = sig_sqrt.*sig_sqrt;
                        coeff_def  = Msq*sig_def;
                        ml(n,nf) = min(sig_def);
                        FA(n,nf) = ...
                            sum(coeff_def(2:end).*coeff_def(2:end)) / ...
                            sum(coeff_def.*coeff_def);
                    end
                end
                % -------------------------------------------------------------
                hf2 = figure(2002);
                hf3 = figure(3003);
                if(C==1)
                    close(hf2);
                    close(hf3);
                    hf2 = figure(2002);
                    hf3 = figure(3003);
                end
                % -------------------------------------------------------------
                figure(hf2);
                xlabel('Actual value of fraction-f','FontSize',FS1);
                ylabel('Generalized fractional anisotropy','FontSize',FS1);
                title(['PSNR=',num2str(PSNR),'; GFA with SH, L=',num2str(L)],'FontSize',FS2);
                hold('on');
                grid('on');
                hb2 = myboxplot(FA,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
                hp2(C) = plot(f,mean(FA),'Color',CLS(C,:),'LineStyle','--'); %#ok<SAGROW>
                % -------------------------------------------------------------
                figure(hf3);
                xlabel('Actual value of fraction-f','FontSize',FS1);
                ylabel('Minimum of the ADC','FontSize',FS1);
                title(['PSNR=',num2str(PSNR),'; minimum of the ADC with SH, L=',num2str(L)],'FontSize',FS2);
                hold('on');
                grid('on');
                hb3 = myboxplot(ml,'colors',CLS(C,:),'positions',f+(C-2)*0.011,'labels',labels(1:20),'symbol','r','widths',0.005,'BoxStyle','filled','notch','on');
                hp3(C) = plot(f,mean(ml),'Color',CLS(C,:),'LineStyle','--'); %#ok<SAGROW>
                % -------------------------------------------------------------
                if(C==3)
                    % ---------------------------------------------------------
                    figure(hf2);
                    hl2 = legend(hp2,'1-Compartment','2-Compartments','3-Compartments');
                    set(hl2,'FontSize',FS1);
                    set(hl2,'Location','SouthEast');
                    set(get(hf2,'CurrentAxes'),'FontSize',FS4);
                    % ---------------------------------------------------------
                    figure(hf3);
                    hl3 = legend(hp3,'1-Compartment','2-Compartments','3-Compartments');
                    set(hl3,'FontSize',FS1);
                    set(hl3,'Location','SouthEast');
                    set(get(hf3,'CurrentAxes'),'FontSize',FS4);
                    % ---------------------------------------------------------
                end
                % -------------------------------------------------------------
                % -------------------------------------------------------------
        end
    end
end

end

function h = myboxplot(varargin)
sf = check_software_platform;
if(sf==1)
    h = boxplot(varargin{:});
else
    [h,eraser] = boxplot(varargin{:});
    delete(eraser.outliers);
    delete(eraser.outliers2);
end
end




