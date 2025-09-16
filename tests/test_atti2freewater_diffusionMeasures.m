sf = check_software_platform;
if(sf==2)
    pkg load statistics;
end

rng(15,'twister');

% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------

SAMPLING = 3;
NOISE    = 'rice';      % 'rice'/'gauss'
switch(SAMPLING)
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
        PSNR  = 10;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAMPLING SCHEME TAKEN FROM HCP-MGH 1007
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
        PSNR  = 16;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
    case 3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SAMPLING SCHEME TAKEN FROM PDD case CO_p07090
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
        PSNR  = 22;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

N = 200; % Number of synthetic voxels to be tested
DISO = 3.0e-3; % isotropic diffusion coefficient of free water
L = 6; % Order of the spherical harmonics
lambda = 0.006; % Tikhonov regularisation for Spherical Harmonics
tau = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precompute random rotations for each compartment
q1  = randn(N,4);
nq1 = sqrt(sum(q1.*q1,2));
q1  = q1./nq1;
%%%
ph2 = (rand(N,1)-1/2)*pi/18;
q2  = [ cos(ph2/2), zeros(N,1), zeros(N,1), sin(ph2/2) ];
%%%
ph3 = (rand(N,1)-1/2)*pi/9;
ax3 = randn(N,3);
nax = sqrt(sum(ax3.*ax3,2));
ax3 = ax3./nax;
q3  = [ cos(ph3/2), ax3.*sin(ph3/2) ];
%%%
q   = { q1, multiply_quatrot(q1,q2), multiply_quatrot(q1,q3) };
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f   = (1:1:20)/20;
hp1 = zeros(1,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau_D = 20.0e-3;
BSH   = GenerateSHMatrix(L,gi); % G x (L+1)(L+2)/2
LSSH  = (BSH'*BSH+lambda*eye(size(BSH,2)))\(BSH'); % (L+1)(L+2)/2 x G
FRT   = GenerateFRTMatrix(L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
apacnst  = 0.4;
rcnst    = sqrt(4*pi*pi*tau_D);
rtopcnst = sqrt(4*pi)*sqrt(pi)/(4*rcnst*rcnst*rcnst);
rtppcnst = sqrt(pi)/rcnst;
rtapcnst = 1/(2*rcnst*rcnst);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for C=1:3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RTOP_true = zeros(N,length(f));
    RTPP_true = zeros(N,length(f));
    RTAP_true = zeros(N,length(f));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PA_estim   = zeros(N,length(f));
    GFA_estim  = zeros(N,length(f));
    RTOP_estim = zeros(N,length(f));
    RTPP_estim = zeros(N,length(f));
    RTAP_estim = zeros(N,length(f));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PA_nowater   = zeros(N,length(f));
    GFA_nowater  = zeros(N,length(f));
    RTOP_nowater = zeros(N,length(f));
    RTPP_nowater = zeros(N,length(f));
    RTAP_nowater = zeros(N,length(f));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PA_water   = zeros(N,length(f));
    GFA_water  = zeros(N,length(f));
    RTOP_water = zeros(N,length(f));
    RTPP_water = zeros(N,length(f));
    RTAP_water = zeros(N,length(f));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PA_corr   = zeros(N,length(f));
    GFA_corr  = zeros(N,length(f));
    RTOP_corr = zeros(N,length(f));
    RTPP_corr = zeros(N,length(f));
    RTAP_corr = zeros(N,length(f));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
            l1 = 1.6 + 0.1*rand(1,N); % 1xN
            l2 = 0.3 + 0.1*rand(1,N); % 1xN
            l3 = 0.3 + 0.1*rand(1,N); % 1xN
            l1 = l1*1.0e-3;
            l2 = l2*1.0e-3;
            l3 = l3*1.0e-3;
            % - Compute the 6 independent components of the diffusion tensor
            %   corresponding to this compartment:
            DT = zeros(6,N); % 6xN
            DT(1,:) = l1*(xx(1).*xx(1)) + l2*(yy(1).*yy(1)) + l3*(zz(1).*zz(1));
            DT(2,:) = l1*(xx(1).*xx(2)) + l2*(yy(1).*yy(2)) + l3*(zz(1).*zz(2));
            DT(3,:) = l1*(xx(1).*xx(3)) + l2*(yy(1).*yy(3)) + l3*(zz(1).*zz(3));
            DT(4,:) = l1*(xx(2).*xx(2)) + l2*(yy(2).*yy(2)) + l3*(zz(2).*zz(2));
            DT(5,:) = l1*(xx(2).*xx(3)) + l2*(yy(2).*yy(3)) + l3*(zz(2).*zz(3));
            DT(6,:) = l1*(xx(3).*xx(3)) + l2*(yy(3).*yy(3)) + l3*(zz(3).*zz(3));
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
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            tens_const = sqrt(4*pi*tau_D);
            RTOP_true(:,nf) = RTOP_true(:,nf) +  (fi(nc,:)./sqrt(l1.*l2.*l3))'/(tens_const*tens_const*tens_const);
            RTPP_true(:,nf) = RTPP_true(:,nf) +  (fi(nc,:)./sqrt(l1))'/(tens_const);
            RTAP_true(:,nf) = RTAP_true(:,nf) +  (fi(nc,:)./sqrt(l2.*l3))'/(tens_const*tens_const);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        % -- Divide by the total volume fraction:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        RTOP_true(:,nf) = RTOP_true(:,nf)./(sum(fi,1)');
        RTPP_true(:,nf) = RTPP_true(:,nf)./(sum(fi,1)');
        RTAP_true(:,nf) = RTAP_true(:,nf)./(sum(fi,1)');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Si = Si./repmat(sum(fi,1),[G,1]); % GxN
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SHADC_estim = atti2shadc( reshape(Si',[size(Si,2),1,1,size(Si,1)]), ...
            gi, bi, ...
            'L', L, 'lambda', lambda, 'tl', 10*eps, 'tu', 1-10*eps ); % N x 1 x 1 x (L+1)(L+2)
        ADC_estim   = sh2signal( SHADC_estim, gi );    % N x 1 x 1 x G
        SHADC_estim = permute(SHADC_estim,[4,1,2,3]);  % (L+1)(L+2)/2 x N
        ADC_estim   = permute(ADC_estim,[4,1,2,3]);    % G x N
        [~ ,r0]     = max(ADC_estim,[],1);
        %%% --------------
        Si_temp = Si;                       % GxN
        Si_temp(Si<10*eps)   = 10*eps;      % GxN
        Si_temp(Si>1-10*eps) = 1-10*eps;    % GxN
        ADC_temp = -log(Si_temp);           % GxN
        if(is_broadcast_available)
            ADC_temp = ADC_temp./bi;                         % GxN
        else
            ADC_temp = bsxfun( @(x,y)(x./y), ADC_temp, bi ); % GxN
        end
        %%% --------------
        GFA_estim(:,nf)  = sqrt( sum(SHADC_estim(2:end,:).*SHADC_estim(2:end,:),1)./sum(SHADC_estim.*SHADC_estim,1) )';
        %%% --------------
        SH_temp = signal2sh( reshape((ADC_temp.^(-3/2))',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda ); % N x 1 x 1 x (L+1)(L+2)/2
        RTOP_estim(:,nf) = rtopcnst*SH_temp(:,1,1,1)';
        %%% --------------
        RTPP_estim(:,nf) = rtppcnst./sqrt(ADC_estim(sub2ind(size(ADC_estim),r0,1:N))');
        %%% --------------
        SH_temp  = signal2sh( reshape((1./ADC_temp)',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );         % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = FRT*permute(SH_temp,[4,1,2,3]); % (L+1)(L+2)/2 x N
        sig_temp = sh2signal( permute(SH_temp,[2,3,4,1]), gi ); % N x 1 x 1 x G
        sig_temp = permute( sig_temp, [4,2,3,1] ); % G x N
        RTAP_estim(:,nf) = rtapcnst.*(sig_temp(sub2ind(size(ADC_estim),r0,1:N))');
        %%% --------------
        SH_temp  = signal2sh( reshape(ADC_temp',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );              % N x 1 x 1 x (L+1)(L+2)/2
        ADCAV    = permute(SH_temp(:,1,1,1),[4,1,2,3]); % 1xN
        ADCAV    = ADCAV/sqrt(4*pi);
        if(is_broadcast_available)
            ADC2 = ADC_temp + ADCAV;                         % GxN
        else
            ADC2 = bsxfun( @(x,y)(x+y), ADC_temp, ADCAV );   % GxN
        end
        ADC2     = ADC2.^(-3/2);
        SH_temp  = signal2sh( reshape(ADC2',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = permute(SH_temp(:,1,1,1),[1,4,2,3]);      % N x 1
        ADC2     = (ADC_temp).^(-3/2);
        SH_temp2 = signal2sh( reshape(ADC2',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp2 = permute(SH_temp2(:,1,1,1),[1,4,2,3]);     % N x 1
        COSPQ    = (4/sqrt(pi))*(SH_temp.*SH_temp)./(SH_temp2.*(ADCAV').^(-3/2));
        SINPQ    = sqrt(1-COSPQ.*COSPQ);
        PA_estim(:,nf) = SINPQ.^(3*apacnst)./(1-3*SINPQ.^apacnst+3*SINPQ.^(2*apacnst));
        %%% --------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -- Compute the overall signal accounting for the free water:
        Si2 = repmat((1-f(nf))*exp(-DISO*bi),[1,N]) + f(nf)*Si; % GxN
        
        % -- Add Rician/Gaussian noise to both the DWI and the baseline
        switch(NOISE)
            case 'rice'
                ni  = sigma*(randn(G,N)+1i*randn(G,N));           % GxN
                n0  = sigma*(randn(1,N)+1i*randn(1,N))/sqrt(Nnb); % GxN
                Si2n = abs(Si2+ni)./repmat(abs(1+n0),[G,1]);      % GxN
                Sin  = abs(Si+ni)./repmat(abs(1+n0),[G,1]);       % GxN
            case 'gauss'
                Si2n = Si2 + sigma*(randn(G,N)+randn(G,N));  % GxN
                Sin  = Si + sigma*(randn(G,N)+randn(G,N));   % GxN
            otherwise
                error('NOISE must be either ''gauss'' or ''rice''')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SHADC_nowater = atti2shadc( reshape(Sin',[size(Sin,2),1,1,size(Sin,1)]), ...
            gi, bi, ...
            'L', L, 'lambda', lambda, 'tl', 10*eps, 'tu', 1-10*eps ); % N x 1 x 1 x (L+1)(L+2)
        ADC_nowater   = sh2signal( SHADC_nowater, gi );    % N x 1 x 1 x G
        SHADC_nowater = permute(SHADC_nowater,[4,1,2,3]);  % (L+1)(L+2)/2 x N
        ADC_nowater   = permute(ADC_nowater,[4,1,2,3]);    % G x N
        [~ ,r0]       = max(ADC_nowater,[],1);
        %%% --------------
        Si_temp = Sin;                      % GxN
        Si_temp(Sin<10*eps)   = 10*eps;     % GxN
        Si_temp(Sin>1-10*eps) = 1-10*eps;   % GxN
        ADC_temp = -log(Si_temp);           % GxN
        if(is_broadcast_available)
            ADC_temp = ADC_temp./bi;
        else
            ADC_temp = bsxfun( @(x,y)(x./y), ADC_temp, bi );
        end
        %%% --------------
        GFA_nowater(:,nf)  = sqrt( sum(SHADC_nowater(2:end,:).*SHADC_nowater(2:end,:),1)./sum(SHADC_nowater.*SHADC_nowater,1) )';
        %%% --------------
        SH_temp = signal2sh( reshape((ADC_temp.^(-3/2))',[size(Sin,2),1,1,size(Sin,1)]), ...
            gi, 'L', L, 'lambda',lambda ); % N x 1 x 1 x (L+1)(L+2)/2
        RTOP_nowater(:,nf) = rtopcnst*SH_temp(:,1,1,1)';
        %%% --------------
        RTPP_nowater(:,nf) = rtppcnst./sqrt(ADC_nowater(sub2ind(size(ADC_nowater),r0,1:N))');
        %%% --------------
        SH_temp  = signal2sh( reshape((1./ADC_temp)',[size(Sin,2),1,1,size(Sin,1)]), ...
            gi, 'L', L, 'lambda',lambda );        % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = FRT*permute(SH_temp,[4,1,2,3]); % (L+1)(L+2)/2 x N
        sig_temp = sh2signal( permute(SH_temp,[2,3,4,1]), gi ); % N x 1 x 1 x G
        sig_temp = permute( sig_temp, [4,2,3,1] ); % G x N
        RTAP_nowater(:,nf) = rtapcnst.*(sig_temp(sub2ind(size(ADC_nowater),r0,1:N))');
        %%% --------------
        SH_temp  = signal2sh( reshape(ADC_temp',[size(Sin,2),1,1,size(Sin,1)]), ...
            gi, 'L', L, 'lambda',lambda );              % N x 1 x 1 x (L+1)(L+2)/2
        ADCAV    = permute(SH_temp(:,1,1,1),[4,1,2,3]); % 1xN
        ADCAV    = ADCAV/sqrt(4*pi);
        if(is_broadcast_available)
            ADC2 = ADC_temp + ADCAV;                         % GxN
        else
            ADC2 = bsxfun( @(x,y)(x+y), ADC_temp, ADCAV );   % GxN
        end
        ADC2     = ADC2.^(-3/2);
        SH_temp  = signal2sh( reshape(ADC2',[size(Sin,2),1,1,size(Sin,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = permute(SH_temp(:,1,1,1),[1,4,2,3]);      % N x 1
        ADC2     = (ADC_temp).^(-3/2);
        SH_temp2 = signal2sh( reshape(ADC2',[size(Sin,2),1,1,size(Sin,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp2 = permute(SH_temp2(:,1,1,1),[1,4,2,3]);     % N x 1
        COSPQ    = (4/sqrt(pi))*(SH_temp.*SH_temp)./(SH_temp2.*(ADCAV').^(-3/2));
        SINPQ    = sqrt(1-COSPQ.*COSPQ);
        PA_nowater(:,nf) = SINPQ.^(3*apacnst)./(1-3*SINPQ.^apacnst+3*SINPQ.^(2*apacnst));
        %%% --------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SHADC_water = atti2shadc( reshape(Si2n',[size(Si2n,2),1,1,size(Si2n,1)]), ...
            gi, bi, ...
            'L', L, 'lambda', lambda, 'tl', 10*eps, 'tu', 1-10*eps ); % N x 1 x 1 x (L+1)(L+2)
        ADC_water   = sh2signal( SHADC_water, gi );    % N x 1 x 1 x G
        SHADC_water = permute(SHADC_water,[4,1,2,3]);  % (L+1)(L+2)/2 x N
        ADC_water   = permute(ADC_water,[4,1,2,3]);    % G x N
        [~ ,r0]     = max(ADC_water,[],1);
        %%% --------------
        Si_temp = Si2n;                      % GxN
        Si_temp(Si2n<10*eps)   = 10*eps;     % GxN
        Si_temp(Si2n>1-10*eps) = 1-10*eps;   % GxN
        ADC_temp = -log(Si_temp);            % GxN
        if(is_broadcast_available)
            ADC_temp = ADC_temp./bi;
        else
            ADC_temp = bsxfun( @(x,y)(x./y), ADC_temp, bi );
        end
        %%% --------------
        GFA_water(:,nf)  = sqrt( sum(SHADC_water(2:end,:).*SHADC_water(2:end,:),1)./sum(SHADC_water.*SHADC_water,1) )';
        %%% --------------
        SH_temp = signal2sh( reshape((ADC_temp.^(-3/2))',[size(Si2n,2),1,1,size(Si2n,1)]), ...
            gi, 'L', L, 'lambda',lambda ); % N x 1 x 1 x (L+1)(L+2)/2
        RTOP_water(:,nf) = rtopcnst*SH_temp(:,1,1,1)';
        %%% --------------
        RTPP_water(:,nf) = rtppcnst./sqrt(ADC_water(sub2ind(size(ADC_water),r0,1:N))');
        %%% --------------
        SH_temp  = signal2sh( reshape((1./ADC_temp)',[size(Si2n,2),1,1,size(Si2n,1)]), ...
            gi, 'L', L, 'lambda',lambda );        % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = FRT*permute(SH_temp,[4,1,2,3]); % (L+1)(L+2)/2 x N
        sig_temp = sh2signal( permute(SH_temp,[2,3,4,1]), gi ); % N x 1 x 1 x G
        sig_temp = permute( sig_temp, [4,2,3,1] ); % G x N
        RTAP_water(:,nf) = rtapcnst.*(sig_temp(sub2ind(size(ADC_water),r0,1:N))');
        %%% --------------
        SH_temp  = signal2sh( reshape(ADC_temp',[size(Si2n,2),1,1,size(Si2n,1)]), ...
            gi, 'L', L, 'lambda',lambda );              % N x 1 x 1 x (L+1)(L+2)/2
        ADCAV    = permute(SH_temp(:,1,1,1),[4,1,2,3]); % 1xN
        ADCAV    = ADCAV/sqrt(4*pi);
        if(is_broadcast_available)
            ADC2 = ADC_temp + ADCAV;                         % GxN
        else
            ADC2 = bsxfun( @(x,y)(x+y), ADC_temp, ADCAV );   % GxN
        end
        ADC2     = ADC2.^(-3/2);
        SH_temp  = signal2sh( reshape(ADC2',[size(Si2n,2),1,1,size(Si2n,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = permute(SH_temp(:,1,1,1),[1,4,2,3]);      % N x 1
        ADC2     = (ADC_temp).^(-3/2);
        SH_temp2 = signal2sh( reshape(ADC2',[size(Si2n,2),1,1,size(Si2n,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp2 = permute(SH_temp2(:,1,1,1),[1,4,2,3]);     % N x 1
        COSPQ    = (4/sqrt(pi))*(SH_temp.*SH_temp)./(SH_temp2.*(ADCAV').^(-3/2));
        SINPQ    = sqrt(1-COSPQ.*COSPQ);
        PA_water(:,nf) = SINPQ.^(3*apacnst)./(1-3*SINPQ.^apacnst+3*SINPQ.^(2*apacnst));
        %%% --------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -- Run the algorithm:
        [f_corr,SH_corr] = atti2freewater( reshape(Si2n',[N,1,1,G]), gi, bi, ...
            'L', L, 'lambda', lambda, 't0', 0.00001, ...
            'verbose', false, ...
            'O1', true, ...
            'O2', true, 'tau', tau, 'fnriters', 100, 'fnrth', 0.000001, 'lmb', 1, ...
            'ADC0', DISO, 'pos', false );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Si2n has size G x N
        % bi has size G x 1
        %%% ----------------------
        SH_corr  = permute(SH_corr,[4,1,2,3]); % (L+1)(L+2)/2 x N
        ADC_corr = BSH*SH_corr;        % G x N
        ADC_corr = ADC_corr.*ADC_corr; % G x N
        %%% ----------------------
        Dm32     = sqrt(ADC_corr);     % G x N
        Dm32     = Dm32.*Dm32.*Dm32;   % G x N
        rtop     = LSSH*(1./(Dm32+1000*eps)); % (L+1)(L+2)/2 x N
        rtop     = rtopcnst*rtop(1,:);        % 1xN
        RTOP_corr(:,nf) = rtop';
        %%% ----------------------
        SH_corr  = LSSH*ADC_corr;
        GFA_corr(:,nf) = sqrt( sum(SH_corr(2:end,:).*SH_corr(2:end,:),1)./sum(SH_corr.*SH_corr,1) )';
        %%% ----------------------
        [~ ,r0] = max(ADC_corr,[],1);
        RTPP_corr(:,nf) = rtppcnst./sqrt(ADC_corr(sub2ind(size(ADC_corr),r0,1:N))');
        %%% ----------------------
        SH_temp  = signal2sh( reshape((1./ADC_corr)',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );        % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = FRT*permute(SH_temp,[4,1,2,3]); % (L+1)(L+2)/2 x N
        sig_temp = sh2signal( permute(SH_temp,[2,3,4,1]), gi ); % N x 1 x 1 x G
        sig_temp = permute( sig_temp, [4,2,3,1] ); % G x N
        RTAP_corr(:,nf) = rtapcnst.*(sig_temp(sub2ind(size(ADC_corr),r0,1:N))');
        %%% --------------
        SH_temp  = signal2sh( reshape(ADC_corr',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );              % N x 1 x 1 x (L+1)(L+2)/2
        ADCAV    = permute(SH_temp(:,1,1,1),[4,1,2,3]); % 1xN
        ADCAV    = ADCAV/sqrt(4*pi);
        if(is_broadcast_available)
            ADC2 = ADC_corr + ADCAV;                         % GxN
        else
            ADC2 = bsxfun( @(x,y)(x+y), ADC_corr, ADCAV );   % GxN
        end
        ADC2     = ADC2.^(-3/2);
        SH_temp  = signal2sh( reshape(ADC2',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp  = permute(SH_temp(:,1,1,1),[1,4,2,3]);      % N x 1
        ADC2     = (ADC_corr).^(-3/2);
        SH_temp2 = signal2sh( reshape(ADC2',[size(Si,2),1,1,size(Si,1)]), ...
            gi, 'L', L, 'lambda',lambda );                   % N x 1 x 1 x (L+1)(L+2)/2
        SH_temp2 = permute(SH_temp2(:,1,1,1),[1,4,2,3]);     % N x 1
        COSPQ    = (4/sqrt(pi))*(SH_temp.*SH_temp)./(SH_temp2.*(ADCAV').^(-3/2));
        if(any(COSPQ>1))
            stophere = true;
        end
        SINPQ    = sqrt(1-COSPQ.*COSPQ);
        PA_corr(:,nf) = SINPQ.^(3*apacnst)./(1-3*SINPQ.^apacnst+3*SINPQ.^(2*apacnst));
        %%% --------------
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    fprintf(1,repmat('\b',[1,txtl]));
    fprintf(1,'\n');
    
    FS1 = 14;
    FS2 = 15;
    FS3 = 14;
    FS4 = 12;
    MSZ = 14;
    Colors = [.0,.5,.0;
        .0,.0,.5;
        .0,.3,.3;
        .5,.0,.0;
        .3,.3,.0];
    hf(1) = figure(1001);
    hf(2) = figure(2002);
    hf(3) = figure(3003);
    hf(4) = figure(4004);
    hf(5) = figure(5005);
    if(C==1)
        close(hf);
        hf(1) = figure(1001);
        hf(2) = figure(2002);
        hf(3) = figure(3003);
        hf(4) = figure(4004);
        hf(5) = figure(5005);
    end
    % -------------------------------------------------------------
    % ------
    figure(hf(1)); subplot(1,3,C);
    title(['GFA: PSNR=',num2str(PSNR),'; L=',num2str(L),'; \lambda=',num2str(lambda)],'Interpreter','tex','FontSize',FS2);
    xlabel('Actual value of fraction-f','FontSize',FS1);
    ylabel('GFA','FontSize',FS1);
    hold('on'); grid('on');
    % ------
    figure(hf(2)); subplot(1,3,C);
    title(['RTOP: PSNR=',num2str(PSNR),'; L=',num2str(L),'; \lambda=',num2str(lambda)],'Interpreter','tex','FontSize',FS2);
    xlabel('Actual value of fraction-f','FontSize',FS1);
    ylabel('RTOP','FontSize',FS1);
    hold('on'); grid('on');
    % ------
    figure(hf(3)); subplot(1,3,C);
    title(['RTPP: PSNR=',num2str(PSNR),'; L=',num2str(L),'; \lambda=',num2str(lambda)],'Interpreter','tex','FontSize',FS2);
    xlabel('Actual value of fraction-f','FontSize',FS1);
    ylabel('RTPP','FontSize',FS1);
    hold('on'); grid('on');
    % ------
    figure(hf(4)); subplot(1,3,C);
    title(['RTAP: PSNR=',num2str(PSNR),'; L=',num2str(L),'; \lambda=',num2str(lambda)],'Interpreter','tex','FontSize',FS2);
    xlabel('Actual value of fraction-f','FontSize',FS1);
    ylabel('RTAP','FontSize',FS1);
    hold('on'); grid('on');
    % ------
    figure(hf(5)); subplot(1,3,C);
    title(['PA: PSNR=',num2str(PSNR),'; L=',num2str(L),'; \lambda=',num2str(lambda)],'Interpreter','tex','FontSize',FS2);
    xlabel('Actual value of fraction-f','FontSize',FS1);
    ylabel('PA','FontSize',FS1);
    hold('on'); grid('on');
    % ------
    % -------------------------------------------------------------
    % ------
    figure(hf(1)); subplot(1,3,C); hold('on');
    hl1(C,1) = plot(f,mean(GFA_estim,1),  'LineWidth',2,'LineStyle',':', 'Color',Colors(1,:),'Marker','o','MarkerSize',MSZ); %#ok<SAGROW>
    hl1(C,2) = plot(f,mean(GFA_nowater,1),'LineWidth',2,'LineStyle','-.','Color',Colors(2,:),'Marker','x','MarkerSize',MSZ); %#ok<SAGROW>
    hl1(C,3) = plot(f,mean(GFA_water,1),  'LineWidth',2,'LineStyle','-.','Color',Colors(3,:),'Marker','^','MarkerSize',MSZ); %#ok<SAGROW>
    hl1(C,4) = plot(f,mean(GFA_corr,1),   'LineWidth',2,'LineStyle','-', 'Color',Colors(4,:),'Marker','v','MarkerSize',MSZ); %#ok<SAGROW>
    hl1 = legend(hl1(C,:),'estimated','measured, no-water','measured','corrected');
    set(hl1,'FontSize',FS1);
    set(get(hf(1),'CurrentAxes'),'FontSize',FS4);
    axis([0,1,0,0.5]);
    % ------
    figure(hf(2)); subplot(1,3,C); hold('on');
    hl2(C,1) = plot(f,mean(RTOP_estim,1),  'LineWidth',2,'LineStyle',':', 'Color',Colors(1,:),'Marker','o','MarkerSize',MSZ); %#ok<SAGROW>
    hl2(C,2) = plot(f,mean(RTOP_nowater,1),'LineWidth',2,'LineStyle','-.','Color',Colors(2,:),'Marker','x','MarkerSize',MSZ); %#ok<SAGROW>
    hl2(C,3) = plot(f,mean(RTOP_water,1),  'LineWidth',2,'LineStyle','-.','Color',Colors(3,:),'Marker','^','MarkerSize',MSZ); %#ok<SAGROW>
    hl2(C,4) = plot(f,mean(RTOP_corr,1),   'LineWidth',2,'LineStyle','-', 'Color',Colors(4,:),'Marker','v','MarkerSize',MSZ); %#ok<SAGROW>
    hl2(C,5) = plot(f,mean(RTOP_true,1),   'LineWidth',2,'LineStyle','-', 'Color',Colors(5,:),'Marker','*','MarkerSize',MSZ); %#ok<SAGROW>
    hl2 = legend(hl2(C,:),'estimated','measured, no-water','measured','corrected','true');
    set(hl2,'FontSize',FS1);
    set(get(hf(2),'CurrentAxes'),'FontSize',FS4);
    axis([0,1,min([min(RTOP_water(:)),min(RTOP_estim(:)),min(RTOP_nowater(:))]),max([max(RTOP_water(:)),max(RTOP_estim(:)),max(RTOP_nowater(:))])]);
    % ------
    figure(hf(3)); subplot(1,3,C); hold('on');
    hl3(C,1) = plot(f,mean(RTPP_estim,1),  'LineWidth',2,'LineStyle',':', 'Color',Colors(1,:),'Marker','o','MarkerSize',MSZ); %#ok<SAGROW>
    hl3(C,2) = plot(f,mean(RTPP_nowater,1),'LineWidth',2,'LineStyle','-.','Color',Colors(2,:),'Marker','x','MarkerSize',MSZ); %#ok<SAGROW>
    hl3(C,3) = plot(f,mean(RTPP_water,1),  'LineWidth',2,'LineStyle','-.','Color',Colors(3,:),'Marker','^','MarkerSize',MSZ); %#ok<SAGROW>
    hl3(C,4) = plot(f,mean(RTPP_corr,1),   'LineWidth',2,'LineStyle','-', 'Color',Colors(4,:),'Marker','v','MarkerSize',MSZ); %#ok<SAGROW>
    hl3(C,5) = plot(f,mean(RTPP_true,1),   'LineWidth',2,'LineStyle','-', 'Color',Colors(5,:),'Marker','*','MarkerSize',MSZ); %#ok<SAGROW>
    hl3 = legend(hl3(C,:),'estimated','measured, no-water','measured','corrected','true');
    set(hl3,'FontSize',FS1);
    set(get(hf(3),'CurrentAxes'),'FontSize',FS4);
    axis([0,1,min([min(RTPP_water(:)),min(RTPP_estim(:)),min(RTPP_nowater(:))]),max([max(RTPP_water(:)),max(RTPP_estim(:)),max(RTPP_nowater(:))])]);
    % ------
    figure(hf(4)); subplot(1,3,C); hold('on');
    hl4(C,1) = plot(f,mean(RTAP_estim,1),  'LineWidth',2,'LineStyle',':', 'Color',Colors(1,:),'Marker','o','MarkerSize',MSZ); %#ok<SAGROW>
    hl4(C,2) = plot(f,mean(RTAP_nowater,1),'LineWidth',2,'LineStyle','-.','Color',Colors(2,:),'Marker','x','MarkerSize',MSZ); %#ok<SAGROW>
    hl4(C,3) = plot(f,mean(RTAP_water,1),  'LineWidth',2,'LineStyle','-.','Color',Colors(3,:),'Marker','^','MarkerSize',MSZ); %#ok<SAGROW>
    hl4(C,4) = plot(f,mean(RTAP_corr,1),   'LineWidth',2,'LineStyle','-', 'Color',Colors(4,:),'Marker','v','MarkerSize',MSZ); %#ok<SAGROW>
    hl4(C,5) = plot(f,mean(RTAP_true,1),   'LineWidth',2,'LineStyle','-' ,'Color',Colors(5,:),'Marker','*','MarkerSize',MSZ); %#ok<SAGROW>
    hl4 = legend(hl4(C,:),'estimated','measured, no-water','measured','corrected','true');
    set(hl4,'FontSize',FS1);
    set(get(hf(4),'CurrentAxes'),'FontSize',FS4);
    axis([0,1,min([min(RTAP_water(:)),min(RTAP_estim(:)),min(RTAP_nowater(:))]),max([max(RTAP_water(:)),max(RTAP_estim(:)),max(RTAP_nowater(:))])]);
    % ------
    figure(hf(5)); subplot(1,3,C); hold('on');
    hl5(C,1) = plot(f,mean(PA_estim,1),  'LineWidth',2,'LineStyle',':', 'Color',Colors(1,:),'Marker','o','MarkerSize',MSZ); %#ok<SAGROW>
    hl5(C,2) = plot(f,mean(PA_nowater,1),'LineWidth',2,'LineStyle','-.','Color',Colors(2,:),'Marker','x','MarkerSize',MSZ); %#ok<SAGROW>
    hl5(C,3) = plot(f,mean(PA_water,1),  'LineWidth',2,'LineStyle','-.','Color',Colors(3,:),'Marker','^','MarkerSize',MSZ); %#ok<SAGROW>
    hl5(C,4) = plot(f,mean(PA_corr,1),   'LineWidth',2,'LineStyle','-', 'Color',Colors(4,:),'Marker','v','MarkerSize',MSZ); %#ok<SAGROW>
    hl5 = legend(hl5(C,:),'estimated','measured, no-water','measured','corrected');
    set(hl1,'FontSize',FS1);
    set(get(hf(1),'CurrentAxes'),'FontSize',FS4);
    axis([0,1,0,1]);
end




