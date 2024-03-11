% test_atti2micro
rng(15,'twister');
if(exist('setup__robotics_toolbox','file')~=2)
    error('This test script uses the quaternion class to produce uniform random rotations: http://www.lpi.tel.uva.es/node/626');
end

SAMPLING = 3;
NOISE    = 'rice';      % 'rice'/'gauss'

switch(SAMPLING)
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SYNTHETIC SAMPLING SCHEME:
        [gi,linit,lend] = icosamplesSphere(4);
        G  = size(gi,1);
        bi = zeros(G,1);
        bi(linit(1):lend(1),:) = 500;  % 6 grads
        bi(linit(2):lend(2),:) = 1000; % 10 grads
        bi(linit(3):lend(3),:) = 1750; % 15 grads
        bi(linit(4):lend(4),:) = 2500; % 30 grads
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        %%%%%%%%%%%
        PSNR  = 100;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SYNTHETIC SAMPLING SCHEME (II):
        load('test_atti2micro.mat');
        bi  = bi_synth2;
        gi  = gi_synth2;
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
        % SAMPLING SCHEME TAKEN FORM HCP-MGH 1007
        load('test_atti2micro.mat');
        bi  = bi_HCP;
        gi  = gi_HCP;
        bpt = (abs(bi)>100);
        Nnb = max(length(find(~bpt)),1);
        bi  = bi(bpt);
        gi  = gi(bpt,:);
        G   = size(gi,1);
        %%%%%%%%%%%
        PSNR  = 16;     % Peak signal to noise ratio in the baseline
        sigma = 1/PSNR; % Standard deviation of noise to be added to the DWI
end

N = 200; % Number of synthetic voxels to be tested
DISO = 3.0e-3; % isotropic diffusion coefficient of free water
L = 2; % Order of the spherical harmonics
lambda = 0.006; % Tikhonov regularisation for Spherical Harmonics
tau = 1;

% Precompute random rotations for each compartment;
q1  = [1,ii,jj,kk]*randn(4,N);     % 1xN
ph2 = (rand(1,N)-1/2)*pi/6;        % 1xN, random angle within the range [-15º,15º]
q2  = cos(ph2/2) + kk.*sin(ph2/2); % 1xN, rotation above z in the range [-15º,15º]
ph3 = (rand(1,N)-1/2)*pi/6;        % 1xN, random angle within the range [-15º,15º]
ax3 = randn(3,N);
ax3 = bsxfun( @(x,y)(x./y), ax3, sqrt(sum(ax3.*ax3,1)) ); % 3xN
q3  = cos(ph3/2) + sin(ph3/2).*([ii,jj,kk]*ax3);          % 1xN
q1  = q1./abs(q1);
q2  = q2./abs(q2);
q3  = q3./abs(q3);
q   = {q1,q1.*q2,q1.*q3};
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
            switch(nc)
                case 1
                    xx = q{nc}.*ii.*conj(q{nc}); % 1xN
                    xx = imag(xx); % 3xN
                    yy = q{nc}.*jj.*conj(q{nc}); % 1xN
                    yy = imag(yy); % 3xN
                    zz = q{nc}.*kk.*conj(q{nc}); % 1xN
                    zz = imag(zz); % 3xN
                case 2
                    xx = q{nc}.*jj.*conj(q{nc}); % 1xN
                    xx = imag(xx); % 3xN
                    yy = q{nc}.*kk.*conj(q{nc}); % 1xN
                    yy = imag(yy); % 3xN
                    zz = q{nc}.*ii.*conj(q{nc}); % 1xN
                    zz = imag(zz); % 3xN
                case 3
                    xx = q{nc}.*kk.*conj(q{nc}); % 1xN
                    xx = imag(xx); % 3xN
                    yy = q{nc}.*ii.*conj(q{nc}); % 1xN
                    yy = imag(yy); % 3xN
                    zz = q{nc}.*jj.*conj(q{nc}); % 1xN
                    zz = imag(zz); % 3xN
            end
            % - Generate random eigenvalues for this compartment:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            l1 = 1.55 + 0.2*rand(1,N); % 1xN
            l2 = 0.25 + 0.2*rand(1,N); % 1xN
            l3 = 0.25 + 0.2*rand(1,N); % 1xN
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % - Compute the 6 independent components of the diffusion tensor
            %   corresponding to this compartment:
            DT = zeros(6,N); % 6xN
            DT(1,:) = l1*(xx(1).*xx(1)) + l2*(yy(1).*yy(1)) + l3*(zz(1).*zz(1));
            DT(2,:) = l1*(xx(1).*xx(2)) + l2*(yy(1).*yy(2)) + l3*(zz(1).*zz(2));
            DT(3,:) = l1*(xx(1).*xx(3)) + l2*(yy(1).*yy(3)) + l3*(zz(1).*zz(3));
            DT(4,:) = l1*(xx(2).*xx(2)) + l2*(yy(2).*yy(2)) + l3*(zz(2).*zz(2));
            DT(5,:) = l1*(xx(2).*xx(3)) + l2*(yy(2).*yy(3)) + l3*(zz(2).*zz(3));
            DT(6,:) = l1*(xx(3).*xx(3)) + l2*(yy(3).*yy(3)) + l3*(zz(3).*zz(3));
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
        [~,~,f2] = atti2micro( reshape(Si2',[N,1,1,G]), gi, bi, ...
            'lambda', lambda, ...
            'verbose', false, ...
            'usef', true, ...
            'fmin', 0.1, 'fmax', 1, ...
            'nmax', 100, 'dl', 1.0e-7, 'dC', 1.0e-6, ...
            'ADC0', DISO ...
            );
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
    hf = figure(1001);
    if(C==1)
        close(hf);
        hf = figure(1001);
        plot(f,f,'Color',[0,0,0],'LineWidth',2,'LineStyle','--');
    end
    xlabel('Actual value of fraction-f','FontSize',FS1);
    ylabel('Estimated value of fraction-f','FontSize',FS1);
    title(['PSNR=',num2str(PSNR)],'Interpreter','tex','FontSize',FS2);
    hold('on');
    grid('on');
    hb1 = boxplot(f2cum,'colors',CLS(C,:),'positions',f,'labels',labels(1:20),'symbol','r');
    hp1(C) = plot(f,mean(f2cum),'Color',CLS(C,:),'LineStyle','--');
    if(C==3)
        hl = legend(hp1,'1-Compartment','2-Compartments','3-Compartments');
        set(hl,'FontSize',FS1);
        set(hl,'Location','SouthEast');
        set(get(hf,'CurrentAxes'),'FontSize',FS4);
        axis([0,1,-0.05,1.05]);
    end
end