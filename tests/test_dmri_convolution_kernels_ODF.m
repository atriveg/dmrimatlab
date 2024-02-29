% test_dmri_convolution_kernels_ODF.m
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf(1,'Design experiment... ');
ref   = 7;
G     = icosamplesSphere(ref,'O1',true);
bi    = [500;1000;2000;3000;5000];
N     = 20;
lpar  = linspace(1.0e-3,3.0e-3,N);
lperp = linspace(0.2e-3,0.9e-3,N);
[lpar,lperp] = meshgrid(lpar,lperp);
lpar  = lpar(:);
lperp = lperp(:);
tau   = 20.0e-3;
L     = 8;
lmb   = 0;
fprintf(1,'Done\n');
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
fprintf(1,'Compute analytical weights... ');
elODF  = dmri_compute_convolution_weights_ODF(bi,lpar,lperp,L);
fprintf(1,'Done\n');
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cuv2   = G(:,3).*G(:,3);   % Gx1
cuv2   = cuv2';            % 1xG
delta  = lpar-lperp;       % N^2x1
rho    = lperp./lpar;      % N^2x1
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
fprintf(1,'Numerically compute ODF weights... ');
for nb=1:length(bi)
    numis = zeros(N*N,L/2+1);
    ptr = 1:N;
    for n=1:N
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        arg = delta(ptr,:)*cuv2;                            % NxG
        if(is_broadcast_available)
            arg = arg + lperp(ptr,:);                       % NxG
        else
            arg = bsxfun( @(x,y)(x+y), arg, lperp(ptr,:) ); % NxG
        end
        arg = exp(-bi(nb)*arg); % NxG
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        arg = arg*SH;           % Nx(L/2+1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numis(ptr,:) = arg;     % Nx(L/2+1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ptr = ptr+N;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    hf = figure(1000+nb);
    close(hf);
    figure(1000+nb);
    grid('on');
    hold('on');
    for nl=1:L/2+1
        plot(delta,numis(:,nl)./exp(-bi(nb)*lperp),'o','Color',cols(nl,:));
        plot(delta,elODF(:,nl,nb)./exp(-bi(nb)*lperp),'*','Color',cols(nl,:));
    end
    xlabel('\delta_{\lambda}');
    ylabel('e_{l,ODF}^i/exp(-b_i\lambda_{\perp})');
    title(['ODF convolution coefficients at b_i = ',num2str(bi(nb)),'   ---   *=analytical; o=numerical']);
end
fprintf(1,'Done\n');
toc;
