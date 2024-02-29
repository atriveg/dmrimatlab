% test_dmri_convolution_kernels_momentsEbis.m
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------

%%% -----------------------------------------------------------------------
gammas = [0,1,2,3,7,8,1/4,pi];
Nx     = 20000;
x1     = (0:Nx-1)/Nx;
x2     = x1.*x1;
N      = 4;
rho    = [0.000001,0.00001,0.0001, ...
    0.001,0.01,0.1,0.2,0.5,0.7, ...
    0.8,0.9,1.0,2.0,3.0,5.0,10.0, ...
    20.0, 30.0, 40.0, 50.0, inf ]';      % Mx1
lperp  = 0.1e-3*ones(size(rho));         % Mx1
lpar   = lperp.*(1+rho)./rho;            % Mx1
lpar(isinf(rho)) = lperp(isinf(rho));    % Mx1
%%% -----------------------------------------------------------------------
MS   = 18;
cols = [ ...
    .7, .0, .0;
    .0, .7, .0;
    .0, .0, .7;
    .5, .5, .0;
    .5, .0, .5;
    .0, .5, .5;
    ];
NC = 6;
%%% -----------------------------------------------------------------------
for g = 1:length(gammas)
    %%% -------------------------------------------------------------------
    gamma = gammas(g);
    %%% -------------------------------------------------------------------
    IOTA = dmri_compute_Ekernel_integrals(lpar,lperp,gamma,N); % MxN
    IOTA = dmri_compute_Legendre_projection(IOTA);             % MxN
    %%% -------------------------------------------------------------------
    hf = figure(1000+g);
    close(hf);
    hf = figure(1000+g);
    hold('on');
    grid('on');
    xlabel('rho_{\lambda}');
    ylabel('I_{\gamma}^n');
    title(['\gamma = ',num2str(gamma)]);
    %%% -------------------------------------------------------------------
    hl    = zeros(1,N+1);
    hlt   = cell(1,N+1);
    rhopl = [ rho(1:end-1); 60 ];
    for n=0:N
        rhof = rho;
        rhof(isinf(rho)) = 0.001/eps;
        f  = (1./rhof)*x2 + 1;  % MxNx
        f  = f.^(-gamma/2);     % MxNx
        LP = legendre(2*n,x1);  % (4*n+1)xNx
        LP = LP(1,:);           % 1xNx
        if(is_broadcast_available)
            f = LP.*f;                         % MxNx
        else
            f = bsxfun( @(x,y)(x.*y), LP, f ); % MxNx
        end
        IOTAn = 4*pi*sum(f,2)*(x1(2)-x1(1)); % Mx1
        figure(hf);
        hl(n+1) = plot( rhopl, IOTAn, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MS, 'Color', cols(rem(n,NC)+1,:) );
        plot( rhopl, IOTA(:,n+1), 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', MS, 'Color', cols(rem(n,NC)+1,:) );
        hlt{n+1} = ['l = ',num2str(n)];
    end
    figure(hf);
    legend(hl,hlt{:});
    hf.CurrentAxes.XTick = [0,10,20,30,40,50,60];
    hf.CurrentAxes.XTickLabel = {'0','10','20','30','40','50','inf'};
    %%% -------------------------------------------------------------------
end
%%% -----------------------------------------------------------------------
