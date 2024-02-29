% test_dmri_convolution_kernels_momentsP.m
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
kappa  = [1.000001,1.00001,1.0001,1.001,1.01,1.1, ...
    1.5,2.0,3.0,5.0,10.0, ...
    20.0,30.0,40.0,50.0]';               % Mx1
lperp  = 0.1e-3*ones(size(kappa));       % Mx1
lpar   = lperp.*kappa./(kappa-1);        % Mx1
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
    iota = dmri_compute_Pkernel_integrals(lpar,lperp,gamma,N); % MxN
    %%% -------------------------------------------------------------------
    hf = figure(1000+g);
    close(hf);
    hf = figure(1000+g);
    hold('on');
    grid('on');
    xlabel('kappa_{\lambda}');
    ylabel('I_{\gamma}^n');
    title(['\gamma = ',num2str(gamma)]);
    %%% -------------------------------------------------------------------
    x   = ones(size(x2)); % 1xNx
    hl  = zeros(1,N+1);
    hlt = cell(1,N+1);
    for n=0:N
        f = 1-(1./kappa)*x2; % MxNx
        f = f.^(-gamma/2);   % MxNx
        if(is_broadcast_available)
            f = x.*f;                         % MxNx
        else
            f = bsxfun( @(x,y)(x.*y), x, f ); % MxNx
        end
        iotan = 2*sum(f,2)*(x1(2)-x1(1)); % Mx1
        figure(hf);
        hl(n+1) = plot( kappa, iotan, 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', MS, 'Color', cols(rem(n,NC)+1,:) );
        plot( kappa, iota(:,n+1), 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', MS, 'Color', cols(rem(n,NC)+1,:) );
        hlt{n+1} = ['n = ',num2str(n)];
        x = x.*x2;
    end
    figure(hf);
    legend(hl,hlt{:});
    %%% -------------------------------------------------------------------
end
%%% -----------------------------------------------------------------------
