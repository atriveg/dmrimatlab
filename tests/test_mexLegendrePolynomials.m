%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_mexLegendrePolynomials.m
clear;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L     = 16;
theta = [ 0, 1.0e-15, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3, 0.01, 0.1, pi/16, pi/8, pi/4, 1, pi/2 ];
phi   = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,p] = meshgrid(theta,phi);
t     = t(:);
p     = p(:);
x     = sin(t).*cos(p);
y     = sin(t).*sin(p);
z     = cos(t);
gi    = [x,y,z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1     = legendre(L,z');
P1     = P1'; % G x (L+1)
[~,P2] = mexGenerateSHMatrix(L,gi); % G x (L+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for m=1:L+1
    close(figure(1));
    hf1 = figure(1);
    hold('on');
    grid('on');
    plot( 1:length(z), P1(:,m), 'LineStyle', 'none', 'Color', [.0,.5,.0], ...
        'Marker', 'o', 'MarkerSize', 14 );
    plot( 1:length(z), P2(:,m), 'LineStyle', 'none', 'Color', [.5,.0,.0], ...
        'Marker', 'x', 'MarkerSize', 14 );
    title(sprintf('L=%i, m=%i',L,m-1));
    ha1 = get(hf1,'CurrentAxes');
    axs = axis(ha1);
    for k=1:length(z)
        text( k, P1(k,m), sprintf('   z = %1.6f',z(k)), ...
            'Rotation', 90, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', 'FontSize', 16, ...
            'FontAngle', 'italic' );
    end
    axs(4) = axs(4) + (axs(4)-axs(3))/2;
    axs(1) = axs(1) - 1;
    axs(2) = axs(2) + 1;
    axis(ha1,axs);
    set(hf1,'Position',[400 300 600 500]);
    set(ha1,'FontSize',16);
    drawnow;
    pause;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
