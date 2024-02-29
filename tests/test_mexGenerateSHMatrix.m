%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_mexGenerateSHMatrix.m
clear;
close all;
format compact;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L     = 16;
K     = (L+1)*(L+2)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gi = spiralPhyllotaxis(K);
gi = optimizeGradients(gi,'verbose',true);
% ---------
% [~,theta,phi] = HARDIG( gi );
% B1 = GenerateSHMatrix(L,[theta,phi]);
% B2 = mexGenerateSHMatrix(L,gi,theta,phi);
% ---------
B1 = GenerateSHMatrix(L,gi);
B2 = mexGenerateSHMatrix(L,gi);
% ---------
close(figure(2));
hf2 = figure(2);
subplot(1,2,1);
imshow([B1,B2],[]);
colormap(jet);
colorbar;
subplot(1,2,2);
imshow( 2*(B1-B2)./(abs(B1+B2)+100*eps), [] );
colormap(jet);
colorbar;
pause;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = [ 0, 1.0e-15, 1.0e-12, 1.0e-9, 1.0e-6, 1.0e-3, 0.01, 0.1, pi/16, pi/8, pi/4, 1, pi/2 ];
phi   = [ 0, pi/4, pi/2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,p] = meshgrid(theta,phi);
t     = t(:);
p     = p(:);
x     = sin(t).*cos(p);
y     = sin(t).*sin(p);
z     = cos(t);
gi    = [x,y,z];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B1 = GenerateSHMatrix(L,gi);
B2 = mexGenerateSHMatrix(L,gi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1;
for l=0:2:L
    for m=-l:1:l
        close(figure(1));
        hf1 = figure(1);
        hold('on');
        grid('on');
        zoom('on');
        plot( 1:length(z), B1(:,k), 'LineStyle', 'none', 'Color', [.0,.5,.0], ...
            'Marker', 'o', 'MarkerSize', 14 );
        plot( 1:length(z), B2(:,k), 'LineStyle', 'none', 'Color', [.5,.0,.0], ...
            'Marker', 'x', 'MarkerSize', 14 );
        title(sprintf('l=%i, m=%i',l,m));
        ha1 = get(hf1,'CurrentAxes');
        axs = axis(ha1);
        for n=1:length(z)
            text( n, B1(n,k), sprintf('   \\theta = %1.2f, \\phi=%1.2f',t(n),p(n)), ...
                'Rotation', 90, 'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle', 'FontSize', 14, ...
                'FontAngle', 'italic', 'Interpreter', 'tex' );
        end
        axs(4) = axs(4) + (axs(4)-axs(3))/2;
        axs(1) = axs(1) - 1;
        axs(2) = axs(2) + 1;
        axis(ha1,axs);
        set(hf1,'Position',[200 300 1400 500]);
        set(ha1,'FontSize',16);
        drawnow;
        pause;
        k=k+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
