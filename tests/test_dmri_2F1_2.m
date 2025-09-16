% test_dmri_2F1__
clear;

sf = check_software_platform;
if(sf==2)
    pkg load gsl;
end

N = 1000;
% ------------
g = -1.13;
x = linspace(-150,150,N)/151;
tic;
if(sf==1)
    f11 = real(hypergeom( [1/2,g/2], 3/2, x ));
else
    f11 = real(gsl_sf_hyperg_2F1( 1/2, g/2, 3/2, x ));
end
T = toc;
fprintf('[gamma=%1.4f] Computation time with native''s: %1.4f seconds\n',g,T);
tic;
f12 = dmri_2F1_( g, x );
T = toc;
fprintf('[gamma=%1.4f] Computation time with mex implementation: %1.4f seconds\n',g,T);
close(figure(1));
h1 = figure(1);
hold('on');
grid('on');
hp = [];
hp(1) = plot( x, f11, 'LineStyle', '--', 'LineWidth', 2, 'Color', [.0,.5,.0] );
hp(2) = plot( x, f12, 'LineStyle', ':', 'LineWidth', 2, 'Color', [.5,.0,.0] );
xlabel('x','FontSize',20);
ylabel(sprintf('_2F_1(  [1/2,%1.3f/2],  3/2,  x  )',g),'FontSize',20);
hl = legend(hp,'GSL''s','Mex implemented');
set(hl,'FontSize',18);
% ------------
g = 1.19;
x = linspace(-150,150,N)/151;
tic;
if(sf==1)
    f21 = real(hypergeom( [1/2,g/2], 3/2, x ));
else
    f21 = real(gsl_sf_hyperg_2F1( 1/2, g/2, 3/2, x ));
end
T = toc;
fprintf('[gamma=%1.4f] Computation time with native''s: %1.4f seconds\n',g,T);
tic;
f22 = dmri_2F1_( g, x );
T = toc;
fprintf('[gamma=%1.4f] Computation time with mex implementation: %1.4f seconds\n',g,T);
close(figure(2));
h2 = figure(2);
hold('on');
grid('on');
hp = [];
hp(1) = plot( x, f21, 'LineStyle', '--', 'LineWidth', 2, 'Color', [.0,.5,.0] );
hp(2) = plot( x, f22, 'LineStyle', ':', 'LineWidth', 2, 'Color', [.5,.0,.0] );
xlabel('x','FontSize',20);
ylabel(sprintf('_2F_1(  [1/2,%1.3f/2],  3/2,  x  )',g),'FontSize',20);
hl = legend(hp,'GSL''s','Mex implemented');
set(hl,'FontSize',18);
% ------------
g = 1.0;
x = linspace(-150,150,N)/151;
tic;
if(sf==1)
    f31 = real(hypergeom( [1/2,g/2], 3/2, x ));
else
    f31 = real(gsl_sf_hyperg_2F1( 1/2, g/2, 3/2, x ));
end
T = toc;
fprintf('[gamma=%1.4f] Computation time with native''s: %1.4f seconds\n',g,T);
tic;
f32 = dmri_2F1_( g, x );
T = toc;
fprintf('[gamma=%1.4f] Computation time with mex implementation: %1.4f seconds\n',g,T);
close(figure(3));
h3 = figure(3);
hold('on');
grid('on');
hp = [];
hp(1) = plot( x, f31, 'LineStyle', '--', 'LineWidth', 2, 'Color', [.0,.5,.0] );
hp(2) = plot( x, f32, 'LineStyle', ':', 'LineWidth', 2, 'Color', [.5,.0,.0] );
xlabel('x','FontSize',20);
ylabel(sprintf('_2F_1(  [1/2,%1.3f/2],  3/2,  x  )',g),'FontSize',20);
hl = legend(hp,'GSL''s','Mex implemented');
set(hl,'FontSize',18);
% ------------
g = 1.0+sqrt(eps)/2;
x = linspace(-150,150,N)/151;
tic;
if(sf==1)
    f41 = real(hypergeom( [1/2,g/2], 3/2, x ));
else
    f41 = real(gsl_sf_hyperg_2F1( 1/2, g/2, 3/2, x ));
end
T = toc;
fprintf('[gamma=%1.4f] Computation time with native''s: %1.4f seconds\n',g,T);
tic;
f42 = dmri_2F1_( g, x );
T = toc;
fprintf('[gamma=%1.4f] Computation time with mex implementation: %1.4f seconds\n',g,T);
close(figure(4));
h4 = figure(4);
hold('on');
grid('on');
hp = [];
hp(1) = plot( x, f41, 'LineStyle', '--', 'LineWidth', 2, 'Color', [.0,.5,.0] );
hp(2) = plot( x, f42, 'LineStyle', ':', 'LineWidth', 2, 'Color', [.5,.0,.0] );
xlabel('x','FontSize',20);
ylabel(sprintf('_2F_1(  [1/2,%1.3f/2],  3/2,  x  )',g),'FontSize',20);
hl = legend(hp,'GSL''s','Mex implemented');
set(hl,'FontSize',18);
% ------------
g = 0.67;
x = linspace(-150,150,N)/151;
tic;
if(sf==1)
    f51 = real(hypergeom( [1/2,g/2], 3/2, x ));
else
    f51 = real(gsl_sf_hyperg_2F1( 1/2, g/2, 3/2, x ));
end
T = toc;
fprintf('[gamma=%1.4f] Computation time with native''s: %1.4f seconds\n',g,T);
tic;
f52 = dmri_2F1_( g, x );
T = toc;
fprintf('[gamma=%1.4f] Computation time with mex implementation: %1.4f seconds\n',g,T);
close(figure(5));
h5 = figure(5);
hold('on');
grid('on');
hp = [];
hp(1) = plot( x, f51, 'LineStyle', '--', 'LineWidth', 2, 'Color', [.0,.5,.0] );
hp(2) = plot( x, f52, 'LineStyle', ':', 'LineWidth', 2, 'Color', [.5,.0,.0] );
xlabel('x','FontSize',20);
ylabel(sprintf('_2F_1(  [1/2,%1.3f/2],  3/2,  x  )',g),'FontSize',20);
hl = legend(hp,'GSL''s','Mex implemented');
set(hl,'FontSize',18);
