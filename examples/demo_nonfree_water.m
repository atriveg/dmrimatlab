% -------------------------------------------------------------------------------------
clear;
close('all');
tic;
load test_data.mat;
T=toc; % This is always a large piece of data
fprintf(1,'It took %f seconds to load data\n',T);
whos -file test_data.mat; % Check the variables loaded, their kinds, and their sizes
% -------------------------------------------------------------------------------------
x   = 12:93;
y   = 10:107;
sl  = [5,7,9,11,13,15]; % These are the slices we will show
atti = atti(x,y,:,:);
mask = mask(x,y,:);
mapType = "high"; % Changing the colormap will trigger the entire demo
switch(mapType)
    case 'high'
        MAP = [0,0,1;1,1,0;1,0,1;0,1,0;0,1,1;1,0,0]; % High contrast color-map
        t0  = (0:5)'/5; t   = (0:511)'/511; % 512 interpolated colors
        MAP = [ interp1(t0,MAP(:,1),t), interp1(t0,MAP(:,2),t), interp1(t0,MAP(:,3),t) ];
    case 'gray'
        MAP = gray(512); % B/W colormap
    case 'default'
        MAP = parula(512);
end
% -------------------------------------------------------------------------------------
mu = 0.0005;%0.00015;
tic;
[lpar,lperp,f] = atti2micro( atti, gi, bi, 'tl', 1.0e-6, 'tu', 1-1.0e-6, ...
    'ADC0', 3.0e-3, 'usef', true, 'mask', mask, 'bth', 100, 'mlperp',0.01e-3, ...
    'mu', mu, 'verbose', true );
T = toc;
fprintf(1,'It took %f seconds to complete\n',T);
tic;
[f2,~] = atti2freewater( atti, gi, bi, 'L', 8, 'lambda', 0.001, 'mask', mask, 'verbose', true, 'fdepth', 3, 'fbins', 2 );
T = toc;
fprintf(1,'It took %f seconds to complete\n',T);
% -------------------------------------------------------------------------------------
F = [ f(:,:,sl(1))',f(:,:,sl(2))';f(:,:,sl(3))',f(:,:,sl(4))';f(:,:, ...
    sl(5))',f(:,:,sl(6))']; % Parallel diffusivity
LPA = [ lpar(:,:,sl(1))',lpar(:,:,sl(2))';lpar(:,:,sl(3))',lpar(:,:, ...
    sl(4))';lpar(:,:,sl(5))',lpar(:,:,sl(6))']; % Parallel diffusivity
LPP = [ lperp(:,:,sl(1))',lperp(:,:,sl(2))';lperp(:,:,sl(3))',lperp(:,:, ...
    sl(4))';lperp(:,:,sl(5))',lperp(:,:,sl(6))']; % Perpendicular diffusivity

close(figure(1));
hf = figure(1);
set(hf,'Position',[120,41,1200,650],'Name','New method based on Spherical-Means');
R=1; C=3;
r=1; c=1; subplot('Position',[(c-1)/C+0.05/C,(r-1)/R+0.05/R,0.9/C,0.9/R]);
imagesc(F,[0,1]); colormap(MAP); colorbar; title('f (non-free water)');
axis('equal');axis('off');axis('tight');
r=1; c=2; subplot('Position',[(c-1)/C+0.05/C,(r-1)/R+0.05/R,0.9/C,0.9/R]);
imagesc(LPA,[0.0,3.0e-3]); colormap(MAP); colorbar; title('\lambda_{||} (mm^2/s)');
axis('equal');axis('off');axis('tight');
r=1; c=3; subplot('Position',[(c-1)/C+0.05/C,(r-1)/R+0.05/R,0.9/C,0.9/R]);
imagesc(LPP,[0,0.3e-3]); colormap(MAP); colorbar; title('\lambda_{\perp} (mm^2/s)');
axis('equal');axis('off');axis('tight');
% -------------------------------------------------------------------------------------
F2 = [ f2(:,:,sl(1))',f2(:,:,sl(2))';f2(:,:,sl(3))',f2(:,:,sl(4))'; ...
    f2(:,:,sl(5))',f2(:,:,sl(6))']; 
close(figure(2));
hf = figure(2);
set(hf,'Position',[120,41,800,650],'Name','Comparison between methods');
R=1; C=2;
r=1; c=1; subplot('Position',[(c-1)/C+0.05/C,(r-1)/R+0.05/R,0.9/C,0.9/R]);
imagesc(1-F,[0,1]); colormap(MAP); colorbar; title('Feee water with msCAD');
axis('equal');axis('off');axis('tight');
r=1; c=2; subplot('Position',[(c-1)/C+0.05/C,(r-1)/R+0.05/R,0.9/C,0.9/R]);
imagesc(1-F2,[0,1]); colormap(MAP); colorbar; title('Free water with tensor model');
axis('equal');axis('off');axis('tight');
% -------------------------------------------------------------------------------------
