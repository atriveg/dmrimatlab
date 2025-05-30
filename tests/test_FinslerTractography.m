%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_fastSweeping.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
close('all');
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ND      = 64;
R       = 200;
C       = 400;
ijk2xyz = [1,0,0;0,0.5,0;0,0,1];
theta   = 2*pi*(0:ND-1)/ND;
gi      = [ cos(theta(:)), sin(theta(:)) ];
dcosts  = ones( R, C, ND );
seeds   = zeros( R, C );
seeds(R/2,C/2) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[costs,dirs] = fastSweeping( dcosts, gi, seeds, ...
    'ijk2xyz', ijk2xyz, 'seedval', 1, 'verbose', true, ...
    'maxiters', 1000, 'change', 1.0e-12 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xit      = [ round(linspace(0,C-1,10)), zeros(1,10), round(linspace(C-1,0,10)), (C-1)*ones(1,10)  ];
yit      = [ (R-1)*ones(1,10), round(linspace(0,R-1,10)), zeros(1,10), round(linspace(R-1,0,10)) ];
zit      = ones(size(xit));
targetsi = [ yit; xit; zit ];
targetsx = ijk2xyz*targetsi;
%%%%%%%%%%%%%%
[paths,~,dirs,~] = backTracing( costs, dirs, targetsx(1:2,:), ...
    'ijk2xyz', ijk2xyz, 'step', 0.1, 'Mlength', 2*(R+C), ...
    'mcurv', 0.1, 'prune', [], 'maxthreads', 8 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = (0:C-1);
jj = (0:R-1);
close(figure(1));
figure(1);
plot_cost_and_velocity_field(ii,jj,costs,dirs,paths,ijk2xyz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ND      = 64;
R       = 200;
C       = 200;
rt      = 30;
CS      = cos((rt/360)*pi);
SS      = sin((rt/360)*pi);
ijk2xyz = [CS,-SS,0;SS,CS,0;0,0,1];
theta   = 2*pi*(0:ND-1)/ND;
gi      = [ cos(theta(:)), sin(theta(:)) ];
dcosts  = ones( R, C, ND );
% ----
fn     = gi(:,1).*gi(:,1);
fn     = fn.*(1-fn);
fn     = 1-(4-2*sqrt(2))*fn.*(1-fn);
dcosts = dcosts.*reshape(fn,[1,1,ND]);
% ----
seeds   = zeros( R, C );
seeds(R/2,C/2) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[costs,dirs] = fastSweeping( dcosts, gi, seeds, ...
    'ijk2xyz', ijk2xyz, 'seedval', 1, 'verbose', true, ...
    'maxiters', 1000, 'change', 1.0e-12 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xit      = [ round(linspace(1,C-2,10)), ones(1,10), round(linspace(C-2,1,10)), (C-2)*ones(1,10)  ];
yit      = [ (R-2)*ones(1,10), round(linspace(1,R-2,10)), ones(1,10), round(linspace(R-2,1,10)) ];
zit      = ones(size(xit));
targetsi = [ yit; xit; zit ];
targetsx = ijk2xyz*targetsi;
%%%%%%%%%%%%%%
[paths,~,dirs,~] = backTracing( costs, dirs, targetsx(1:2,:), ...
    'ijk2xyz', ijk2xyz, 'step', 0.1, 'Mlength', 2*(R+C), ...
    'mcurv', 0.1, 'prune', [], 'maxthreads', 8 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = (0:C-1);
jj = (0:R-1);
close(figure(2));
figure(2);
plot_cost_and_velocity_field(ii,jj,costs,dirs,paths,ijk2xyz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ND      = 64;
R       = 200;
C       = 200;
ijk2xyz = [1,0,0;0,1,0;0,0,1];
theta   = 2*pi*(0:ND-1)/ND;
gi      = [ cos(theta(:)), sin(theta(:)) ];
dcosts  = ones( R, C, ND );
seeds   = zeros( R, C );
seeds(1,1) = 1;
seeds(R,C) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask    = true(R,C);
NB      = 6;
h       = round(R/(NB+1));
for b=1:NB
    if(rem(b,2)==0)
        st = round(0.1*C);
        ed = C;
    else
        st = 1;
        ed = round(0.9*C);
    end
    mask(b*h-1:b*h+1,st:ed) = false;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[costs,dirs] = fastSweeping( dcosts, gi, seeds, ...
    'ijk2xyz', ijk2xyz, 'seedval', 1, 'verbose', true, ...
    'maxiters', 1000, 'change', 1.0e-12, 'mask', mask );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xit       = 90:5:110;
[xit,yit] = meshgrid(xit,xit);
xit       = xit(:)';
yit       = yit(:)';
zit       = ones(size(xit));
targetsi  = [ yit; xit; zit ];
targetsx  = ijk2xyz*targetsi;
%%%%%%%%%%%%%%
[paths,~,dirs,~] = backTracing( costs, dirs, targetsx(1:2,:), ...
    'mask', mask, ...
    'ijk2xyz', ijk2xyz, 'step', 0.1, 'Mlength', 2*(R+C), ...
    'mcurv', 0.01, 'prune', [], 'maxthreads', 8 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = (0:C-1);
jj = (0:R-1);
close(figure(3));
figure(3);
plot_cost_and_velocity_field(ii,jj,costs,dirs,paths,ijk2xyz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ND      = 64;
R       = 200;
C       = 200;
ijk2xyz = [0,-1,0;1,0,0;0,0,1];
theta   = 2*pi*(0:ND-1)/ND;
gi      = [ cos(theta(:)), sin(theta(:)) ];
dcost   = gi(:,1).*gi(:,1) + 1;
dcost   = reshape(dcost,[1,1,ND]);
[xx,yy] = meshgrid(1:C,1:R);
scost   = sqrt(xx.*xx+yy.*yy);
scost   = scost./max(scost(:)) + 0.5;
dcosts  = scost.*dcost;
seeds   = zeros( R, C );
seeds(R,C) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[costs,dirs] = fastSweeping( dcosts, gi, seeds, ...
    'ijk2xyz', ijk2xyz, 'seedval', 1, 'verbose', true, ...
    'maxiters', 1000, 'change', 1.0e-12 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xit      = [ round(linspace(0,C-1,10)), zeros(1,10) ];
yit      = [ (R-1)*ones(1,10), round(linspace(0,R-1,10)) ];
zit      = ones(size(xit));
targetsi = [ yit; xit; zit ];
targetsx = ijk2xyz*targetsi;
%%%%%%%%%%%%%%
[paths,~,dirs,~] = backTracing( costs, dirs, targetsx(1:2,:), ...
    'ijk2xyz', ijk2xyz, 'step', 0.1, 'Mlength', 2*(R+C), ...
    'mcurv', 0.1, 'prune', [], 'maxthreads', 8 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ii = (0:C-1);
jj = (0:R-1);
close(figure(4));
figure(4);
ha = plot_cost_and_velocity_field(ii,jj,costs,dirs,paths,ijk2xyz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ND      = 64;
R       = 50;
C       = R;
S       = 30;
rt      = 30;
CS      = cos((rt/360)*pi);
SS      = sin((rt/360)*pi);
ijk2xyz = [
    CS,-SS,0,0;
    SS,CS,0,0;
    0,0,1,0;
    0,0,0,1];
gi      = designGradients(ND/2);
gi      = [gi;-gi];
dcosts  = ones(R,C,S,ND);
dcost   = gi(:,1).*gi(:,1) + 0.01;
dcost   = reshape(dcost,[1,1,1,ND]);
dcosts  = dcosts.*dcost;
seeds   = zeros( R, C, S );
x = 2*(0:R-1)/R - 1;
y = 2*(0:C-1)/C - 1;
z = 2*(0:S-1)/S - 1;
[x,y,z] = meshgrid(x,y,z);
ellipse = (x.*x)/(0.5*0.5) + (y.*y)/(0.4*0.4) + (z.*z)/(0.6*0.6);
seeds(round(ellipse)==1) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[costs,dirs] = fastSweeping( dcosts, gi, seeds, ...
    'ijk2xyz', ijk2xyz, 'seedval', 1, 'verbose', true, ...
    'maxiters', 1000, 'change', 1.0e-12 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xit      = [2,C/2-1,C-3];
yit      = [2,R/2-1,R-3];
zit      = [2,S/2-1,S-3];
[xit,yit,zit] = meshgrid(xit,yit,zit);
targetsi = [ yit(:)'; xit(:)'; zit(:)'; ones(1,numel(xit)) ];
targetsx = ijk2xyz*targetsi;
%%%%%%%%%%%%%%
[paths,~,dirs,~] = backTracing( costs, dirs, targetsx(1:3,:), ...
    'ijk2xyz', ijk2xyz, 'step', 0.1, ...
    'Mlength', 3*(R+C+S), 'mcurv', 0.01, 'prune', [], 'maxthreads', 8 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(5));
figure(5);
ii = (1:R);
jj = (1:C);
kk = (1:S);
plot_cost_and_velocity_field3D(ii,jj,kk,costs,dirs,paths,ijk2xyz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ND      = 64;
R       = 40;
C       = 40;
S       = 60;
%%%%%%%%%%%%
rt      = 30;
CS      = cos((rt/360)*pi);
SS      = sin((rt/360)*pi);
ijk2xyz = [
    CS,-SS,0,0;
    SS,CS,0,0;
    0,0,1,0;
    0,0,0,1];
rt      = 15;
CS      = cos((rt/360)*pi);
SS      = sin((rt/360)*pi);
ijk2xyz = ijk2xyz*[
    CS, 0, -SS, -12;
    0,  1,   0, 21;
    SS, 0,  CS, 43;
    0,  0,   0, 1   ];
%%%%%%%%%%%%
gi      = designGradients(ND/2);
gi      = [gi;-gi];
dcosts  = ones( R, C, S, ND );
seeds   = zeros( R, C, S );
seeds(1:4,1:4,1:4) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask    = true(R,C,S);
NB      = 4;
h       = round(R/(NB+1));
for b=1:NB
    if(rem(b,2)==0)
        st = round(0.2*C);
        ed = C;
    else
        st = 1;
        ed = round(0.8*C);
    end
    mask(b*h-1:b*h+1,st:ed,:) = false;
end
%%%%%%%%%%%%%%%%%%%
hz       = round(S/(NB+1));
for b=1:NB
    if(rem(b,2)==0)
        wr = 1:hz-2;
        wc = 1:hz-2;
    else
        wr = R-hz+2:R;
        wc = C-hz+2:C;
    end
    mask(:,:,b*hz-1:b*hz+1)     = false;
    mask(wr,wc,b*hz-1:b*hz+1)   = true;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[costs,dirs] = fastSweeping( dcosts, gi, seeds, ...
    'ijk2xyz', ijk2xyz, 'seedval', 1, 'verbose', true, ...
    'maxiters', 1000, 'change', 1.0e-12, 'mask', mask );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xit       = (NB*h+h/2) + (-round(h/4):round(h/4));
yit       = C-3:C-2;
zit       = (NB*hz+hz/2+1) + (-round(hz/4):round(hz/4));
[xit,yit,zit] = meshgrid(xit,yit,zit);
xit       = xit(:)';
yit       = yit(:)';
zit       = zit(:)';
oit       = ones(size(xit));
targetsi  = [ yit; xit; zit; oit ];
targetsx  = ijk2xyz*targetsi;
%%%%%%%%%%%%%%
[paths,~,dirs,~] = backTracing( costs, dirs, targetsx(1:3,:), ...
    'mask', mask, ...
    'ijk2xyz', ijk2xyz, 'step', 0.1, 'Mlength', 3*(R+C+S), ...
    'mcurv', 0.001, 'prune', [], 'maxthreads', 8 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(6));
figure(6);
ii = (1:R);
jj = (1:C);
kk = (1:S);
plot_cost_and_velocity_field3D(ii,jj,kk,costs,dirs,paths,ijk2xyz,[h-3,h,NB*hz]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = plot_cost_and_velocity_field(ii,jj,costs,dirs,paths,ijk2xyz)
if(nargin<5)
    ijk2xyz = eye(3);
end
% -------------------------------------------------------------------------
ha(1) = subplot(1,2,1);
hold('on');
imagesc(ii,jj,costs);
axis('equal');
axis('xy');
colormap('gray');
colorbar;
% ----
for n=1:length(paths)
    streamline = paths{n};
    if(~isempty(streamline))
        NP = size(streamline,2);
        streamline = ijk2xyz\[streamline;ones(1,NP)];
        streamline = streamline(1:2,:) + 1/2; % Matlab's 1-based indexing
        plot( ha(1), streamline(2,1), streamline(1,1), ...
            'LineStyle', 'none', 'Marker', 'o', ...
            'MarkerSize', 6, 'Color', [.5,.0,.0] );
        plot( ha(1), streamline(2,1), streamline(1,1), ...
            'LineStyle', 'none', 'Marker', 'x', ...
            'MarkerSize', 6, 'Color', [.5,.0,.0] );
        plot( ha(1), streamline(2,1), streamline(1,1), ...
            'LineStyle', 'none', 'Marker', '+', ...
            'MarkerSize', 6, 'Color', [.5,.0,.0] );
        plot( ha(1), streamline(2,:), streamline(1,:), ...
            'LineStyle', '-', 'LineWidth', 2, ...
            'Marker', 'none', 'Color', [.0,.7,.0] );
    end
end
% -------------------------------------------------------------------------
ha(2) = subplot(1,2,2);
hold('on');
dcols = dirs(:,:,[2,1])/2+0.5;
dcols(:,:,3) = 0.5;
image(ii,jj,dcols);
axis('equal');
% ----
dirs(:,:,3) = 1;
[P,Q,R] = size(dirs);
dirs  = reshape( dirs, [P*Q,R] )';
dirs  = ijk2xyz\dirs - ijk2xyz\[0;0;1];
dirs  = reshape(dirs',[P,Q,R]);
dirsx = dirs(:,:,2);
dirsy = dirs(:,:,1);
NA    = 20;
dx    = round(Q/NA);
dy    = round(P/NA);
sx    = max(1,round(dx/2));
sy    = max(1,round(dy/2));
quiver( ii(sx:dx:end), jj(sy:dy:end), ...
    dirsx(sy:dy:end,sx:dx:end), dirsy(sy:dy:end,sx:dx:end), ...
    'Color', [0,0,0] );
% ----
axis('equal');
axis('xy');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ha = plot_cost_and_velocity_field3D(ii,jj,kk,costs,dirs,paths,ijk2xyz,slices)
[R,C,S] = size(costs);
if(nargin<8)
    slices = [ round(R/2), round(C/2), round(S/2) ];
end
% -----------------------------------
x = jj;
y = ii;
% ---------------
ha(1) = subplot(2,2,1);
hold('on');
costs2D = costs(:,:,slices(3));
dirs2D  = squeeze(dirs(:,:,slices(3),:));
imagesc(x,y,costs2D);
colormap('gray');
axis('equal');
axis('xy');
xlabel('j');
ylabel('i');
% ---------------
dirs2D(:,:,4) = 1;
dirs2D = reshape( dirs2D, [R*C,4] )';
dirs2D = ijk2xyz\dirs2D - ijk2xyz\[0;0;0;1];
dirs2D = reshape( dirs2D', [R,C,4] );
dirsx  = dirs2D(:,:,2);
dirsy  = dirs2D(:,:,1);
% ---------------
NA    = 20;
dx    = round(C/NA);
dy    = round(R/NA);
sx    = max(1,round(dx/2));
sy    = max(1,round(dy/2));
quiver( x(sx:dx:end), y(sy:dy:end), ...
    dirsx(sy:dy:end,sx:dx:end), dirsy(sy:dy:end,sx:dx:end), ...
    'Color', [1,0,0] );
% -----------------------------------
x = ii;
y = kk;
% ---------------
ha(2) = subplot(2,2,2);
hold('on');
costs2D = squeeze(costs(:,slices(2),:));
dirs2D  = squeeze(dirs(:,slices(2),:,:));
imagesc(x,y,costs2D');
colormap('gray');
axis('equal');
axis('xy');
xlabel('i');
ylabel('k');
% ---------------
dirs2D(:,:,4) = 1;
dirs2D = reshape( dirs2D, [R*S,4] )';
dirs2D = ijk2xyz\dirs2D - ijk2xyz\[0;0;0;1];
dirs2D = reshape( dirs2D', [R,S,4] );
dirsx  = dirs2D(:,:,1)';
dirsy  = dirs2D(:,:,3)';
% ---------------
NA    = 20;
dx    = round(S/NA);
dy    = round(R/NA);
sx    = max(1,round(dx/2));
sy    = max(1,round(dy/2));
quiver( x(sx:dx:end), y(sy:dy:end), ...
    dirsx(sy:dy:end,sx:dx:end), dirsy(sy:dy:end,sx:dx:end), ...
    'Color', [1,0,0] );
% -----------------------------------
x = jj;
y = kk;
% ---------------
ha(3) = subplot(2,2,3);
hold('on');
costs2D = squeeze(costs(slices(1),:,:));
dirs2D  = squeeze(dirs(slices(1),:,:,:));
imagesc(x,y,costs2D');
colormap('gray');
axis('equal');
axis('xy');
xlabel('j');
ylabel('k');
% ---------------
dirs2D(:,:,4) = 1;
dirs2D = reshape( dirs2D, [C*S,4] )';
dirs2D = ijk2xyz\dirs2D - ijk2xyz\[0;0;0;1];
dirs2D = reshape( dirs2D', [C,S,4] );
dirsx  = dirs2D(:,:,2)';
dirsy  = dirs2D(:,:,3)';
% ---------------
NA    = 20;
dx    = round(S/NA);
dy    = round(C/NA);
sx    = max(1,round(dx/2));
sy    = max(1,round(dy/2));
quiver( x(sx:dx:end), y(sy:dy:end), ...
    dirsx(sy:dy:end,sx:dx:end), dirsy(sy:dy:end,sx:dx:end), ...
    'Color', [1,0,0] );
% -----------------------------------
ha(4) = subplot(2,2,4);
hold('on');
axial    = squeeze(costs(:,:,slices(3)));
sagittal = squeeze(costs(:,slices(2),:));
coronal  = squeeze(costs(slices(1),:,:));
% ---
axial = ( axial - min(axial(:)) )/( max(axial(:)) - min(axial(:)) );
sagittal = ( sagittal - min(sagittal(:)) )/( max(sagittal(:)) - min(sagittal(:)) );
coronal = ( coronal - min(coronal(:)) )/( max(coronal(:)) - min(coronal(:)) );
% ---
FA = 0.9;
% ---
[JJ,II] = meshgrid(jj,ii);
KK      = ones(size(II))*slices(3);
surf( II, JJ, KK, axial, 'EdgeColor', 'none', 'FaceAlpha', FA ); colormap(gray);
% ---
[KK,II] = meshgrid(kk,ii);
JJ      = ones(size(II))*slices(2);
surf( II, JJ, KK, sagittal, 'EdgeColor', 'none', 'FaceAlpha', FA ); colormap(gray);
% ---
[KK,JJ] = meshgrid(kk,jj);
II      = ones(size(JJ))*slices(1);
surf( II, JJ, KK, coronal, 'EdgeColor', 'none', 'FaceAlpha', FA ); colormap(gray);
% -----------------------------------
for n=1:length(paths)
    streamline = paths{n};
    if(~isempty(streamline))
        NP = size(streamline,2);
        streamline = ijk2xyz\[streamline;ones(1,NP)];
        streamline = streamline(1:3,:) + 1/2; % Matlab's 1-based indexing
        plot3( ha(4), ...
            streamline(1,1), streamline(2,1), streamline(3,1), ...
            'LineStyle', 'none', 'Marker', 'o', ...
            'MarkerSize', 6, 'Color', [.5,.0,.0] );
        plot3( ha(4), ...
            streamline(1,1), streamline(2,1), streamline(3,1), ...
            'LineStyle', 'none', 'Marker', 'x', ...
            'MarkerSize', 6, 'Color', [.5,.0,.0] );
        plot3( ha(4), ...
            streamline(1,1), streamline(2,1), streamline(3,1), ...
            'LineStyle', 'none', 'Marker', '+', ...
            'MarkerSize', 6, 'Color', [.5,.0,.0] );
        plot3( ha(4), ...
            streamline(1,:), streamline(2,:), streamline(3,:), ...
            'LineStyle', '-', 'LineWidth', 2, ...
            'Marker', 'none', 'Color', [.0,.5,.0] );
    end
end
% -----------------------------------
xlabel('i');
ylabel('j');
zlabel('k');
view(30,45);
rotate3d('on');
axis('equal');
axis('xy');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
