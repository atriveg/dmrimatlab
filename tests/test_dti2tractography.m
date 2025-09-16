function test_dti2tractography
% -------------------------------------------------------------------------
clear;
close('all');
clc;
rng(131210);
% -------------------------------------------------------------------------
noise   = 20e-5;
NI      = 110;
NJ      = 100;
NK      = 70;
ijk2xyz = eye(4);
[u1,u2,u3,II,JJ] = createEigVecs(NI,NJ,NK);
% -------------------------------------------------------------------------
% -----------
u1x = squeeze(u1(:,:,round(NK*0.9),1));
u1y = squeeze(u1(:,:,round(NK*0.9),2));
% -----------
u2x = squeeze(u2(:,:,round(NK*0.9),1));
u2y = squeeze(u2(:,:,round(NK*0.9),2));
% -----------
u3x = squeeze(u3(:,:,round(NK*0.9),1));
u3y = squeeze(u3(:,:,round(NK*0.9),2));
% -----------
close(figure(1));
figure(1);
hold('on');
d = 8;
quiver( JJ(1:d:end,1:d:end), II(1:d:end,1:d:end), u1y(1:d:end,1:d:end), u1x(1:d:end,1:d:end), 'Color', [.0,.0,.5] );
quiver( JJ(1:d:end,1:d:end), II(1:d:end,1:d:end), u2y(1:d:end,1:d:end), u2x(1:d:end,1:d:end), 'Color', [.0,.5,.0] );
quiver( JJ(1:d:end,1:d:end), II(1:d:end,1:d:end), u3y(1:d:end,1:d:end), u3x(1:d:end,1:d:end), 'Color', [.5,.0,.0] );
axis('equal');
axis('xy');
title('Eigenvectors map');
xlabel('j');
ylabel('i');
% -------------------------------------------------------------------------
l1 = 1.8e-3 * ones(NI,NJ,NK);
l2 = 0.7e-3 * ones(NI,NJ,NK);
l3 = 0.3e-3 * ones(NI,NJ,NK);
[l1,l2,l3] = borderFA(l1,l2,l3);
dti = gatherTensor(u1,u2,u3,l1,l2,l3,noise);
% -------------------------------------------------------------------------
[u1c,~,~,l1c,l2c,l3c] = dti2spectrum( dti );
FA  = spectrum2scalar( l1c, l2c, l3c, 'scalar', 'fa' );
RGB = spectrum2colorcode( u1c, l1c, l2c, l3c );
% -------------------------------------------------------------------------
close(figure(2));
figure(2);
subplot(1,2,1);
imshow(FA(:,:,round(NK*0.9)),[0,1]); colormap; colorbar; axis('xy');
xlabel('j');
ylabel('i');
subplot(1,2,2);
imshow(squeeze(RGB(:,:,round(NK*0.9),:))); axis('xy');
xlabel('j');
ylabel('i');
% -------------------------------------------------------------------------
rad    = min(2*NI/5,2*NJ/5);
NS     = 5;
seedsi = zeros(3,NS);
for n=1:NS
    seedsi(1,n) = round( NI/2 + rad*cos(2*pi*(n-1)/NS) );
    seedsi(2,n) = round( NJ/2 + rad*sin(2*pi*(n-1)/NS) );
    seedsi(3,n) = round( NK/2 );
end
run_single_test(3,dti,ijk2xyz,seedsi,FA,u1,II,JJ);
% -------------------------------------------------------------------------
l1 = 0.7e-3 * ones(NI,NJ,NK);
l2 = 0.3e-3 * ones(NI,NJ,NK);
l3 = 1.8e-3 * ones(NI,NJ,NK);
[l1,l2,l3] = borderFA(l1,l2,l3);
dti = gatherTensor(u1,u2,u3,l1,l2,l3,noise);
run_single_test(4,dti,ijk2xyz,seedsi,FA,u3,II,JJ);
% -------------------------------------------------------------------------
l1 = 0.7e-3 * ones(NI,NJ,NK);
l2 = 1.8e-3 * ones(NI,NJ,NK);
l3 = 0.3e-3 * ones(NI,NJ,NK);
[l1,l2,l3] = borderFA(l1,l2,l3);
dti = gatherTensor(u1,u2,u3,l1,l2,l3,noise);
% -----------------
rad = min(NI/4,NJ/4);
NS  = 7;
seedsi = zeros(3,NS);
for n=1:NS
    seedsi(1,n) = round( NI/2 + rad*cos(2*pi*(n-1)/NS) );
    seedsi(2,n) = round( NJ/2 + rad*sin(2*pi*(n-1)/NS) );
    seedsi(3,n) = round( NK/2 );
end
run_single_test(5,dti,ijk2xyz,seedsi,FA,u2,II,JJ);
% -------------------------------------------------------------------------
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function run_single_test(figid,dti,ijk2xyz,seedsi,FA,uMAX,II,JJ)
seedsx       = ijk2xyz*[seedsi;ones(1,size(seedsi,2))];
seedsx       = seedsx(1:3,:);
[NI,NJ,NK,~] = size(dti);
% ------------
[paths,~,stopconds] = dti2tractography( dti, seedsx, 'scalar', 'fa', ...
    'threshold', 0.4, 'ijk2xyz', ijk2xyz, 'step', 0.1, ...
    'Mlength', (NI+NJ+NK)*4, 'mcurv', 0.5, 'prune', [], 'maxthreads', 1 );
% ------------
close(figure(figid));
figure(figid);
hold('on');
H = 10;
for s=0:H:NK-1
    surf( II, JJ, s*ones(size(II)), FA(:,:,s+1), 'EdgeColor', 'none', ...
        'FaceAlpha', 0.1 );
end
% ------------
for n = 1:size(seedsi,2)
    v1 = squeeze( uMAX( seedsi(1,n)+1, seedsi(2,n)+1, seedsi(3,n)+1, : ) );
    dd = [ seedsi(:,n)+1, seedsi(:,n)+1+v1*10 ];
    plot3( dd(2,:), dd(1,:), dd(3,:), 'LineWidth', 2, ...
        'Color', [.0,.0,.0], 'LineStyle', '-', 'Marker', '.', ...
        'MarkerSize', 18 );
    plot3( seedsi(2,n)+1, seedsi(1,n)+1, seedsi(3,n)+1, 'Marker', '.', ...
        'MarkerSize', 32, 'Color', [.0,.0,.0], 'LineStyle', 'none' );
end
axis('equal');
xlabel('j');
ylabel('i');
zlabel('k');
view(45,30);
rotate3d('on');
% ------------
for n = 1:length(paths)/2
    % --------------
    path = paths{2*n-1};
    path = ijk2xyz\[path;ones(1,size(path,2))] + 1;
    plot3( path(2,:), path(1,:), path(3,:), 'LineWidth', 2, ...
        'LineStyle', '-', 'Color', [.0,.5,.5] );
    if(stopconds(2*n-1)==1)
        plot3( path(2,end), path(1,end), path(3,end), ...
            'LineStyle', 'none', 'Color', [.0,.8,.0], ...
            'Marker', '.', 'MarkerSize', 15 );
    else
        plot3( path(2,end), path(1,end), path(3,end), ...
            'LineStyle', 'none', 'Color', [.5,.0,.0], ...
            'Marker', 'o', 'MarkerSize', 15 );
        plot3( path(2,end), path(1,end), path(3,end), ...
            'LineStyle', 'none', 'Color', [.5,.0,.0], ...
            'Marker', 'x', 'MarkerSize', 15 );
    end
    % --------------
    path = paths{2*n};
    path = ijk2xyz\[path;ones(1,size(path,2))] + 1;
    plot3( path(2,:), path(1,:), path(3,:), 'LineWidth', 2, ...
        'LineStyle', '-', 'Color', [.5,.5,.0] );
    if(stopconds(2*n)==1)
        plot3( path(2,end), path(1,end), path(3,end), ...
            'LineStyle', 'none', 'Color', [.0,.8,.0], ...
            'Marker', '.', 'MarkerSize', 15 );
    else
        plot3( path(2,end), path(1,end), path(3,end), ...
            'LineStyle', 'none', 'Color', [.5,.0,.0], ...
            'Marker', 'o', 'MarkerSize', 15 );
        plot3( path(2,end), path(1,end), path(3,end), ...
            'LineStyle', 'none', 'Color', [.5,.0,.0], ...
            'Marker', 'x', 'MarkerSize', 15 );
    end
    % --------------
end
% ------------
end
% -------------------------------------------------------------------------
function [u1,u2,u3,II,JJ] = createEigVecs(NI,NJ,NK)

ii = (0:NI-1);
jj = (0:NJ-1);
[JJ,II] = meshgrid(jj,ii);

phi = atan2(JJ-NJ/2,II-NI/2);

h   = linspace(0.05,0.995,NK);
h   = permute(h,[1,3,2]);
u1x = -sin(phi).*h;
u1y = cos(phi).*h;
u1z = ones(size(u1x)).*sqrt(1-h.*h);
u1(:,:,:,1) = u1x;
u1(:,:,:,2) = u1y;
u1(:,:,:,3) = u1z;

u2x = cos(phi);
u2y = sin(phi);
u2z = zeros(size(u2x));
u2(:,:,1) = u2x;
u2(:,:,2) = u2y;
u2(:,:,3) = u2z;
u2 = reshape(u2,[NI,NJ,1,3]);
u2 = repmat(u2,[1,1,NK,1]);

u3 = cross(u1,u2,4);

end
% -------------------------------------------------------------------------
function dti = gatherTensor(u1,u2,u3,l1,l2,l3,noise)
D11 = l1.*u1(:,:,:,1).*u1(:,:,:,1) + l2.*u2(:,:,:,1).*u2(:,:,:,1) + l3.*u3(:,:,:,1).*u3(:,:,:,1);
D12 = l1.*u1(:,:,:,1).*u1(:,:,:,2) + l2.*u2(:,:,:,1).*u2(:,:,:,2) + l3.*u3(:,:,:,1).*u3(:,:,:,2);
D13 = l1.*u1(:,:,:,1).*u1(:,:,:,3) + l2.*u2(:,:,:,1).*u2(:,:,:,3) + l3.*u3(:,:,:,1).*u3(:,:,:,3);
D22 = l1.*u1(:,:,:,2).*u1(:,:,:,2) + l2.*u2(:,:,:,2).*u2(:,:,:,2) + l3.*u3(:,:,:,2).*u3(:,:,:,2);
D23 = l1.*u1(:,:,:,2).*u1(:,:,:,3) + l2.*u2(:,:,:,2).*u2(:,:,:,3) + l3.*u3(:,:,:,2).*u3(:,:,:,3);
D33 = l1.*u1(:,:,:,3).*u1(:,:,:,3) + l2.*u2(:,:,:,3).*u2(:,:,:,3) + l3.*u3(:,:,:,3).*u3(:,:,:,3);
dti = zeros( [size(l1),6] );
dti(:,:,:,1) = D11;
dti(:,:,:,2) = D12;
dti(:,:,:,3) = D13;
dti(:,:,:,4) = D22;
dti(:,:,:,5) = D23;
dti(:,:,:,6) = D33;
if(nargin>6)
    dti = dti + noise*randn(size(dti));
end
end
% -------------------------------------------------------------------------
function [l1,l2,l3] = borderFA(l1,l2,l3)
[NI,NJ,NK] = size(l1);
% -----------
l1(1:2,:,:) = 1.0e-3; l1(NI-1:NI,:,:) = 1.0e-3;
l2(1:2,:,:) = 1.0e-3; l2(NI-1:NI,:,:) = 1.0e-3;
l3(1:2,:,:) = 1.0e-3; l3(NI-1:NI,:,:) = 1.0e-3;
% -----------
l1(:,1:2,:) = 1.0e-3; l1(:,NJ-1:NJ,:) = 1.0e-3;
l2(:,1:2,:) = 1.0e-3; l2(:,NJ-1:NJ,:) = 1.0e-3;
l3(:,1:2,:) = 1.0e-3; l3(:,NJ-1:NJ,:) = 1.0e-3;
% -----------
l1(:,:,1:2) = 1.0e-3; l1(:,:,NK-1:NK) = 1.0e-3;
l2(:,:,1:2) = 1.0e-3; l2(:,:,NK-1:NK) = 1.0e-3;
l3(:,:,1:2) = 1.0e-3; l3(:,:,NK-1:NK) = 1.0e-3;
end
% -------------------------------------------------------------------------

