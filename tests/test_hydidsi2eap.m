function test_hydidsi2eap

sf = check_software_platform;
if(sf==2)
    pkg load statistics;
end

%%% test_hydidsi2eap.m
clc;
close('all');
tau = 40e-3;
Rth = 0.03;
lattice = [2;3;4];
[atti,gi,bi] = generate_atti_signal;
atti = reshape( atti, [1,1,1,size(atti,1)] );
tic;
[eap,dti,lattice,Qx,Qy,Qz,~,~,~] = atti2hydidsi( atti, gi, bi, ...
    'tl', 1.0e-7, 'tu', 1-1.0e-7, 'lambda', 0.4, 'tau', tau, ...
    'Rth', Rth, 'lattice', lattice, 'const', true, 'usemex', true, ...
    'miters', 10000, 'otol', 1.0e-6, 'stol', 1.0e-6, 'ctol', 1.0e-8 );
T = toc;
fprintf(1,'It took %1.5f seconds to fit the EAP\n',T);
eap  = reshape( eap, [1,numel(eap)] );
dti  = reshape( dti,  [1,6] );
% Note: spline interpolation does not seem to work with
%       interp3 in GNU Octave
intp = 'linear'; % 'nearest' / 'linear' / 'spline'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ogrid = zeros( 2*lattice(2)+1, 2*lattice(1)+1, 2*lattice(3)+1 );
NPG   = numel(ogrid);
ogrid((NPG+1)/2:end)  = eap;
ogrid((NPG-1)/2:-1:1) = eap(2:end);
egrid = zeros( 2*lattice(2)+5, 2*lattice(1)+5, 2*lattice(3)+5 );
egrid(3:end-2,3:end-2,3:end-2) = ogrid;
xgl = (-lattice(1)-2:lattice(1)+2)/Qx;
ygl = (-lattice(2)-2:lattice(2)+2)/Qy;
zgl = (-lattice(3)-2:lattice(3)+2)/Qz;
[xg,yg,zg] = meshgrid(xgl,ygl,zgl);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u1,u2,u3] = dti2spectrum( reshape(dti,[1,1,1,6]) ); % In descending order
u1 = u1(:);
u2 = u2(:);
u3 = u3(:);
U  = [u3,u2,u1]; % We need them in ascending order!!
U(:,1) = cross(U(:,2),U(:,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RG   = 80;
NG   = 2*RG+1;
rx   = 0.04*(-RG:RG)/RG;
ry   = 0.04*(-RG:RG)/RG;
rz   = 0.04*(-RG:RG)/RG;
hx   = rx(2)-rx(1);
hy   = ry(2)-ry(1);
hz   = rz(2)-rz(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rxg,ryg,rzg] = meshgrid(rx,ry,rz);
szo  = size(rxg);
rxyz = (U')*[ rxg(:)'; ryg(:)'; rzg(:)' ];
rxg  = reshape( rxyz(1,:), szo );
ryg  = reshape( rxyz(2,:), szo );
rzg  = reshape( rxyz(3,:), szo );
tic;
eap_th = interp3(xg,yg,zg,egrid,rxg,ryg,rzg,intp(1),0);
%eap_th(eap_th<0) = 0;
T = toc;
fprintf(1,'It took %1.5f seconds to interpolate the EAP\n',T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rxg,ryg,rzg] = meshgrid(rx,ry,rz);
ri  = sqrt(rxg.*rxg+ryg.*ryg+rzg.*rzg);
ri  = ri(:)';
iri = ri;
iri(abs(ri)>eps) = 1./ri(abs(ri)>eps);
iri(abs(ri)<eps) = 0;
ui = [rxg(:)';ryg(:)';rzg(:)'].*iri;
tic;
eap = hydidsi2eap( ...
    reshape(eap,[1,1,1,numel(eap)]), ...
    reshape(dti,[1,1,1,6]), ...
    lattice, Qx, Qy, Qz, ui', ri', 'interp', intp(1) );
eap = reshape(eap,[NG,NG,NG]);
T = toc;
fprintf(1,'It took %1.5f seconds to render the EAP\n',T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_X = sum(sum(eap_th,1),3)*hy*hz;
f_Y = sum(sum(eap_th,2),3)*hx*hz;
f_Z = sum(sum(eap_th,1),2)*hx*hy;

f_X = reshape(f_X,[1,numel(f_X)]);
f_Y = reshape(f_Y,[1,numel(f_Y)]);
f_Z = reshape(f_Z,[1,numel(f_Z)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_X2 = sum(sum(eap,1),3)*hy*hz;
f_Y2 = sum(sum(eap,2),3)*hx*hz;
f_Z2 = sum(sum(eap,1),2)*hx*hy;

f_X2 = reshape(f_X2,[1,numel(f_X2)]);
f_Y2 = reshape(f_Y2,[1,numel(f_Y2)]);
f_Z2 = reshape(f_Z2,[1,numel(f_Z2)]);
% -------------------------------------------------------------------------
close(figure(3003));
hf = figure(3003);
hp = [];
% -------------
figure(hf);
ha = subplot(1,3,1); hold('on'); grid('on');
hp(1) = plot(rx,f_X, 'LineStyle',':','LineWidth',2,'Color',[.0,.5,.0]);
hp(2) = plot(rx,f_X2,'LineStyle','-.','LineWidth',2,'Color',[.5,.0,.0]);
legend(hp,'Interpolated','Computed','Location','NorthWest','FontSize',16); legend('boxoff');
xlabel('x [mm]','FontSize',18);
ylabel('f_X(x) [mm^{-1}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinate X','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,2); hold('on'); grid('on');
hp(1) = plot(ry,f_Y, 'LineStyle',':','LineWidth',2,'Color',[.0,.5,.0]);
hp(2) = plot(ry,f_Y2,'LineStyle','-.','LineWidth',2,'Color',[.5,.0,.0]);
legend(hp,'Interpolated','Computed','Location','NorthWest','FontSize',16); legend('boxoff');
xlabel('y [mm]','FontSize',18);
ylabel('f_Y(y) [mm^{-1}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinate Y','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,3); hold('on'); grid('on');
hp(1) = plot(rz,f_Z, 'LineStyle',':','LineWidth',2,'Color',[.0,.5,.0]);
hp(2) = plot(rz,f_Z2,'LineStyle','-.','LineWidth',2,'Color',[.5,.0,.0]);
legend(hp,'Interpolated','Computed','Location','NorthWest','FontSize',16); legend('boxoff');
xlabel('z [mm]','FontSize',18);
ylabel('f_Z(z) [mm^{-1}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinate Z','FontSize',18);
drawnow;
% -------------
% -------------------------------------------------------------------------

f_XY = sum(eap_th,3)*hz;
f_XZ = sum(eap_th,1)*hy; f_XZ = reshape(f_XZ,[size(f_XZ,2),size(f_XZ,3)])';
f_YZ = sum(eap_th,2)*hx; f_YZ = reshape(f_YZ,[size(f_YZ,1),size(f_YZ,3)])';

f_XY2 = sum(eap,3)*hz;
f_XZ2 = sum(eap,1)*hy; f_XZ2 = reshape(f_XZ2,[size(f_XZ2,2),size(f_XZ2,3)])';
f_YZ2 = sum(eap,2)*hx; f_YZ2 = reshape(f_YZ2,[size(f_YZ2,1),size(f_YZ2,3)])';
% -------------------------------------------------------------------------
close(figure(4004));
hf = figure(4004);
% -------------
figure(hf);
ha = subplot(1,3,1); hold('on'); grid('on'); rotate3d('on');
% --
surf(rx,ry,f_XY,ones(size(f_XY)),'EdgeColor',[0,0,0],'EdgeAlpha',0.2,'FaceAlpha',0.5,'FaceColor',[.0,.5,.0]);
surf(rx,ry,f_XY2,ones(size(f_XY2)),'EdgeColor',[0,0,0],'EdgeAlpha',0.2,'FaceAlpha',0.5,'FaceColor',[.5,.0,.0]);
% --
xlabel('x [mm]','FontSize',18);
ylabel('y [mm]','FontSize',18);
zlabel('f_{XY}(x,y) [mm^{-2}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinates X and Y [rotate]','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,2); hold('on'); grid('on'); rotate3d('on');
% --
surf(rx,rz,f_XZ,ones(size(f_XZ)),'EdgeColor',[0,0,0],'EdgeAlpha',0.2,'FaceAlpha',0.5,'FaceColor',[.0,.5,.0]);
surf(rx,rz,f_XZ2,ones(size(f_XZ2)),'EdgeColor',[0,0,0],'EdgeAlpha',0.2,'FaceAlpha',0.5,'FaceColor',[.5,.0,.0]);
% --
xlabel('x [mm]','FontSize',18);
ylabel('z [mm]','FontSize',18);
zlabel('f_{XZ}(x,z) [mm^{-2}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinates X and Z [rotate]','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,3); hold('on'); grid('on'); rotate3d('on');
% --
surf(ry,rz,f_YZ,ones(size(f_YZ)),'EdgeColor',[0,0,0],'EdgeAlpha',0.2,'FaceAlpha',0.5,'FaceColor',[.0,.5,.0]);
surf(ry,rz,f_YZ2,ones(size(f_YZ2)),'EdgeColor',[0,0,0],'EdgeAlpha',0.2,'FaceAlpha',0.5,'FaceColor',[.5,.0,.0]);
% --
xlabel('y [mm]','FontSize',18);
ylabel('z [mm]','FontSize',18);
zlabel('f_{YZ}(y,z) [mm^{-2}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinates Y and Z [rotate]','FontSize',18);
drawnow;
% -------------------------------------------------------------------------
end

function [atti,gi,bi] = generate_atti_signal
[gi,bi] = generate_ms_sampling;
Diso = 3.0e-3;
f1 = 0.6;
% -------------------------------------
% Fiber 1
dir   = [1.1,0.2,0.1];
dir   = dir/norm(dir);
lpar  = 2.0e-3;
lperp = 0.3e-3;
ffw   = 0.1;
fin   = 0.4;
uTv   = sum(gi.*dir,2);
uTv   = uTv.*uTv;
% ---
cfw   = exp( -bi*Diso );
cin   = exp( -bi.*(lpar*uTv) );
cen   = exp( -bi.*((lpar-lperp)*uTv+lperp) );
% ---
S1    = ffw*cfw + fin*cin + (1-ffw-fin)*cen;
% -------------------------------------
% Fiber 2
dir   = [0.04,0.06,0.9];
dir   = dir/norm(dir);
lpar  = 2.1e-3;
lperp = 0.15e-3;
ffw   = 0.075;
fin   = 0.6;
uTv   = sum(gi.*dir,2);
uTv   = uTv.*uTv;
% ---
cfw   = exp( -bi*Diso );
cin   = exp( -bi.*(lpar*uTv) );
cen   = exp( -bi.*((lpar-lperp)*uTv+lperp) );
% ---
S2    = ffw*cfw + fin*cin + (1-ffw-fin)*cen;
% -------------------------------------
atti = f1*S1 + (1-f1)*S2;
end

function [gi,bi] = generate_ms_sampling
bs = [500,1000,1500,2000,3000,5000];
Ns = [15,30,30,40,50,60];
bi = [];
gi = [];
for k=1:length(bs)
    bi = [ bi; bs(k)*ones(Ns(k),1) ];    %#ok<AGROW>
    gi = [ gi; designGradients(Ns(k)) ]; %#ok<AGROW>
end
end

