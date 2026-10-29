function test_hydidsi2samples

sf = check_software_platform;
if(sf==2)
    pkg load statistics;
end

%%% test_hydidsi2samples.m
clc;
close('all');
N  = 1e6;
n_th = 1;
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
intp = 'linear'; % 'nearest' / 'linear' / 'spline'
c    = [];
seed = [13,7225,190318];
tic;
[rsx,rsy,rsz,c] = hydidsi2samples( eap, dti, lattice', ...
    [Qx(:),Qy(:),Qz(:)], N, intp(1), c, seed, n_th );
T = toc;
fprintf(1,'It took %1.4f seconds to generate %d samples\n',T,N);
fprintf(1,'The optimal value found for c is %1.4f\n',c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lx   = max(abs(rsx));
Ly   = max(abs(rsy));
Lz   = max(abs(rsz));
RG   = 40;
NG   = 2*RG+1;
rx   = Lx*(-RG:RG)/RG;
ry   = Ly*(-RG:RG)/RG;
rz   = Lz*(-RG:RG)/RG;
hx   = rx(2)-rx(1);
hy   = ry(2)-ry(1);
hz   = rz(2)-rz(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
[rxg,ryg,rzg] = meshgrid(rx,ry,rz);
ri  = sqrt(rxg.*rxg+ryg.*ryg+rzg.*rzg);
ri  = ri(:)';
iri = ri;
iri(abs(ri)>eps) = 1./ri(abs(ri)>eps);
iri(abs(ri)<eps) = 0;
ui = [rxg(:)';ryg(:)';rzg(:)'].*iri;
tic;
eap_th = hydidsi2eap( ...
    reshape(eap,[1,1,1,numel(eap)]), ...
    reshape(dti,[1,1,1,6]), ...
    lattice, Qx, Qy, Qz, ui', ri', 'interp', intp(1) );
eap_th = reshape(eap_th,[NG,NG,NG]);
T = toc;
fprintf(1,'It took %1.5f seconds to render the EAP\n',T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_X = sum(sum(eap_th,1),3)*hy*hz;
f_Y = sum(sum(eap_th,2),3)*hx*hz;
f_Z = sum(sum(eap_th,1),2)*hx*hy;

f_X = reshape(f_X,[1,numel(f_X)]);
f_Y = reshape(f_Y,[1,numel(f_Y)]);
f_Z = reshape(f_Z,[1,numel(f_Z)]);

nnan = sum(isnan(rsx));
ninf = sum(isinf(rsx));
mmin = min(rsx);
mmax = max(rsx);
fprintf(1,'Generated x coordinates range [%1.3e,%1.3e], with %d infs and %d nans\n',mmin,mmax,ninf,nnan);
nnan = sum(isnan(rsy));
ninf = sum(isinf(rsy));
mmin = min(rsy);
mmax = max(rsy);
fprintf(1,'Generated y coordinates range [%1.3e,%1.3e], with %d infs and %d nans\n',mmin,mmax,ninf,nnan);
nnan = sum(isnan(rsz));
ninf = sum(isinf(rsz));
mmin = min(rsz);
mmax = max(rsz);
fprintf(1,'Generated z coordinates range [%1.3e,%1.3e], with %d infs and %d nans\n',mmin,mmax,ninf,nnan);

% -------------------------------------------------------------------------
close(figure(3003));
hf = figure(3003);
hp = [];
% -------------
figure(hf);
ha = subplot(1,3,1); hold('on'); grid('on');
hp(1) = histogram(rsx,40,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',[.5,.0,.0],'LineWidth',2);
hp(2) = plot(rx,f_X,'LineStyle','--','LineWidth',2,'Color',[.0,.5,.0]);
legend(hp,'Histogram of data','Target PDF','Location','NorthWest','FontSize',16); legend('boxoff');
xlabel('x [mm]','FontSize',18);
ylabel('f_X(x) [mm^{-1}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinate X','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,2); hold('on'); grid('on');
hp(1) = histogram(rsy,40,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',[.5,.0,.0],'LineWidth',2);
hp(2) = plot(ry,f_Y,'LineStyle','--','LineWidth',2,'Color',[.0,.5,.0]);
legend(hp,'Histogram of data','Target PDF','Location','NorthWest','FontSize',16); legend('boxoff');
xlabel('y [mm]','FontSize',18);
ylabel('f_Y(y) [mm^{-1}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinate Y','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,3); hold('on'); grid('on');
hp(1) = histogram(rsz,40,'Normalization','pdf','DisplayStyle','stairs','EdgeColor',[.5,.0,.0],'LineWidth',2);
hp(2) = plot(rz,f_Z,'LineStyle','--','LineWidth',2,'Color',[.0,.5,.0]);
legend(hp,'Histogram of data','Target PDF','Location','NorthWest','FontSize',16); legend('boxoff');
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

% -------------------------------------------------------------------------
close(figure(4004));
hf = figure(4004);
% -------------
figure(hf);
ha = subplot(1,3,1); hold('on'); grid('on'); rotate3d('on');
hs = surf(rx,ry,f_XY,ones(size(f_XY)),'EdgeColor','none','FaceAlpha',0.5);
set(hs,'FaceColor','flat');
%set(hs,'CDataMode','manual');
set(hs,'CData',ones(size(f_XY)));
colormap([.5,.0,.0]);
histogram2(rsx,rsy,40,'Normalization','pdf');
xlabel('x [mm]','FontSize',18);
ylabel('y [mm]','FontSize',18);
zlabel('f_{XY}(x,y) [mm^{-2}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinates X and Y [rotate]','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,2); hold('on'); grid('on'); rotate3d('on');
hs = surf(rx,rz,f_XZ,ones(size(f_XZ)),'EdgeColor','none','FaceAlpha',0.5);
set(hs,'FaceColor','flat');
%set(hs,'CDataMode','manual');
set(hs,'CData',ones(size(f_XZ)));
colormap([.5,.0,.0]);
histogram2(rsx,rsz,40,'Normalization','pdf');
xlabel('x [mm]','FontSize',18);
ylabel('z [mm]','FontSize',18);
zlabel('f_{XZ}(x,z) [mm^{-2}]','FontSize',18);
set(ha,'FontSize',14);
title('Marginal PDF in coordinates X and Z [rotate]','FontSize',18);
drawnow;
% -------------
figure(hf);
ha = subplot(1,3,3); hold('on'); grid('on'); rotate3d('on');
hs = surf(ry,rz,f_YZ,ones(size(f_YZ)),'EdgeColor','none','FaceAlpha',0.5);
set(hs,'FaceColor','flat');
%set(hs,'CDataMode','manual');
set(hs,'CData',ones(size(f_YZ)));
colormap([.5,.0,.0]);
histogram2(rsy,rsz,40,'Normalization','pdf');
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

