% test_att_modules
clear;
close('all');
rng(25042020);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lpar  = [ 3.0e-3, 2.0e-3,  1.0e-3 ];
lperp = [ 0.1e-3, 0.45e-3, 0.9e-3 ];
[lpar,lperp] = meshgrid(lpar,lperp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[gi,~,~] = icosamplesSphere(7,'O1',true); % 481 directions
L = 8;
% Compute the LS matix for SH fitting:
B   = GenerateSHMatrix( L, gi ); % G x (L+1)(L+2)/2
WLS = (B'*B)\(B');               % (L+1)*(L+2)/2 x G
%%% -----------------------------------------------------------------------
l1 = 2.0e-3;  u1 = randn(3,1);       u1 = u1./norm(u1);
l2 = 0.5e-3;  u2 = [u1(2);-u1(1);0]; u2 = u2./norm(u2);
l3 = 0.6e-3;  u3 = cross(u1,u2);     u3 = u3./norm(u3);
D1  = l1*(u1*u1') + l2*(u2*u2') + l3*(u3*u3');
D1  = D1\eye(3);
T1  = D1(1,1)*(gi(:,1).*gi(:,1)) + 2*D1(1,2)*(gi(:,1).*gi(:,2)) + ...
    2*D1(1,3)*(gi(:,1).*gi(:,3)) +   D1(2,2)*(gi(:,2).*gi(:,2)) + ...
    2*D1(2,3)*(gi(:,2).*gi(:,3)) +   D1(3,3)*(gi(:,3).*gi(:,3));
sg1 = 1./sqrt(T1.*T1.*T1); % G x 1
% -----------------
SH1 = WLS*sg1;             % (L+1)*(L+2)/2 x 1
SH1 = SH1*(1/sqrt(4*pi))/SH1(1);
%%% -----------------------------------------------------------------------
l1 = 2.4e-3;  u1 = randn(3,1);       u1 = u1./norm(u1);
l2 = 0.2e-3;  u2 = [u1(2);-u1(1);0]; u2 = u2./norm(u2);
l3 = 0.3e-3;  u3 = cross(u1,u2);     u3 = u3./norm(u3);
D1  = l1*(u1*u1') + l2*(u2*u2') + l3*(u3*u3');
D1  = D1\eye(3);
T1  = D1(1,1)*(gi(:,1).*gi(:,1)) + 2*D1(1,2)*(gi(:,1).*gi(:,2)) + ...
    2*D1(1,3)*(gi(:,1).*gi(:,3)) +   D1(2,2)*(gi(:,2).*gi(:,2)) + ...
    2*D1(2,3)*(gi(:,2).*gi(:,3)) +   D1(3,3)*(gi(:,3).*gi(:,3));
sg1 = 0.55./sqrt(T1.*T1.*T1); % G x 1
% -----------------
u10 = u1; u20 = u2; u30 = u3;
% -----------------
l1 = 2.3e-3;  u1 = u20;
l2 = 0.4e-3;  u2 = u10;
l3 = 0.3e-3;  u3 = u30;
D1  = l1*(u1*u1') + l2*(u2*u2') + l3*(u3*u3');
D1  = D1\eye(3);
T1  = D1(1,1)*(gi(:,1).*gi(:,1)) + 2*D1(1,2)*(gi(:,1).*gi(:,2)) + ...
    2*D1(1,3)*(gi(:,1).*gi(:,3)) +   D1(2,2)*(gi(:,2).*gi(:,2)) + ...
    2*D1(2,3)*(gi(:,2).*gi(:,3)) +   D1(3,3)*(gi(:,3).*gi(:,3));
sg1 = sg1 + 0.45./sqrt(T1.*T1.*T1); % G x 1
% -----------------
SH2 = WLS*sg1;             % (L+1)*(L+2)/2 x 1
SH2 = SH2*(1/sqrt(4*pi))/SH2(1);
%%% -----------------------------------------------------------------------
l1 = 2.7e-3;  u1 = randn(3,1);       u1 = u1./norm(u1);
l2 = 0.2e-3;  u2 = [u1(2);-u1(1);0]; u2 = u2./norm(u2);
l3 = 0.1e-3;  u3 = cross(u1,u2);     u3 = u3./norm(u3);
D1  = l1*(u1*u1') + l2*(u2*u2') + l3*(u3*u3');
D1  = D1\eye(3);
T1  = D1(1,1)*(gi(:,1).*gi(:,1)) + 2*D1(1,2)*(gi(:,1).*gi(:,2)) + ...
    2*D1(1,3)*(gi(:,1).*gi(:,3)) +   D1(2,2)*(gi(:,2).*gi(:,2)) + ...
    2*D1(2,3)*(gi(:,2).*gi(:,3)) +   D1(3,3)*(gi(:,3).*gi(:,3));
sg1 = 0.25./sqrt(T1.*T1.*T1); % G x 1
% -----------------
u10 = u1; u20 = u2; u30 = u3;
% -----------------
l1 = 2.8e-3;   u1 = u20;
l2 = 0.25e-3;  u2 = u10;
l3 = 0.25e-3;  u3 = u30;
D1  = l1*(u1*u1') + l2*(u2*u2') + l3*(u3*u3');
D1  = D1\eye(3);
T1  = D1(1,1)*(gi(:,1).*gi(:,1)) + 2*D1(1,2)*(gi(:,1).*gi(:,2)) + ...
    2*D1(1,3)*(gi(:,1).*gi(:,3)) +   D1(2,2)*(gi(:,2).*gi(:,2)) + ...
    2*D1(2,3)*(gi(:,2).*gi(:,3)) +   D1(3,3)*(gi(:,3).*gi(:,3));
sg1 = sg1 + 0.3./sqrt(T1.*T1.*T1); % G x 1
% -----------------
l1 = 2.9e-3;   u1 = u30;
l2 = 0.1e-3;   u2 = u20;
l3 = 0.15e-3;  u3 = u10;
D1  = l1*(u1*u1') + l2*(u2*u2') + l3*(u3*u3');
D1  = D1\eye(3);
T1  = D1(1,1)*(gi(:,1).*gi(:,1)) + 2*D1(1,2)*(gi(:,1).*gi(:,2)) + ...
    2*D1(1,3)*(gi(:,1).*gi(:,3)) +   D1(2,2)*(gi(:,2).*gi(:,2)) + ...
    2*D1(2,3)*(gi(:,2).*gi(:,3)) +   D1(3,3)*(gi(:,3).*gi(:,3));
sg1 = sg1 + 0.45./sqrt(T1.*T1.*T1); % G x 1
% -----------------
SH3 = WLS*sg1;             % (L+1)*(L+2)/2 x 1
SH3 = SH3*(1/sqrt(4*pi))/SH3(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph = createIcosahedron;
for n=1:4
    ph = refinePolyhedron(ph);
end
vts = ph.vertices;
NF  = size(ph.facets);
p1  = vts(ph.facets(:,1),:);
p2  = vts(ph.facets(:,2),:);
p3  = vts(ph.facets(:,3),:);
a   = sqrt(sum((p1-p2).*(p1-p2),2));
b   = sqrt(sum((p2-p3).*(p2-p3),2));
c   = sqrt(sum((p3-p1).*(p3-p1),2));
s   = (a+b+c)/2;
A   = sqrt(s.*(s-a).*(s-b).*(s-c));
B2  = GenerateSHMatrix( L, vts );
% -----------------------
close(figure(1));
figure(1);
% -----------------------
subplot(1,3,1);
odf1 = B2*SH1;
% ---------
P = mean([odf1(ph.facets(:,1)),odf1(ph.facets(:,2)),odf1(ph.facets(:,3))],2);
P = sum(P.*A);
% ---------
vts1  = vts.*repmat(odf1,[1,3]);
patch( 'Vertices', vts1, ...
    'Faces', ph.facets, 'FaceColor', [0,0.5,0], ...
    'FaceAlpha', 0.8, 'EdgeColor', 'none' );
axis('equal');
light;
lighting(phong);
title(['Estimated prob. :',num2str(P)]);
% -----------------------
subplot(1,3,2);
odf2 = B2*SH2;
% ---------
P = mean([odf2(ph.facets(:,1)),odf2(ph.facets(:,2)),odf2(ph.facets(:,3))],2);
P = sum(P.*A);
% ---------
vts2  = vts.*repmat(odf2,[1,3]);
patch( 'Vertices', vts2, ...
    'Faces', ph.facets, 'FaceColor', [0,0.5,0], ...
    'FaceAlpha', 0.8, 'EdgeColor', 'none' );
axis('equal');
light;
lighting(phong);
title(['Estimated prob. :',num2str(P)]);
% -----------------------
subplot(1,3,3);
odf3 = B2*SH3;
% ---------
P = mean([odf3(ph.facets(:,1)),odf3(ph.facets(:,2)),odf3(ph.facets(:,3))],2);
P = sum(P.*A);
% ---------
vts3  = vts.*repmat(odf3,[1,3]);
patch( 'Vertices', vts3, ...
    'Faces', ph.facets, 'FaceColor', [0,0.5,0], ...
    'FaceAlpha', 0.8, 'EdgeColor', 'none' );
axis('equal');
light;
lighting(phong);
title(['Estimated prob. :',num2str(P)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lpar  = repmat(lpar,[1,1,3]);     % 3 x 3 x 3
lperp = repmat(lperp,[1,1,3]);    % 3 x 3 x 3
K     = (L+1)*(L+2)/2;
SH    = zeros(3,3,3,K);           % 3 x 3 x 3 x (L+1)*(L+2)/2
SH(:,:,1,:) = repmat(reshape(SH1,[1,1,1,K]),[3,3,1,1]);
SH(:,:,2,:) = repmat(reshape(SH2,[1,1,1,K]),[3,3,1,1]);
SH(:,:,3,:) = repmat(reshape(SH3,[1,1,1,K]),[3,3,1,1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------
tau = 35.0e-3;
b   = [0:0.5:100,101:1:200,202:2:500,505:5:1000,1010:10:2000,2020:20:3000,3050:50:5000,5100:100:10000,10200:200:15000,15500:500:20000,21000:1000:30000];
q   = sqrt(b/(4*pi*pi*tau));
qi  = ( q(2:end).*q(2:end).*q(2:end) - q(1:end-1).*q(1:end-1).*q(1:end-1) )/3;
qc  = ( q(2:end) + q(1:end-1) )/2;
bc  = ( b(2:end) + b(1:end-1) )/2;
% -----------
NG  = size(gi,1);
l1  = 1.7e-3;
l2  = 0.4e-3;
l3  = 0.2e-3;
u1  = [1;-1;0]/sqrt(2);
u2  = [1;1;0]/sqrt(2);
u3  = [0;0;1];
DT  = l1*(u1*u1') + l2*(u2*u2') + l3*(u3*u3');
QF  = zeros(NG,1);
for n=1:NG
    QF(n) = gi(n,:)*DT*gi(n,:)';
end
% -----------
modE = zeros(3,3,3);
modG = 0;
modT = 0;
cpEG = zeros(3,3,3);
AVM  = GenerateSHMatrix( L, gi ); % NG x K
AVM  = (AVM'*AVM)\(AVM');         % K x NG
AVM  = AVM(1,:)*sqrt(4*pi);       % 1 x NG
% -----------
NQ      = length(qi);
CHUNKSZ = 4;
Nstr    = 0;
for ck=1:ceil(NQ/CHUNKSZ)
    % ---------------------------------------------
    idi = (ck-1)*CHUNKSZ+1;
    idf = min(ck*CHUNKSZ,NQ);
    NC  = idf-idi+1;
    % ---------------------------------------------
    Gi  = repmat(gi,[NC,1]);                   % NC*NG x 3
    bs  = repmat(bc(idi:idf),[size(gi,1),1]);  % NG x NC
    bi  = bs(:);                               % NG*NC x 1 
    % ---------------------------------------------
    attiE = micro2atti(SH,lpar,lperp,[],Gi,bi,'bth',0.5); % 3 x 3 x 3 x (NG*NC)
    attiE = permute(attiE,[4,1,2,3]);                  % (NG*NC) x 3 x 3 x 3
    attiE = reshape(attiE,[NG,NC,3,3,3]);              % NG x NC x 3 x 3 x 3
    % ---------------------------------------------
    attiG = exp(-QF*bc(idi:idf));              % NG x NC
    attiT = ones(NG,1)*exp(-qc(idi:idf)/7);    % NG x NC 
    % ---------------------------------------------
    attiE2 = attiE.*attiE; % NG x NC x 3 x 3 x 3 
    attiG2 = attiG.*attiG; % NG x NC
    for x=1:3
        for y=1:3
            for z=1:3
                modE(x,y,z) = modE(x,y,z) + qi(idi:idf)*(AVM*attiE2(:,:,x,y,z))';
                cpEG(x,y,z) = cpEG(x,y,z) + ...
                    qi(idi:idf)*(AVM*(attiE(:,:,x,y,z).*attiG))';
            end
        end
    end
    modG = modG + qi(idi:idf)*(AVM*attiG2)';
    modT = modT + qi(idi:idf)*(AVM*attiT)';
    % ---------------------------------------------
    fprintf(1,repmat('\b',[1,Nstr]));
    str  = sprintf('Processed bs up to %1.1f of a maximum of %1.1f',bc(idf),bc(end));
    fprintf(1,str);
    Nstr = length(str);
    % ---------------------------------------------
end
fprintf(1,'\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modG0 = (4*pi*tau);
modG0 = modG0*modG0*modG0;
modG0 = modG0*det(2*DT);
modG0 = 1/sqrt(modG0);
modT0 = 4*pi*686;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(exist('dmri_PA_weights.mat','file')==2)
    S = load('dmri_PA_weights.mat');
    if(S.L<20)
        dmri_compute_PA_weights(L);
        S = load('dmri_PA_weights.mat');
    end
else
    dmri_compute_PA_weights( max(L,20) );
    S = load('dmri_PA_weights.mat');
end
close(figure(2));
figure(2);
cols = [ ...
    .50,.00,0.00; ...
    .00,.50,0.00; ...
    .00,.00,0.50; ...
    .30,.30,0.00; ...
    .30,.00,0.30; ...
    .00,.30,0.30; ...
    .40,.20,0.20; ...
    .20,.40,0.20; ...
    .20,.20,0.40; ...
    .60,.10,0.30; ...
    .10,.60,0.30   ];
styles = {'-','--','-.','-','--','-.','-','--','-.','-','--'};
hold('on');
grid('on');
LMAX = 7;
hp = zeros(1,LMAX);
for n=1:LMAX
    hp(n) = plot(S.rhol,S.wPAl(n,:),'Color',cols(n,:),'LineStyle',styles{n},'LineWidth',2);
end
xlabel('$\log_{10}(\rho_{\lambda})$','FontSize',24,'Interpreter','latex','FontName','serif');
ylabel('$\log_{10}({\cal I}(\rho_{\lambda}))$','FontSize',24,'Interpreter','latex','FontName','serif');
legtxt = {'l=0','l=2','l=4','l=6','l=8','l=10','l=12','l=14','l=16','l=18','l=20'};
hl = legend(hp,legtxt{1:LMAX});
set(hl,'Location','SouthWest','FontSize',18,'FontAngle','italic','FontName','serif');
set(gca,'FontSize',14);
axis([-28,12,-30,16]);

delta = (lpar-lperp); % 3 x 3 x 3
rho   = lperp./delta; % 3 x 3 x 3
rhol  = log(rho);     % 3 x 3 x 3
rholi = S.rhol;       % 1 x NI
rhol(rhol<rholi(1)) = rholi(1);
wPAli = S.wPAl;       % (LI/2+1) x NI
wPA   = zeros([size(lpar),L/2+1]);  % 3 x 3 x 3 x (L/2+1)
for l=1:L/2+1
    wPA(:,:,:,l) = exp( interp1(rholi,wPAli(l,:),rhol) ); % 3 x 3 x 3 x 1
end
wPA   = reshape( wPA, [numel(lpar),L/2+1] ); % Q x (L/2+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH    = reshape( SH, [numel(lpar),K] );  % Q x K
rho   = rho(:);                          % Q x 1
delta = delta(:);                        % Q x 1
nfact = 4*pi*tau*delta;                  % Q x 1
nfact = pi./sqrt(nfact.*nfact.*nfact);   % Q x 1
wghts = wPA(:,dmri_sh_expand_coeffs(L)); % Q x K
wghts = wghts.*repmat(nfact,[1,K]);      % Q x K
modE0 = sum(SH.*SH.*wghts,2);            % Q x 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lev = log2(max((2*K-16)/15,1)) + 2;
vi  = icosamplesSphere(ceil(lev),...
    'O1',true,'iters',20,'verbose',false); % G x 3
l1   = l1+lperp(:); % Q x 1
l2   = l2+lperp(:); % Q x 1
l3   = l3+lperp(:); % Q x 1
u1   = repmat(u1',[numel(lpar),1]);
u2   = repmat(u2',[numel(lpar),1]);
u3   = repmat(u3',[numel(lpar),1]);
prj1 = ( u1(:,1)*vi(:,1)' + u1(:,2)*vi(:,2)' + u1(:,3)*vi(:,3)' ); % Q x G
prj2 = ( u2(:,1)*vi(:,1)' + u2(:,2)*vi(:,2)' + u2(:,3)*vi(:,3)' ); % Q x G
prj3 = ( u3(:,1)*vi(:,1)' + u3(:,2)*vi(:,2)' + u3(:,3)*vi(:,3)' ); % Q x G
if(is_broadcast_available_test)
    prj1 = prj1.*prj1./l1; % Q x G
    prj2 = prj2.*prj2./l2; % Q x G
    prj3 = prj3.*prj3./l3; % Q x G
else
    prj1 = bsxfun( @(x,y)(x./y), prj1.*prj1, l1 ); % Q x G
    prj2 = bsxfun( @(x,y)(x./y), prj2.*prj2, l2 ); % Q x G
    prj3 = bsxfun( @(x,y)(x./y), prj3.*prj3, l3 ); % Q x G
end
detp  = prj1 + prj2 + prj3;    % Q x G
if(is_broadcast_available_test)
    detp = delta.*detp + 1;    % Q x G
    detp = detp.*(l1.*l2.*l3); % Q x G
else
    detp = bsxfun(@(x,y)(x.*y),delta,detp) + 1;  % Q x G
    detp = bsxfun(@(x,y)(x.*y),detp,l1.*l2.*l3); % Q x G
end
detp = 1./sqrt(detp);               % Q x G
sh2 = signal2sh( reshape(detp,[numel(lpar),1,1,size(vi,1)]), vi, ...
    'L', L, 'lambda', 1.0e-9, 'chunksz', 100 );
sh2 = reshape(sh2,[numel(lpar),K]); % Q x K
cpEG0 = sum(SH.*sh2,2);             % Q x 1
cpEG0 = cpEG0./sqrt(64*pi*pi*pi*tau*tau*tau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'[DUMB INTEGRAL] Theoretical: %1.4f; Numerical: %1.4f\n',modT0,modT);
fprintf(1,'[TESOR MODEL] Theoretical: %1.4f; Numerical: %1.4f\n',modG0,modG);
% -------------------
close(figure(3));
figure(3);
hold('on');
grid('on');
hl(1) = plot(1:numel(lpar),modE0',  'r--x','MarkerSize',20);
hl(2) = plot(1:numel(lpar),modE(:)','g:o', 'MarkerSize',20);
legend(hl,'Theoretical','Numerical');
title('||E(q)||^2');
for n=1:numel(lpar)
    x0 = n;
    y0 = modE0(n);
    text(x0,y0+6e4,['\rho_{\lambda}=',sprintf('%1.3f',rho(n))],...
        'FontSize',8,'FontWeight','bold','Interpreter','tex',...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
    text(x0,y0+4e4,['\lambda_{||}=',sprintf('%1.3f',1000*lpar(n))],...
        'FontSize',8,'FontWeight','bold','Interpreter','tex',...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
    text(x0,y0+2e4,['\lambda_{\perp}=',sprintf('%1.3f',1000*lperp(n))],...
        'FontSize',8,'FontWeight','bold','Interpreter','tex',...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
end
axis([1,numel(lpar),0,3.5e5]);
% -------------------
close(figure(4));
figure(4);
hold('on');
grid('on');
hl(1) = plot(1:numel(lpar),cpEG0',  'r--x','MarkerSize',20);
hl(2) = plot(1:numel(lpar),cpEG(:)','g:o', 'MarkerSize',20);
legend(hl,'Theoretical','Numerical');
title('<E(q),E_Gauss(q)>');
for n=1:numel(lpar)
    x0 = n;
    y0 = cpEG0(n);
    text(x0,y0+1.5e4,['\rho_{\lambda}=',sprintf('%1.3f',rho(n))],...
        'FontSize',8,'FontWeight','bold','Interpreter','tex',...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
    text(x0,y0+1.0e4,['\lambda_{||}=',sprintf('%1.3f',1000*lpar(n))],...
        'FontSize',8,'FontWeight','bold','Interpreter','tex',...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
    text(x0,y0+0.5e4,['\lambda_{\perp}=',sprintf('%1.3f',1000*lperp(n))],...
        'FontSize',8,'FontWeight','bold','Interpreter','tex',...
        'VerticalAlignment','bottom','HorizontalAlignment','center');
end
axis([1,numel(lpar),3e4,16e4]);
