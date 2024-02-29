% test_PosODFsGrads.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
rng(13061981);
Ni     = [31,31,31];
bi     = [1000,2000,3000];
N      = sum(Ni);
L      = 6;
K      = (L+1)*(L+2)/2;
Kp     = (2*L+1)*(2*L+2)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E      = 2+rand(N,1);
% ---
psi    = rand(K,1);
psi(1) = 1;
psi    = sqrt(0.8793)*psi./sqrt(sum(psi.*psi));
% ---
lambda = abs( ones(1,N)./((0:2:2*L)'+1) + 0.03*rand(L+1,N) );
pshell = (0:N-1)';
% ---
mu = 237;
% ---
nu = 0.1;
% ---
gi =[];
for s=1:length(Ni)
    gi = [ gi; ...
        designGradients(Ni(s),'verbose', false,'plot',false) ]; %#ok<AGROW>
end
% ---
[lagrangian,gradient,hessian,residual,jacobian,Ylm3] = ...
    mexTestPosODFsGrads(E,psi,lambda,pshell,mu,nu,gi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ylm  = GenerateSHMatrix( L, gi );
Ylm2 = GenerateSHMatrix( 2*L, gi );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi  = (sh2squaredsh(psi'))';
s1 = Ylm*psi;
s1 = s1.*s1;
s2 = Ylm2*phi;
phi2 = zeros(Kp,N);
pos = 1;
lap = 0;
for ll=0:2:2*L
    for mm=-ll:ll
        phi2(pos,:) = phi(pos)*lambda(ll/2+1,:);
        tmp = -ll*(ll+1)*phi(pos);
        lap = lap + tmp*tmp; 
        pos = pos+1;
    end
end
reconst = sum(Ylm2.*(phi2'),2);
cost = E-reconst;
cost = sum(cost.*cost)/2 + nu*lap/2;
lagrangianN = cost + mu*(sum(psi.*psi)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradientN = zeros(K+1,1);
for k=1:K
    psin = psi;
    psip = psi;
    dp   = max(1000*eps,abs(psi(k))/10000000);
    psin(k) = psi(k) + dp/2;
    psip(k) = psi(k) - dp/2;
    lagrangiann = mexTestPosODFsGrads(E,psin,lambda,pshell,mu,nu,gi);
    lagrangianp = mexTestPosODFsGrads(E,psip,lambda,pshell,mu,nu,gi);
    gradientN(k) = ( lagrangiann - lagrangianp )/dp;
end
dmu = max( 1000000*eps, abs(mu)/1000000 );
lagrangiann = mexTestPosODFsGrads(E,psi,lambda,pshell,mu+dmu/2,nu,gi);
lagrangianp = mexTestPosODFsGrads(E,psi,lambda,pshell,mu-dmu/2,nu,gi);
gradientN(K+1) = ( lagrangiann - lagrangianp )/dmu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hessianN = zeros(K+1,K+1);
for k=1:K
    psin = psi;
    psip = psi;
    dp   = max(1000*eps,abs(psi(k))/10000000);
    psin(k) = psi(k) + dp/2;
    psip(k) = psi(k) - dp/2;
    [~,gradientn] = mexTestPosODFsGrads(E,psin,lambda,pshell,mu,nu,gi);
    [~,gradientp] = mexTestPosODFsGrads(E,psip,lambda,pshell,mu,nu,gi);
    hessianN(:,k) = ( gradientn - gradientp )/dp;
end
dmu = max( 1000000*eps, abs(mu)/1000000 );
[~,gradientn] = mexTestPosODFsGrads(E,psi,lambda,pshell,mu+dmu/2,nu,gi);
[~,gradientp] = mexTestPosODFsGrads(E,psi,lambda,pshell,mu-dmu/2,nu,gi);
hessianN(:,K+1) = ( gradientn - gradientp )/dmu;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lagrangian,lagrangianN], %#ok<NOPTS>
[gradient';gradientN'],   %#ok<NOPTS>
% -----
close(figure(1));
hf1 = figure(1);
subplot(1,2,1);
imshow([hessianN,hessian],[]);
colormap(jet);
colorbar;
subplot(1,2,2);
imshow((hessian-hessianN)./(abs(hessian)+10000*eps),[]);
colormap(jet);
colorbar;
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,residualN] = computeVectorApprox(E,psi,lambda,Ylm2,nu);
jacobianN = zeros( N+Kp, K-1 );
for k=2:K
    psin = psi;
    psip = psi;
    dp   = max(1000*eps,abs(psi(k))/10000000);
    psin(k) = psi(k) + dp/2;
    psip(k) = psi(k) - dp/2;
    appn = computeVectorApprox(E,psin,lambda,Ylm2,nu);
    appp = computeVectorApprox(E,psip,lambda,Ylm2,nu);
    jacobianN(:,k-1) = (appn-appp)/dp;
end
% -----
[residual,residualN], %#ok<NOPTS>
% -----
close(figure(2));
hf2 = figure(2);
subplot(1,2,1);
imshow([jacobianN(1:N,:),jacobian(1:N,:)],[]);
colormap(jet);
colorbar;
subplot(1,2,2);
imshow((jacobianN(1:N,:)-jacobian(1:N,:))./(abs(jacobianN(1:N,:))+10000*eps),[]);
colormap(jet);
colorbar;
drawnow;
% -----
close(figure(3));
hf3 = figure(3);
subplot(1,2,1);
imshow([jacobianN(N+1:end,:),jacobian(N+1:end,:)],[]);
colormap(jet);
colorbar;
subplot(1,2,2);
imshow((jacobianN(N+1:end,:)-jacobian(N+1:end,:))./(abs(jacobianN(N+1:end,:))+10000*eps),[]);
colormap(jet);
colorbar;
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REP    = 4000;
% htrace = zeros(1,REP);
% jnorms = zeros(1,REP);
% gnorms = zeros(1,REP);
% resids = zeros(1,REP);
% for r=1:REP
%     [lagrangian,gradient,hessian,residual,jacobian,Ylm3] = ...
%         mexTestPosODFsGrads(E,psi,pshell,lambda,mu,nu,gi);
%     htrace(1,r) = trace(hessian);
%     jnorms(1,r) = trace(jacobian'*jacobian);
%     gnorms(1,r) = norm(gradient);
%     resids(1,r) = norm(residual);
% end
% close(figure(4));
% hf4 = figure(4);
% subplot(2,2,1);
% grid('on');
% plot(1:REP,htrace);
% xlabel('iters');
% ylabel('htrace');
% subplot(2,2,2);
% grid('on');
% plot(1:REP,jnorms);
% xlabel('iters');
% ylabel('jnorms');
% subplot(2,2,3);
% grid('on');
% plot(1:REP,gnorms);
% xlabel('iters');
% ylabel('gnorms');
% subplot(2,2,4);
% grid('on');
% plot(1:REP,resids);
% xlabel('iters');
% ylabel('resids');

function [approx,Q] = computeVectorApprox(E,psi,lambda,Y,nu)
% E:      N x 1
% psi:    K x 1
% lambda: (L+1) x N
% Y:      N x Kp
% nu:     1 x 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L    = size(lambda,1)-1;
N    = size(lambda,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi(1) = sqrt( 1 - sum(psi(2:end).*psi(2:end)) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi  = sh2squaredsh(psi'); % 1 x Kp
Kp   = size(phi,2);
phi2 = zeros(N,Kp);
pos  = 1;
for l=0:2:2*L
    phi2(:,pos:pos+(2*l+1)-1) = (lambda(l/2+1,:)').*phi(1,pos:pos+(2*l+1)-1);
    pos = pos + (2*l+1);
end
approx1 = sum(phi2.*Y,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ML      = GenerateSHEigMatrix( 2*L ); % Kp x Kp
approx2 = sqrt(nu)*ML*(phi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
approx = [approx1;approx2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
approx1 = E-approx1;
Q       = ( sum(approx1.*approx1) + sum(approx2.*approx2) )/2;
end
