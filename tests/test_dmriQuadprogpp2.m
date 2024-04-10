N = 11;
a = (0:N-1)/N*2*pi+0.13;
x = [-1,4];
y = [-1,4];
close(figure(1));
figure(1);
hold('on');
grid('on');
dist = 0.6;
alpha = zeros(1,length(a));
beta  = zeros(1,length(a));
term  = zeros(1,length(a));
for n=1:length(a)
    alpha(n) = cos(a(n));
    beta(n)  = sin(a(n));
    x0    = 1 - dist*alpha(n);
    y0    = 1 - dist*beta(n);
    y = y0 - alpha(n)*(x-x0)/beta(n);
    term(n) = -alpha(n)*x0-beta(n)*y0;
    plot(x,y);
end
axis([-1,4,-1,4]);
axis('equal');

A  = [-alpha(:),-beta(:)];
a  = term(:);
B  = [1,1];
b  = 2;
lb = [-10;1.0];

plot(x,(b-B(1)*x)/B(2),'k--');
Q = [2,4;4,10];
f = [-6;-12];

x = quadprog(Q,f,A,a,B,b,lb,[],[],optimoptions(@quadprog,'Display','none'));
plot(x(1),x(2),'ro');

x2 = -Q\f;
plot(x2(1),x2(2),'k*');

fprintf(1,'------------------------------------\n');

[x,res,iters,ll] = dmriQuadprog2(...
    Q,f,A,a,B,b,lb,[]);

fprintf(1,'Matlab''s code converged within %d iters, residual %1.5f\n',iters,res);

plot(x(1),x(2),'k+');

fprintf(1,'------------------------------------\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CC = [A;B];
uu = [a;b];
Qi = Q\eye(size(Q,1));
HH = CC*Qi*CC';
ff = ((f')*Qi*(CC')+uu')';
step = sqrt(sum(ff.*ff)) / (trace(HH)/size(HH,1));
fprintf(1,'Trying step: %1.4f\n',step);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[x,status,iters,cost,lambda] = mexDMRIquadprog(Q,f,A,a,B,b,lb,[],1000,step);


fprintf(1,'Returned status %i with %d iters and cost %1.5f\n',status,iters,cost);


CI = size(A,1);
CE = size(B,1);
CL = 2;
CU = 0;

pos = 0;
fprintf(1,'lamdba_ineq: %s\n',num2str(lambda(pos+1:pos+CI)'));
pos = pos+CI;
fprintf(1,'lamdba_eq: %s\n',num2str(lambda(pos+1:pos+CE)'));
pos = pos+CE;
fprintf(1,'lamdba_lb: %s\n',num2str(lambda(pos+1:pos+CL)'));
pos = pos+CL;
fprintf(1,'lamdba_ub: %s\n',num2str(lambda(pos+1:pos+CU)'));
pos = pos+CU;
fprintf(1,'------------------------------------\n');

plot(x(1),x(2),'go','MarkerSize',16);
