Q   = [2.0,-1.0;-1.0,2.0];
%f   = [6.0;0.0];
f   = [-3.0;0.0];
B   = [1.0,1.0];
b   = 5.0;
%A   = -[1.0,0.0;0.0,-1.0;-1.0,-2.0];
A   = -[1.0,0.0;0.0,1.0;-1.0,-1.0];
a   = -[0.0;1.0;-4.0];
lb  = [];
ub  = [];

close(figure(1))
hf = figure(1);
hold('on');
grid('on');

x = quadprog(Q,f,A,a,B,b,lb,ub,[],optimoptions(@quadprog,'Display','none'));
plot(x(1),x(2),'ro');

x2 = -Q\f;
plot(x2(1),x2(2),'k*');

fprintf(1,'------------------------------------\n');

[x,res,iters,ll] = dmriQuadprog2(...
    Q,f,A,a,B,b,[],[]);

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

[x,status,iters,cost,lambda] = mexDMRIquadprog(Q,f,A,a,B,b,lb,ub,1000,step);


fprintf(1,'Returned status %i with %d iters and cost %1.5f\n',status,iters,cost);


CI = size(A,1);
CE = size(B,1);
CL = 0;
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
