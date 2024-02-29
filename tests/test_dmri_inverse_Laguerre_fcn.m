% test_dmri_inverse_Laguerre_fcn.m
N        = 1000000;
Amax     = 6;
A        = [linspace(0,Amax,N),linspace(10,50,N),linspace(50,150,N)];
A        = [A(1),exp(linspace(log(1.0e-9),log(A(2)),N)),A(2:end)];

%A        = [0.0001,0.001,0.0099,0.15];
X        = laguerreL_1_2(A);

[Ae,nit] = dmri_inverse_Laguerre_fcn(X);

close(figure(1003));
hf = figure(1003);
subplot(1,2,1); hold('on'); grid('on');
xlabel('Actual A/sigma');
ylabel('Estimated A/sigma');
plot(A,A,'k--');
plot(A,Ae,'g-');
subplot(1,2,2); hold('on'); grid('on');
xlabel('Actual A/sigma');
ylabel('Number of iterations');
plot(A,nit,'g-');
