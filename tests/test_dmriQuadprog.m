% test_dmriQuadprog.m
clear;
rng(13061981);
N = 32;
k = (0:N-1);
m = (0:N-1)';
A = cos((2*pi/N)*(m*k))/sqrt(N);
A = A(1:N/2,1:N/2);

xt  = rand(N/2,1);
xt  = -log(xt);
xt  = xt./sum(xt);
yt  = A*xt;
sig = mean(abs(yt))/10;
y   = yt + sig*randn(size(yt));

H = A'*A;
b = -A'*y;

x0 = A\y;
x0(x0<0) = 0;
x0 = x0./sum(x0);

H00 = H(1,1);
f0  = 1;
b0  = b(1);

for k=1:100
    [x,res] = dmriQuadprog(H,b,ones(size(xt)),x0,...
        k-1,1.0e-12,1.0e-12);
    ress(k)  = res + H00/f0/f0/2 + b0/f0 + 0.5*(y'*y); %#ok<SAGROW>
    ress2(k) = 0.5*x'*H*x + b'*x + 0.5*(y'*y);         %#ok<SAGROW>
    ress3(k) = 0.5*((y-A*x)')*(y-A*x);                 %#ok<SAGROW>
    c1(k) = sum(x); %#ok<SAGROW>
    c2(k) = min(x); %#ok<SAGROW>
end
close(figure(1));
hf1 = figure(1);
hold('on');
grid('on');
hp(1) = plot( 0:length(ress)-1,  ress,  'LineStyle', '-', 'Color', [0.5,0.0,0.0], 'LineWidth', 2 );
hp(2) = plot( 0:length(ress2)-1, ress2, 'LineStyle', ':', 'Color', [0.0,0.0,0.5], 'LineWidth', 2 );
hp(3) = plot( 0:length(ress3)-1, ress3, 'LineStyle', '-.', 'Color', [0.0,0.5,0.0], 'LineWidth', 2 );
xlabel('iteration');
ylabel('residual');
legend(hp,'Residual in dmriQuadrog','Residual in original quadprog','Residual in the original problem');

close(figure(2));
hf2 = figure(2);
subplot(2,1,1); plot(0:length(c1)-1,c1); xlabel('iteration'); ylabel('constraint 1');
subplot(2,1,2); plot(0:length(c2)-1,c2); xlabel('iteration'); ylabel('constraint 2');