% test_sh2signal2sh

rng(13,'twister');

L  = 8;
lambda = 1.0e-12;
K  = (L+1)*(L+2)/2;
fc = ones(1,K);
i1 = 1;
i2 = 1;
for l=2:2:L
    i1 = i2 + 1;
    i2 = i1 + 2*l;
    fc(i1:i2) = 1/(l*l/4);
end

X  = 11;
SH = rand(X*X*X,1);
SH = [SH,rand(X*X*X,K-1)-1/2];
SH = SH./repmat(fc,[X*X*X,1]);
SH = reshape(SH,[X,X,X,K]);

gi = icosamplesSphere(4,'O1',true);

signal = sh2signal( SH, gi );

SH2 = signal2sh( signal, gi, 'L', L, 'lambda', lambda );

signal2 = sh2signal( SH2, gi );

err1 = sum((SH-SH2).*(SH-SH2),4);
err2 = sum((signal-signal2).*(signal-signal2),4);

fprintf(1,'Errors with lambda = %1.4g:\n',lambda);
fprintf(1,'In the SH coefficients;      min: %1.4g; max: %1.4g; mean: %1.4g\n',min(err1(:)),max(err1(:)),mean(err1(:)));
fprintf(1,'In the reconstructed signal; min: %1.4g; max: %1.4g; mean: %1.4g\n',min(err2(:)),max(err2(:)),mean(err2(:)));
