% test_sh2squaredsh.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(13061981);
L = 8;
N = 16;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sh1 = zeros(N,(L+1)*(L+2)/2);
sh1(:,1) = 1;
pos = 2;
sig = 1;
for l=2:2:L
    sh1(:,pos:pos+2*l) = sig*randn(N,2*l+1);
    sig = sig/2;
    pos = pos + (2*l+1);
end
sh1 = sh1./sqrt(sum(sh1.*sh1,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask = ones(N,1);
mask(3,:) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sh2  = sh2squaredsh(sh1,mask);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Gi = designGradients( 7 );
B1 = GenerateSHMatrix( L,   Gi );
B2 = GenerateSHMatrix( 2*L, Gi );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sig1 = sh1*(B1');
sig2 = sh2*(B2');
sig1.*sig1,
sig2,
sum(sh1.*sh1,2)',
sqrt(4*pi)*sh2(:,1)',
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%