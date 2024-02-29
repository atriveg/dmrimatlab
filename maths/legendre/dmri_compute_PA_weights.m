function dmri_compute_PA_weights(L,rho)
% function dmri_compute_PA_weights(L[,rho])
%
%    Computes the weights associated to the calculation of the norm of a
%    diffusion signal (or its isotropic equivalent) given by a linear
%    convolutional model expressed in terms of Spherical Harmonics
%
%       rho: 1xN double, the values off rho = lperp/(lpa-lperp) that will
%            provide the interpolation values.
%       L: 1x1 non-negative, even integer. The maximum order of the
%            Spherical harmonics expansion that will be used.
%
%    Note this function is meant to be run only once. There is no output
%    value because it stores the precomputed values in a matfile named
%    dmri_PA_weights.mat, so that this function will be run only if a value
%    of l greater than the one stored in the matfile is requested.

% -------------------------------------------------------------------------
assert(abs(L-round(L))<100*eps,'l must be integer');
assert(L>=2,'l must be >= 2');
assert(abs(L/2-round(L/2))<0.1,'l must be even');
% -------------------------------------------------------------------------
use_parallel = use_parallel_test;
% -------------------------------------------------------------------------
if(nargin<2)
    rho = (1:0.5:9.5);
    rho = [ ...
        rho/1.0e12,rho/1.0e11,rho/1.0e10,rho/1.0e9, ...
        rho/1.0e8,rho/1.0e7,rho/1.0e6,rho/1.0e5, ...
        rho/1.0e4,rho/1.0e3,rho/1.0e2,rho/1.0e1, ...
        rho*1.0e0,rho*1.0e1,rho*1.0e2,rho*1.0e3,rho*1.0e4];
else
    rho = rho(:)';
end
% -------------------------------------------------------------------------
TH  = 1000*eps;
N   = size(rho,2);
wPA = zeros(L/2+1,N);
for l=0:L/2
    weightsl = zeros(1,N);
    if(use_parallel)
        for n=1:N
            funct       = @(x1,x2)(integrand(x1,x2,rho(n),2*l));
            weightsl(n) = 4*integral2( funct, 0, 1, 0, 1, ...
                'AbsTol', 1.0e-12, 'RelTol', 1.0e-8 );
        end
    else 
        for n=1:N
            funct       = @(x1,x2)(integrand(x1,x2,rho(n),2*l));
            weightsl(n) = 4*integral2( funct, 0, 1, 0, 1, ...
                'AbsTol', 1.0e-12, 'RelTol', 1.0e-8 );
        end
    end
    p  = find(weightsl<TH,1,'first');
    if(~isempty(p))
        weightsl(p:end) = 0;
    end
    wPA(l+1,:) = weightsl;
end
rhol = log(rho);
wPAl = log(wPA);
% -------------------------------------------------------------------------
% Save the results:
file = mfilename('fullpath');
path = strrep(file,'dmri_compute_PA_weights','dmri_PA_weights');
path = [path,'.mat'];
save(path,'rho','wPA','rhol','wPAl','L');
% -------------------------------------------------------------------------

function f = integrand(x1,x2,rho,l)
sz = size(x1);
f1 = legendre(l,x1(:));
f1 = reshape(f1(1,:),sz);
f2 = legendre(l,x2(:));
f2 = reshape(f2(1,:),sz);
f3 = ( 2*rho + x1.*x1 + x2.*x2 );
f3 = sqrt(f3.*f3.*f3);
f  = (f1.*f2)./f3;


