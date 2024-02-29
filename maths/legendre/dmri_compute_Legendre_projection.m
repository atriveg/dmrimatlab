function Cl0 = dmri_compute_Legendre_projection(iota)
% function Cl0 = dmri_compute_Legendre_projection(iota)
%
%    Computes the projection of a kernel function of the type
%    1/(1+-x^2/K)^(gamma/2) on the basis of Legendre polynomials. This is
%    used together with dmri_compute_Ekernel_integrals and
%    dmri_compute_Pkernel_integrals to find convolution kernels to compute
%    moments. For example:
%
%        IOTA = dmri_compute_Ekernel_integrals(lpar,lperp,gamma,N);
%        Cl0  = dmri_compute_Legendre_projection(IOTA);
%
%    will compute the C_l^0 even SH convolution factors for the
%    rotation-invariant kernel:
%
%        k(theta,phi) = k(theta) = 1/(1+cos(theta)^2/rho)^(nu/2),
%                 rho = lperp/(lpar-lperp)
%
%    up to degree N/2
Cl0 = zeros(size(iota));
L   = (size(iota,2)-1)*2;
Pl  = dmri_legendre_polynomials(L);
for l=1:size(iota,2)
    for n=1:l
        Cl0(:,l) = Cl0(:,l) + Pl(l,n).*iota(:,n);
    end
end
Cl0 = 2*pi*Cl0;
