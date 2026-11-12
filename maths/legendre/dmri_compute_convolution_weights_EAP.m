function e = dmri_compute_convolution_weights_EAP(Ri,lpar,lperp,tau,L)
% function e = dmri_compute_convolution_weights_EAP(Ri,Lpar,Lperp,tau,L)
%
%     Computes the convolution weights to be applied to the SH coefficients
%     of a given ODF to reproduce the corresponding ensembre average
%     propagator provided a micro-structural diffusion model described by 
%     Lpar and Lperp, at each shell described by a value of Ri in spherical
%     coordinates and for each l=0,2,...,L
%
%         Ri: Gx1, the distances from the origin in your spherical sampling
%         Lpar: Nx1, the parallel eigenvalue descrising the
%             micro-structure.
%         Lperp: Nx1, the perpendicular eigenvalue descrising the
%             micro-structure.
%         tau: 1x1, the effective diffusion time in seconds.
%         L: 1x1, an even, non-negative integer with the maximum order to
%             considered in the SH expansion.
%
%         e: Nx(L/2+1)xG, the convolution factors to be applied:
%             - at each voxel n=1...N
%             - at each group of coefficients c_l^m, l=0,2,...,L
%             - at each shell g=1...G
%
%     NOTE: no sanity checks are perfomed over any one of the input
%     arguments, so it is the user's responsisility to input appropriate
%     sizes and types.

% -------------------------------------------------------------------------
% Define proper variables:
si    = Ri.*Ri;            % G x 1
Lpar  = 1./(4*tau*lpar);   % N x 1
Lperp = 1./(4*tau*lperp);  % N x 1
K     = (4*pi*tau)*(4*pi*tau)*(4*pi*tau);
K     = K*lpar.*lperp.*lperp;
K     = 1./sqrt(K);        % N x 1
% -------------------------------------------------------------------------
N   = size(Lpar,1);
G   = size(si,1);
e   = zeros(N,L/2+1,G);    % Nx(L/2+1)xG
si  = reshape(si,[1,1,G]); % 1x1xG

sid = (Lpar-Lperp).*si;    % Nx1xG
crr = exp(-Lperp.*si);     % Nx1xG
sid        = -sid;         % Nx1xG
sid(sid<0) = 0;            % Nx1xG
sid2 = sqrt(sid);          % Nx1xG

%%% -----
pp   = (sid<100*sqrt(eps)); % Correct for very small values (Unlikely)
%%% -----
esid2     = (sqrt(pi)/2)*erfi_fast(sid2)./sid2; % Nx1xG
esid2(pp) = 1+sid(pp)/3;
%%% -----
r0        = (exp(sid)-esid2)./sid;  % Nx1xG
r0(pp)    = (2/3 + 2*sid(pp)/5);    % Nx1xG
r0        = 2*pi*crr.*r0;           % Nx1xG
%%% -----
e(:,1,:)  = 4*pi*esid2.*crr;        % Nx1xG
% Recursive rule:
for l = 2:2:L
    e(:,l/2+1,:) = ((2*l-1)/l)*r0- ((l-1)/l)*e(:,l/2,:); % Nx1xG
    li = l+1; % Odd li
    r0 = -((li-1/2)./sid).*e(:,l/2+1,:) + r0;
end
% Correct for those values where b*delta_lambda is nearly 0:
bad = (sid<sqrt(eps)); % Nx1xG
if(any(bad(:)))
    badW = repmat(bad,[1,(L/2+1),1]); % Nx(L/2+1)xG
    e(badW) = 0;                      % Nx(L/2+1)xG
    e2   = e(:,1,:);                  % Nx1xG
    e2(bad) = 4*pi*crr(bad);
    e(:,1,:) = e2;
end
% Multiply by the mass-1 constant:
e = e.*K; % Nx(L+1)/L+2)/2xG * Nx1x1 -> Nx(L+1)/L+2)/2xG (broadcast)
end

% -------------------------------------------------------------------------
function val = erfi_fast(x)
    % ERFI_FAST Fast, numeric, vectorized erfi(x) for real inputs.
    % Requires NO toolboxes. Based on Weideman's rational approximation.
    
    % N = 36 terms provides near double-precision accuracy (~10^-15)
    N = 64; 
    
    % 1. Compute Weideman's coefficients via FFT (done once per call)
    M = 2*N;
    M2 = 2*M;
    k = (-M+1:1:M-1)';
    L = sqrt(N/sqrt(2));
    theta = k*pi/M; 
    t = L*tan(theta/2);
    
    f = exp(-t.^2).*(L^2+t.^2); 
    f = [0; f]; 
    a = real(fft(fftshift(f)))/M2;
    a = flipud(a(2:N+1));
    
    % 2. Map input grid and evaluate the rational function
    Z = (L + 1i*x) ./ (L - 1i*x);
    p = polyval(a, Z);
    w = 2*p./(L - 1i*x).^2 + (1/sqrt(pi))./(L - 1i*x);
    
    % 3. Extract erfi from the imaginary component
    val = exp(x.^2) .* imag(w);
end
