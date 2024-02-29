function e = dmri_compute_convolution_weights_ODF(bi,lpar,lperp,L)
% function e = dmri_compute_convolution_weights_ODF(bi,lpar,lperp,L)
%
%     Computes the convolution weights to be applied to the SH coefficients
%     of a given ODF to reproduce the corresponding attenuation signal
%     provided a micro-structural diffusion model described by lpar and
%     lperp, at each shell with b-value bi and for each l=0,2,...,L
%
%         bi: Gx1, the b-values describing each shell
%         lpar: Nx1, the parallel eigenvalue describing the
%             micro-structure.
%         lperp: Nx1, the perpendicular eigenvalue describing the
%             micro-structure.
%         L: 1x1, an even, non-negative integer with the maximum order to
%             considered in the SH expansion.
%
%         e: Nx(L/2+1)xG, the convolution factors to be applied:
%             - at each voxel n=1...N
%             - at each group of coefficients c_l^m, l=0,2,...,L
%             - at each shell g=1...G
%
%     NOTE: no sanity checks are perfomed over any one of the input
%     arguments, so it is the user's responsibility to input appropriate
%     sizes and types.

% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------
N   = size(lpar,1);
G   = size(bi,1);
e   = zeros(N,L/2+1,G);    % Nx(L/2+1)xG
bi  = reshape(bi,[1,1,G]); % 1x1xG
if(is_broadcast_available)
    bid = (lpar-lperp).*bi;                         % Nx1xG
    crr = exp(-lperp.*bi);                          % Nx1xG
else
    bid = bsxfun( @(x,y)(x.*y), (lpar-lperp), bi ); % Nx1xG
    crr = bsxfun( @(x,y)(x.*y), lperp, bi );        % Nx1xG
    crr = exp(-crr);
end
bid(bid<0) = 0;            % Nx1xG
bid2 = sqrt(bid);          % Nx1xG
% bid  = b_i*\delta_{\lambda}
% bid2 = \sqrt{ b_i*\delta_{\lambda} }
% crr  = exp(-b_i*\lambda_{\perp})
%%% -----
pp   = (bid<100*sqrt(eps)); % Correct for very small values (Unlikely)
%%% -----
ebid2     = (sqrt(pi)/2)*erf(bid2)./bid2; % Nx1xG
ebid2(pp) = 1-bid(pp)/3;
%%% -----
r0       = (ebid2-exp(-bid))./bid;  % Nx1xG
r0(pp)   = (2/3 - 2*bid(pp)/5);     % Nx1xG
r0       = 2*pi*crr.*r0;            % Nx1xG
%%% -----
e(:,1,:) = 4*pi*ebid2.*crr;         % Nx1xG
% Recursive rule:
for l = 2:2:L
    e(:,l/2+1,:) = ((2*l-1)/l)*r0- ((l-1)/l)*e(:,l/2,:); % Nx1xG
    li = l+1; % Odd li
    r0 = ((li-1/2)./bid).*e(:,l/2+1,:) + r0;
end
% Correct for those values where b*delta_lambda is nearly 0:
bad = (bid<sqrt(eps)); % Nx1xG
if(any(bad(:)))
    badW = repmat(bad,[1,(L/2+1),1]); % Nx(L/2+1)xG
    e(badW) = 0;                      % Nx(L/2+1)xG
    e2   = e(:,1,:);                  % Nx1xG
    e2(bad) = 4*pi*crr(bad);
    e(:,1,:) = e2;
end
