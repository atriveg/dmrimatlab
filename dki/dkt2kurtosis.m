function [mk,ak,rk] = dkt2kurtosis( dkt, dti, varargin )
% function [MK,AK,RK] = dkt2kurtosis( dkt, dti, ...
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Computes the mean kurtosis MK, the axial kurtosis AK and the radial
%   kurtosis RK from the 4th order kurtosis tensor dkt and the 2nd order
%   diffusion tensor dti. According to the signal model:
%
%     log(atti) = log(S) - log(S_0)
%       = - b * sum_{i=1,j=1}^3 ni·nj·D_{ij}
%         + b^2/6*MD^2 * sum_{i=1,j=1,k=1,l=1}^3 ni·nj·nk·nl·W_{ijkl}
% 
%   where D_{ij} is the 2nd order diffusion tensor, W_{ijkl} is the 4th 
%   order kurtosis tensor and MD is the mean diffusivity computed from 
%   D_{ij}, the directional kurtosis for each unit direction n=[n1,n2,n3]^T 
%   reads:
%
%     K(n) = MD^2 * (sum_{i=1,j=1}^3 ni·nj·D_{ij})^{-2}
%                            * sum_{i=1,j=1,k=1,l=1}^3 ni·nj·nk·nl·W_{ijkl}
%
%   If u1 is the main eigenvector of D_{ij}, we calculate:
%
%       AK = K(u1)
%       RK = \int_{u \perp u_1} K(u) du / (2*pi)
%       MK = \int\int_{u} K(u) du / (4*pi^2)
%
%   These expressions are analytically computed using A. Tabesh's
%   developments in Tabesh et al. (2011) Mag. Res. Med. 65: 823-836.
%
%   MANDATORY INPUTS:
%
%      dkt: MxNxPx15, where the last dimension contains the 15 unique
%         components of the 4th order kurtosis tensor, W_{ijkl}. The
%         ordering of these components is consistent with that used in
%         HOT-related functions. It is returned by atti2dkt.
%      dti: MxNxPx6, where the last dimension contains the 6 unique 
%         components of the 2nd order diffusion tensor, D_{ik}. It is
%         returned by atti2dkt (or atti2dti).
%
%   OUTPUTS:
%
%      MK, AK, RK: MxNxP each, the respective values of the mean kurtosis,
%      the axial kurtosis and the radial kurtosis.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      chunksz: the computation internally uses matrix products with HOT
%         coefficients with chunks of data. This is the size of such chunks
%         (default: 1000).
%      maxthreads: the algorithm interanlly resorts to dti2spectrum, which 
%         is run with multiple threads. This is the maximum allowed number 
%         of threads, which can indeed be reduced if it exceeds the number 
%         of logical cores (default: the number of logical cores in the 
%         machine).

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<2)
    error('At lest the dkt and the dti must be supplied');
end
% ---------
[M,N,P,K]        = size(dkt);
[M2,N2,P2,K2,L2] = size(dti);
% ---------
assert(K==15,'The 4-th dimension of a proper dkt volume must have size 15');
if(K2==3)
    if(L2==3)
        dti = cat( 4, dti(:,:,:,1,1), dti(:,:,:,1,2), dti(:,:,:,1,3), ...
            dti(:,:,:,2,2), dti(:,:,:,2,3), dti(:,:,:,3,3) );
    else
        error('The dti volume provided does not look like a diffusion tensor volume');
    end
elseif(K2==6)
    if(L2~=1)
        error('The dti volume provided does not look like a diffusion tensor volume');
    end
else
    error('The dti volume provided does not look like a diffusion tensor volume');
end
% ---------
assert(isequal([M,N,P],[M2,N2,P2]),'The FoVs of dkt and dti do not match');
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true]; % always 1x1 char
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% The spectrum of the DT is required for many internal computations:
[u1,u2,u3,l1,l2,l3] = dti2spectrum( dti, 'clip', true, ...
    'mask', opt.mask, 'maxthreads', opt.maxthreads );
% We will reshape all volumes and keep only masked values to work
% comfortably:
dkt = reshape( dkt, [M*N*P,15] );
dkt = dkt( opt.mask, : );
u1  = reshape( u1, [M*N*P,3] );
u1  = u1( opt.mask, : );
u2  = reshape( u2, [M*N*P,3] );
u2  = u2( opt.mask, : );
u3  = reshape( u3, [M*N*P,3] );
u3  = u3( opt.mask, : );
l1  = reshape( l1, [M*N*P,1] );
l1  = l1( opt.mask, : );
l2  = reshape( l2, [M*N*P,1] );
l2  = l2( opt.mask, : );
l3  = reshape( l3, [M*N*P,1] );
l3  = l3( opt.mask, : );
% Start by computing the required kurtosis tensor components in the rotated
% frame aligned with the eigenvectors of the diffusion tensor:
[W1111,W2222,W3333,W1122,W1133,W2233] = W_rotated_frame(dkt,u1,u2,u3,opt.chunksz);
% Compute mean kutosis...
mk_ = ...
    tabesh_F1(l1,l2,l3,opt.maxthreads).*W1111 + ...
    tabesh_F1(l2,l1,l3,opt.maxthreads).*W2222 + ...
    tabesh_F1(l3,l2,l1,opt.maxthreads).*W3333 + ...
    tabesh_F2(l1,l2,l3,opt.maxthreads).*W2233 + ...
    tabesh_F2(l2,l1,l3,opt.maxthreads).*W1133 + ...
    tabesh_F2(l3,l2,l1,opt.maxthreads).*W1122;
% ... the axial kurtosis...
ak_ = tabesh_H1(l1,l2,l3).*W1111;
% ... and the radial kurtosis...
rk_ = ...
    tabesh_G1(l1,l2,l3).*W2222 + ...
    tabesh_G1(l1,l3,l2).*W3333 + ...
    tabesh_G2(l1,l2,l3).*W2233;
% Cast the results to the proper size:
mk = zeros(M,N,P);
ak = zeros(M,N,P);
rk = zeros(M,N,P);
mk(opt.mask) = mk_;
ak(opt.mask) = ak_;
rk(opt.mask) = rk_;
end

% Auxiliar functions implemented right from Tabesh's paper:
% -------------------------------------------------------------------------
% Rotation of HOT
function [W1111,W2222,W3333,W1122,W1133,W2233] = W_rotated_frame(dkt,u1,u2,u3,chunksz)
% dkt: N x 15
% u1:  N x 3
% u2:  N x 3
% u3:  N x 3
% W1111' = sum_{ijkl=1}^3 u1_i u1_j u1_k u1_l W_{ijkl}
% W2222' = sum_{ijkl=1}^3 u2_i u2_j u2_k u2_l W_{ijkl}
% W3333' = sum_{ijkl=1}^3 u3_i u3_j u3_k u3_l W_{ijkl}
% W1122' = sum_{ijkl=1}^3 u1_i u1_j u2_k u2_l W_{ijkl}
% W1133' = sum_{ijkl=1}^3 u1_i u1_j u3_k u3_l W_{ijkl}
% W2233' = sum_{ijkl=1}^3 u2_i u2_j u3_k u3_l W_{ijkl}
[idx,u] = unroll_hot_order2_idx;
W1111 = compute_rotatedW_component(dkt,u1,u1,idx,u,chunksz);
W2222 = compute_rotatedW_component(dkt,u2,u2,idx,u,chunksz);
W3333 = compute_rotatedW_component(dkt,u3,u3,idx,u,chunksz);
W1122 = compute_rotatedW_component(dkt,u1,u2,idx,u,chunksz);
W1133 = compute_rotatedW_component(dkt,u1,u3,idx,u,chunksz);
W2233 = compute_rotatedW_component(dkt,u2,u3,idx,u,chunksz);
end

function Wiijj = compute_rotatedW_component(dkt,ui,uj,idx,u,chunksz)
% dkt:     N x 15
% ui:      N x 3
% uj:      N x 3
% idx:     81 x 1
N     = size(dkt,1);
Wiijj = zeros(N,1);
for ck=1:ceil(N/chunksz)
    % -------------------------------
    idi = (ck-1)*chunksz+1;
    idf = min(ck*chunksz,N);
    % -------------------------------
    DKT = dkt(idi:idf,idx);             % chunksz x 81
    % -------------------------------
    Ui1 = ui( idi:idf, u(:,1) );        % chunksz x 81
    Ui2 = ui( idi:idf, u(:,2) );        % chunksz x 81
    Uj1 = uj( idi:idf, u(:,3) );        % chunksz x 81
    Uj2 = uj( idi:idf, u(:,4) );        % chunksz x 81
    % -------------------------------
    Wiijj(idi:idf,1) = sum(DKT.*Ui1.*Ui2.*Uj1.*Uj2,2); % chunksz x 1
end
end

function [idx,u] = unroll_hot_order2_idx
% Returns the 81 x 1 array idx that can be used to
% obtain the "unrolled" version of the packed HOT
% representation up top order 4 (idx contains
% values in [1,15], so that you can:
%    hot_unrolled = hot(idx)
% Finally, u is 81 x 4 and contains the indices of
% the 3 x 1 input vector that must be multiplied to
% compute each of the 81 addends
[~,~,pwrs]   = hot2signal( zeros(1,1,1,15), [1,0,0] );
[u4,u3,u2,u1] = ndgrid(1:3,1:3,1:3,1:3);
u1 = u1(:);          % 81 x 1
u2 = u2(:);          % 81 x 1
u3 = u3(:);          % 81 x 1
u4 = u4(:);          % 81 x 1
u  = [u1,u2,u3,u4];  % 81 x 4
n1 = sum(u==1,2);    % 81 x 1, how many 1's?
n2 = sum(u==2,2);    % 81 x 1, how many 2's?
n3 = sum(u==3,2);    % 81 x 1, how many 3's?
% id is 81 x 15, each row contains one and only one "true" value
% at column c(r). This c(r) is precisely the index of the packed
% HOT representation combination r corresponds to:
id  = ( n1==(pwrs(:,1))' ) & ( n2==(pwrs(:,2))' ) & ( n3==(pwrs(:,3))' );
% Use find to retrieve this index:
[r,c] = find(id);
% Note find will not return r,c in the order we need:
[~,idx] = sort(r);
idx     = c(idx);
end

% -------------------------------------------------------------------------
% Special functions
function g1 = tabesh_G1(l1,l2,l3)
th    = sqrt(eps(1.0e-3))*10;
% This is not managed in the original paper, but seems
% necessary since the asymptotic expressions become
% singular if l3->0:
sing1 = ( (abs(l1-l2)<th) & (abs(l2-l3)<th) );
sing2 = ( (abs(l2-l3)<th) & ~sing1 );
reg   = ( ~sing1 & ~sing2 );
g1    = zeros(size(l1));
% ----
g1(sing1) = 3/8;
% ----
g1(sing2) = (l1(sing2)+2*l3(sing2));
g1(sing2) = g1(sing2).*g1(sing2)/24;
g1(sing2) = g1(sing2)  - (l2(sing2)-l3(sing2)).*((l1(sing2)+2*l3(sing2)).*(5*l1(sing2)+4*l3(sing2)))./l3(sing2)/72;
g1(sing2) = g1(sing2)./(l3(sing2).*l3(sing2));
% ----
g1(reg) = (l1(reg)+l2(reg)+l3(reg))./(l2(reg)-l3(reg));
g1(reg) = g1(reg).*g1(reg)./l2(reg)/18;
g1(reg) = g1(reg).*( 2*l2(reg) + (l3(reg).*l3(reg)-3*l2(reg).*l3(reg))./sqrt(l2(reg).*l3(reg)) );
end

function g2 = tabesh_G2(l1,l2,l3)
th   = sqrt(eps(1.0e-3))*10;
% This is not managed in the original paper, but seems
% necessary since the asymptotic expressions become
% singular if l3->0:
sing1 = ( (abs(l1-l2)<th) & (abs(l2-l3)<th) );
sing2 = ( (abs(l2-l3)<th) & ~sing1 );
reg   = ( ~sing1 & ~sing2 );
g2    = zeros(size(l1));
% ----
g2(sing1) = 3/4;
% ----
g2(sing2) = (l1(sing2)+2*l3(sing2));
g2(sing2) = g2(sing2).*g2(sing2)/12;
g2(sing2) = g2(sing2) + (l2(sing2)-l3(sing2)).*(l1(sing2)+2*l3(sing2))/6;
g2(sing2) = g2(sing2)./(l3(sing2).*l3(sing2));
% ----
g2(reg) = (l1(reg)+l2(reg)+l3(reg))./(l2(reg)-l3(reg));
g2(reg) = g2(reg).*g2(reg)/3;
g2(reg) = g2(reg).*( (l2(reg)+l3(reg))./sqrt(l2(reg).*l3(reg)) - 2 );
end

function h1 = tabesh_H1(l1,l2,l3)
th   = sqrt(eps(1.0e-3))*10;
% This is not managed in the original paper, but seems
% necessary
sing = ( (abs(l1-l2)<th) & (abs(l2-l3)<th) );
h1   = zeros(size(l1));
h1(~sing) = (l1(~sing)+l2(~sing)+l3(~sing))./(3*l1(~sing));
h1(~sing) = h1(~sing).*h1(~sing);
h1(sing)  = 1;
end

function f1 = tabesh_F1_core(l1,l2,l3,mt)
f1 = (3*l1.*l1-l1.*l2-l1.*l3-l2.*l3)./(3*l1.*sqrt(l2.*l3));
f1 = f1.*carlsonRd(l1./l2,l1./l3,1,mt);
f1 = f1 + sqrt(l2.*l3)./l1.*carlsonRf(l1./l2,l1./l3,1,mt) - 1;
f1 = f1.*(l1+l2+l3).*(l1+l2+l3)./((l1-l2).*(l1-l3)*18);
end

function f2 = tabesh_F2_core(l1,l2,l3,mt)
f2 = (2*l1-l2-l3)./(3*sqrt(l2.*l3));
f2 = f2.*carlsonRd(l1./l2,l1./l3,1,mt);
f2 = f2 + (l2+l3)./sqrt(l2.*l3).*carlsonRf(l1./l2,l1./l3,1,mt) - 2;
f2 = f2.*(l1+l2+l3).*(l1+l2+l3)./((l2-l3).*(l2-l3)*3);
end

function f1 = tabesh_F1(l1,l2,l3,mt)
th    = sqrt(eps(1.0e-3))*10;
sing1 = ( (abs(l1-l2)<th) & (abs(l1-l3)<th) );
sing2 = ( ~sing1 & (abs(l1-l2)<th) );
sing3 = ( ~sing1 & (abs(l1-l3)<th) );
reg   = ( ~sing1 & ~sing2 & ~sing3 );
f1    = zeros(size(l1));
% In the first case, use constant values:
f1(sing1) = 1/5;
% In the second case, A_1111(l1,l1,l3) = 3*A_2233(l3,l1,l1)
% so that F1(l1,l1,l3) = 3*(1/6)*F2(l3,l1,l1). Use this to 
% compute the value at exactly l1=l2
x0 = tabesh_F2( l3(sing2), l1(sing2), l1(sing2), mt );
% ---
% Compute the value at l1 = l2 + 2*th:
x1 = tabesh_F1_core( l2(sing2)+2*th, l2(sing2), l3(sing2), mt );
% ---
% Compute the final value by interpolation:
f1(sing2) = x0 + (l1(sing2)-l2(sing2)).*(x1-x0)/th/2;
% In the third case, A_1111(l1,l2,l1) = 3*A_2233(l2,l1,l1)
% so that F1(l1,l2,l1) = 3*(1/6)*F2(l2,l1,l1). Use this to 
% compute the value at exactly l1=l2
x0 = tabesh_F2( l2(sing3), l1(sing3), l1(sing3), mt );
% ---
% Compute the value at l1 = l3 + 2*th:
x1 = tabesh_F1_core( l3(sing3)+2*th, l2(sing3), l3(sing3), mt );
% ---
% Compute the final value by interpolation:
f1(sing3) = x0 + (l1(sing3)-l3(sing3)).*(x1-x0)/th/2;
% In the fourth case, use the default implementation
f1(reg) = tabesh_F1_core( l1(reg), l2(reg), l3(reg), mt );
end

function f2 = tabesh_F2(l1,l2,l3,mt)
% Note in this case we cannot assume l1>=l2>=l3,
% since tabesh_F1 calls this function with exchanged
% roles for l1 and l3.
th    = sqrt(eps(1.0e-3))*10;
sing1 = ( (abs(l1-l2)<th) & (abs(l2-l3)<th) );
sing2 = ( ~sing1 & (abs(l2-l3)<th) );
reg   = ( ~sing1 & ~sing2 );
f2    = zeros(size(l1));
% In the first case, use constant values:
f2(sing1) = 2/5;
% In the second case, use eq. [A10] from Tabesh's paper to
% compute the value at l2=l3:
% ---
t1 = l1(sing2);
t3 = l3(sing2);
x1 = (t1+2*t3)./(t1-t3)./t3/12;
x1 = x1.*x1;
x2 = t3.*(t1+2*t3) + t1.*(t1-4*t3).*tabesh_alpha(1-t1./t3);
x0 = x1.*x2;
% ---
% Compute the value at l2 = l3 + 2*th:
x1 = tabesh_F2_core( t1, t3+2*th, t3, mt );
% ---
% Compute the final value by interpolation:
f2(sing2) = x0 + (l2(sing2)-l3(sing2)).*(x1-x0)/th/2;
% In the third case, use the default implementation:
f2(reg) = tabesh_F2_core( l1(reg), l2(reg), l3(reg), mt );
end

function a = tabesh_alpha(x)
p     = ( abs(x)> 10*eps );
a     = zeros(size(x));
a(p)  = atan(sqrt(abs(x(p))))./sqrt(abs(x(p)));
a(~p) = 1 - abs(x(~p))/3;
end
