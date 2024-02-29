function [lpar,lperp,nit] = compute_lambda_perp_from_E(E,b,mu,nu,mlp,lpar,nmax,dl,dC,chunksz)
% function [lpar,lperp,nit] = compute_lambda_perp_from_E(E,b,mu,nu,mlp,lpar,nmax,dl,dC,chunksz)
%
%  Computes the parallel diffusivity in a simplified scenario where the
%  perpendicular diffusivity is known (or fixed) beforehand. This means we
%  just need to invert the log-sqrt-over-erf function.
%
%    E: N1 x N2 x N3 x ... x Nd x NB, the signal in the natural domain.
%    b: NB x 1, the set of b-values for each channel in the last dimension
%       of E
%    mu: 1 x 1, penalty weight to avoid lperp being extremely small
%      (default 0.001)
%    nu: 1 x 1, penalty weight to avoid lperp being extremely close to lpar
%      (default 0.001)
%    mlp: 1 x 1, minimum allowed value for lperp (default: 0.01e-3)
%    lpar: N x 1, non-optimized value of lpar, which
%       can be different at each voxel. NOTE: no sanity checks will be 
%       performed on lpar to ensure it lies within a physically meaningful
%       range (default: all 2.0e-3)
%    nmax: 1 x 1, maximum number of Newton-Raphson iterations 
%       (default: 100).
%    dl: maximum allowed chage in lambda from one iteration to the next
%       (default: 1.0e-6).
%    dC: maximum allowed change in the (squared) norm of the vector
%       objective funcion (default: 1.0e-3)
%    chunksz: 1x1, data are processed chunk-by-chunk (default: 100)
%
%    lpar: N1 x N2 x N3 x ... x Nd, the diffusivity in the parallel
%       direction.
%    lperp: N1 x N2 x N3 x ... x Nd, the diffusivity in the perpendicular
%       direction.
%    nit: N1 x N2 x N3 x ... x Nd, the number of iterations needed to
%       obtain the solution
%

%%% =======================================================================
% Reshape to work with arbitrary dimensions:
sz   = size(E);
NB   = sz(end);
sz2  = sz(1:end-1);
M    = prod(sz2);
E    = reshape(E,[M,NB]); % M x NB
%%% =======================================================================
if( nargin<3 )
    mu = 0.001;
end
if( nargin<4 )
    nu = 0.001;
end
if( nargin<5 )
    mlp = 0.01e-3;
end
if( nargin<6 )
    lpar = 2.0e-3*ones(M,1);
else
    if(isscalar(lpar))
        lpar = lpar*ones(M,1);
    end
end
if( nargin<7 )
    nmax = 100;
end
if( nargin<8 )
    dl = 1.0e-6;
end
if( nargin<9 )
    dC = 1.0e-3;
end
if( nargin<10 )
    chunksz = 10000;
end
%%% =======================================================================
lperp = zeros(M,1);
nit   = zeros(M,1);
%%% =======================================================================
for ck=1:ceil(M/chunksz)
    % -----------
    idi       = (ck-1)*chunksz+1;
    idf       = min(ck*chunksz,M);
    % -----------
    [lpec,nitc] = compute_lambda_perp_from_E_chunk( ...
        E(idi:idf,:), ...
        b, ...
        mu, ...
        nu, ...
        mlp, ...
        lpar(idi:idf,:), ...
        nmax, ...
        dl, ...
        dC );
    % -----------
    lperp(idi:idf) = lpec;
    nit(idi:idf)   = nitc;
    % -----------
end
%%% =======================================================================
% Reshape and scale back:
if(numel(sz2)>1)
    lpar  = reshape(lpar,sz2);
    lperp = reshape(lperp,sz2);
    nit   = reshape(nit,sz2);
end


%%% =======================================================================
function [lpe,nit] = compute_lambda_perp_from_E_chunk( E, b, mu, nu, mlp, lpa, nmax, dl, dC )
%%% =======================================================================
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
%%% =======================================================================
% Normalize b to avoid numerical issues (this scales diffusivities 1000
% times)
b    = b(:)/1000; % NB x 1
lpa  = 1000*lpa;
%%% =======================================================================
% Check sizes and compute the log-signal:
M  = size(E,1);
E(E<eps)   = eps;
E(E>1-eps) = 1-eps;
El = log(E); % M x NB
%%% =======================================================================
% Find a first iteration via greedy search
ND    = 5;                   % 1 x 1
NL    = 4;                   % 1 x 1
lpem  = mlp*ones(M,1);       % M x 1
lpeM  = lpa;                 % M x 1
lstep = (lpeM-lpem)/(ND-1);  % M x 1
cost  = zeros(M,ND);         % M x ND
for nl=1:NL
    for nd=1:ND
        % -----------------------------------------------------------------
        lpe   = lpem + lstep*(nd-1);                        % M x 1
        delta = lpa - lpe;
        cost_ = El + lpe*(b') + sqrt_over_erfL(delta*(b')); % M x NB
        cost_ = sum(cost_.*cost_,2);                        % M x 1
        cost(:,nd) = cost_;
        % -----------------------------------------------------------------  
    end
    % -----------------------------------------------------------------
    [~,idx] = min( cost, [], 2 ); % idx has size M x 1
    lpe     = lpem + lstep.*(idx-1);
    % -----------------------------------------------------------------
    lpem  = max(mlp,lpe-lstep/2); % M x 1
    lpeM  = min(lpa,lpe+lstep/2); % M x 1
    lstep = (lpeM-lpem)/(ND-1);   % M x 1
end
%%% =======================================================================
% Constrained Newton-Raphson iterations
doIt  = true(M,1);
tau   = 0.001*ones(M,1);
nit   = zeros(M,1);
n     = 0;
bnd1  = false(size(lpa));
bnd2  = false(size(lpa));
ins   = true(size(lpa));
while( any(doIt) && (n<nmax) )
    % ---------------------------------------------------------------------
    % Process only the necessary data:
    El_    = El(doIt,:);   % P x NB
    lpe_   = lpe(doIt);    % P x 1
    lpa_   = lpa(doIt);    % P x 1
    delta_ = lpa_ - lpe_;  % P x 1
    tau_   = tau(doIt);    % P x 1
    bnd1_  = bnd1(doIt);   % P x 1
    bnd2_  = bnd2(doIt);   % P x 1
    ins_   = ins(doIt);    % P x 1
    % ---------------------------------------------------------------------
    % Compute the current cost function:
    argF   = delta_*(b');             % P x NB
    func_  = sqrt_over_erfL(argF);    % P x NB, the log-sqrt-erf decay
    cost_  = lpe_*(b') + func_ + El_; % P x NB
    % ---------------------------------------------------------------------
    % Compute the current jacobian:
    d2_ = b';                                  % 1 x NB
    d3_ = sqrt_over_erfL_d1(argF);             % P x NB
    if(is_broadcast_available)
        d3_ = d3_.*(b');                       % P x NB
        d2_ = d2_ - d3_;                       % P x NB
    else
        d3_ = bsxfun( @(x,y)(x.*y), d3_, b' ); % P x NB
        d2_ = bsxfun( @(x,y)(x-y), d2_, d3_ ); % P x NB
    end
    % ---------------------------------------------------------------------
    % Compute the current penalty term:
    [PV_,Px_,Pxx_] = compute_penalty_term(lpe_,delta_,mu,nu);
    % ---------------------------------------------------------------------
    J22  = sum(d2_.*d2_,2) + Pxx_; % P x 1
    LM   = tau_.*J22;              % P x 1
    J22  = J22 + LM;               % P x 1
    % ---------------------------------------------------------------------
    % Invert:
    Ji22 = 1./(J22+LM);            % P x 1
    % ---------------------------------------------------------------------
    % Compute the step:
    f2 = cost_*b + Px_;          % P x 1
    f2 = f2 - sum(d3_.*cost_,2); % P x 1
    u2 = -Ji22.*f2;              % P x 1
    % ---------------------------------------------------------------------
    % Manage boundary points.
    % ---------------------------------------------------------------------
    % What if lpe = mlp?
    if( any(bnd1_) )
        % ---------
        pi_ = ( bnd1_ & (u2>0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd1_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd1_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        % This means the gradient points outside the feasible region, and
        % will always do. No other way than keeping the value and exiting
        forceb1 = po_;                % P x 1
        % ---------
    else
        forceb1 = false(size(bnd1_)); % P x 1
    end
    % ---------------------------------------------------------------------
    % What if lpe = lpa?
    if( any(bnd2_) )
        % ---------
        pi_ = ( bnd2_ & (u2<0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd2_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd2_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        % This means the gradient points outside the feasible region, and
        % will always do. No other way than keeping the value and exiting
        forceb2 = po_;                % P x 1
        % ---------
    else
        forceb2 = false(size(bnd2_)); % P x 1
    end
    % ---------------------------------------------------------------------
    % Compute the new parameter:
    lpeN = lpe_ + u2;
    % ---------------------------------------------------------------------
    % Recheck if the solutions are in bounds or otherwise they are upon the
    % frontiers of the feasible region:
    bnd1N = (lpeN<mlp);
    bnd2N = (lpeN>lpa_);
    bnd1N(bnd2N) = false;
    lpeN(bnd1N) = mlp;
    lpeN(bnd2N) = lpa_(bnd2N);
    insN = ~(bnd1N|bnd2N);
    deltaN = lpa_ - lpeN;
    % ---------------------------------------------------------------------
    % Compute the new cost...
    costN  = lpeN*(b') + sqrt_over_erfL(deltaN*(b')) + El_; % P x NB
    % ... and the new penalty:
    PVN    = compute_penalty_term(lpeN,deltaN,mu,nu);       % P x 1
    % ---------------------------------------------------------------------
    % Check if the iteration actually succeeded:
    cost_         = sum(cost_.*cost_,2)/2 + PV_; % P x 1
    costN         = sum(costN.*costN,2)/2 + PVN; % P x 1
    fails         = costN>cost_;
    tau_(fails)   = 10*tau_(fails);
    tau_(~fails)  = tau_(~fails)/10;
    lpeN(fails)   = lpe_(fails);
    bnd1N(fails)  = bnd1_(fails);
    bnd2N(fails)  = bnd2_(fails);
    insN(fails)   = ins_(fails);
    % ---------------------------------------------------------------------
    % Prepare for the new step
    lpe(doIt)   = lpeN;
    bnd1(doIt)  = bnd1N;
    bnd2(doIt)  = bnd2N;
    ins(doIt)   = insN;
    tau(doIt)   = tau_;
    nit(doIt)   = nit(doIt) + 1;
    doIt(doIt)  = (   (abs(cost_-costN)>(dC*dC) ) | ...
        (abs(lpe_-lpeN)>dl*1000) | fails   ) & ...
        (~(forceb1|forceb2));
    n           = n+1;
    % ---------------------------------------------------------------------
end
lpe  = lpe/1000;
%%% =======================================================================

%%% =======================================================================
% Compute the value and the first derivatives of the function describing
% the (logarithmic) decay of the signal due to the parallel component
function f = sqrt_over_erfL(x)
xs   = real(sqrt(x));               % sqrt(x), PxNB
x2   = x.*x;                        % x^2, PxB
p    = (x<100*sqrt(eps));           % Avoid numerical issues: series expansions at x=0
% ----------------------------------------------------
% The function itself:
f    = (2/sqrt(pi))*xs./erf(xs); % P x NB, the erf-decay
f    = log(f);                   % P x NB, logartihmic decay
f(p) = x(p)/3 - x2(p)/(45/2);    % P x NB, small values

%%% =======================================================================
% Compute the value and the first derivatives of the function describing
% the (logarithmic) decay of the signal due to the parallel component
function fd = sqrt_over_erfL_d1(x)
xs   = real(sqrt(x));               % sqrt(x), PxNB
x2   = x.*x;                        % x^2, PxB
p    = (x<100*sqrt(eps));           % Avoid numerical issues: series expansions at x=0
% ----------------------------------------------------
aux = exp(-x)./(sqrt(pi).*xs.*erf(xs));    % PxNB
% First derivative:
fd    = 1./(2*x) - aux;                    % PxNB
fd(p) = 1/3 - x(p)/(45/4) + x2(p)/(945/8); % PxNB, small values

%%% =======================================================================
function [PV,Px,Pxx] = compute_penalty_term(lperp,delta,mu,nu)
% ---------------------------
pp        = (lperp<100*sqrt(eps));
lperp(pp) = 100*sqrt(eps);
pp        = (delta<100*sqrt(eps));
delta(pp) = 100*sqrt(eps);
lpar      = lperp+delta;
% ---------------------------
PV  = nu*( lperp./delta ) + mu*( delta./lperp );
% ---------------------------
Px  = lpar.*( nu./(delta.*delta) - mu./(lperp.*lperp) );
% ---------------------------
Pxx = 2*lpar.*( nu./(delta.*delta.*delta) + mu./(lperp.*lperp.*lperp) );
