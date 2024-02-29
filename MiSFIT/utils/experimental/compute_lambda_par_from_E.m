function [lpar,lperp,nit] = compute_lambda_par_from_E(E,b,lperp,nmax,dl,dC,ADC0,chunksz)
% function [lpar,lperp,nit] = compute_lambda_par_from_E(E,b,lperp,nmax,dl,dC,ADC0,chunksz)
%
%  Computes the parallel diffusivity in a simplified scenario where the
%  perpendicular diffusivity is known (or fixed) beforehand. This means we
%  just need to invert the log-sqrt-over-erf function.
%
%    E: N1 x N2 x N3 x ... x Nd x NB, the signal in the natural domain.
%    b: NB x 1, the set of b-values for each channel in the last dimension
%       of E
%    lperp: N1 x N2 x N3 x ... x Nd, non-optimized value of lperp, which
%       can be different at each voxel. NOTE: no sanity checks will be 
%       performed on lperp to ensure it lies within a physically meaningful
%       range (default: all zeros)
%    nmax: 1 x 1, maximum number of Newton-Raphson iterations 
%       (default: 100).
%    dl: maximum allowed chage in lambda from one iteration to the next
%       (default: 1.0e-6).
%    dC: maximum allowed change in the (squared) norm of the vector
%       objective funcion (default: 1.0e-3)
%    ADC0: 1x1, free water diffusivity at human body (default: 3.0e-3)
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
    lperp = zeros(M,1);
else
    if(isscalar(lperp))
        lperp = lperp*ones(M,1);
    end
end
if( nargin<4 )
    nmax = 100;
end
if( nargin<5 )
    dl = 1.0e-6;
end
if( nargin<6 )
    dC = 1.0e-3;
end
if( nargin<7 )
    ADC0 = 3.0e-3;
end
if( nargin<8 )
    chunksz = 10000;
end
%%% =======================================================================
lpar  = zeros(M,1);
nit   = zeros(M,1);
%%% =======================================================================
for ck=1:ceil(M/chunksz)
    % -----------
    idi       = (ck-1)*chunksz+1;
    idf       = min(ck*chunksz,M);
    % -----------
    [lpac,nitc] = compute_lambda_par_from_E_chunk( ...
        E(idi:idf,:), ...
        b, ...
        lperp(idi:idf,:), ...
        nmax, ...
        dl, ...
        dC, ...
        ADC0 );
    % -----------
    lpar(idi:idf) = lpac;
    nit(idi:idf)  = nitc;
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
function [lpac,nit] = compute_lambda_par_from_E_chunk( E, b, lpe, nmax, dl, dC, ADC0 )
%%% =======================================================================
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
%%% =======================================================================
% Normalize b to avoid numerical issues (this scales diffusivities 1000
% times)
b    = b(:)/1000; % NB x 1
lpe  = 1000*lpe;
ADC0 = 1000*ADC0;
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
dmin  = zeros(M,1);          % M x 1
dmax  = ADC0 - lpe;          % M x 1
dstep = (dmax-dmin)/(ND-1);  % M x 1
cost  = zeros(M,ND);         % M x ND
for nl=1:NL
    for nd=1:ND
        % -----------------------------------------------------------------
        delta = dmin + dstep*(nd-1);                        % M x 1
        cost_ = El + lpe*(b') + sqrt_over_erfL(delta*(b')); % M x NB
        cost_ = sum(cost_.*cost_,2);                        % M x 1
        cost(:,nd) = cost_;
        % -----------------------------------------------------------------  
    end
    % -----------------------------------------------------------------
    [~,idx] = min( cost, [], 2 ); % idx has size M x 1
    delta   = dmin + dstep.*(idx-1);
    % -----------------------------------------------------------------
    dmin  = max(dmin,delta-dstep/2); % M x 1
    dmax  = min(dmax,delta+dstep/2); % M x 1
    dstep = (dmax-dmin)/(ND-1);      % M x 1
end
%%% =======================================================================
% Constrained Newton-Raphson iterations
doIt = true(M,1);
tau  = 0.001*ones(M,1);
nit  = zeros(M,1);
n    = 0;
lpac = lpe+delta;
bnd1 = false(size(lpac));
bnd2 = false(size(lpac));
ins  = true(size(lpac));
while( any(doIt) && (n<nmax) )
    % ---------------------------------------------------------------------
    % Process only the necessary data:
    El_    = El(doIt,:);   % P x NB
    lpe_   = lpe(doIt);    % P x 1
    delta_ = delta(doIt);  % P x 1
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
    d1_ = sqrt_over_erfL_d1(argF);  % P x NB
    if(is_broadcast_available)
        d1_ = d1_.*(b');                       % P x NB
    else
        d1_ = bsxfun( @(x,y)(x.*y), d1_, b' ); % P x NB
    end
    % ---------------------------------------------------------------------
    J11  = sum(d1_.*d1_,2); % P x 1
    LM   = tau_.*J11;       % P x 1
    Ji11 = 1./(J11+LM);     % P x 1
    % ---------------------------------------------------------------------
    % Compute the step:
    u1 = sum(d1_.*cost_,2); % P x 1
    u1 = Ji11.*u1;          % P x 1
    % ---------------------------------------------------------------------
    % Manage boundary points.
    % ---------------------------------------------------------------------
    % What if delta=0?
    if( any(bnd1_) )
        % ---------
        pi_ = ( bnd1_ & (u1>0) );  % P x 1, this takes the voxel into the feasible region
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
    % What if lpe+delta=ADC0?
    if( any(bnd2_) )
        % ---------
        pi_ = ( bnd2_ & (u1<0) );  % P x 1, this takes the voxel into the feasible region
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
    deltaN = delta_ + u1;             % P x 1
    % ---------------------------------------------------------------------
    % Recheck if the solutions are in bounds or otherwise they are upon the
    % frontiers of the feasible region:
    bnd1N = (deltaN<0);
    bnd2N = (lpe_+deltaN>ADC0);
    bnd1N(bnd2N) = false;
    deltaN(bnd1N) = 0;
    deltaN(bnd2N) = ADC0 - lpe_(bnd2N);
    insN = ~(bnd1N|bnd2N);
    % ---------------------------------------------------------------------
    % Compute the new cost:
    costN  = lpe_*(b') + sqrt_over_erfL(deltaN*(b')) + El_; % P x NB
    % ---------------------------------------------------------------------
    % Check if the iteration actually succeeded:
    cost_         = sum(cost_.*cost_,2)/2; % P x 1
    costN         = sum(costN.*costN,2)/2; % P x 1
    fails         = costN>cost_;
    tau_(fails)   = 10*tau_(fails);
    tau_(~fails)  = tau_(~fails)/10;
    deltaN(fails) = delta_(fails);
    bnd1N(fails)  = bnd1_(fails);
    bnd2N(fails)  = bnd2_(fails);
    insN(fails)   = ins_(fails);
    % ---------------------------------------------------------------------
    % Prepare for the new step
    delta(doIt) = deltaN;
    bnd1(doIt)  = bnd1N;
    bnd2(doIt)  = bnd2N;
    ins(doIt)   = insN;
    tau(doIt)   = tau_;
    nit(doIt)   = nit(doIt) + 1;
    doIt(doIt)  = (   (abs(cost_-costN)>(dC*dC) ) | ...
        (abs(delta_-deltaN)>dl*1000) | fails   ) & ...
        (~(forceb1|forceb2));
    n           = n+1;
    % ---------------------------------------------------------------------
end
lpac = lpe + delta; % M x 1
% Scale back:
lpac  = lpac/1000;
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

