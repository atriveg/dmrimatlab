function [lpar,lperp,f,nit] = compute_lambda_perp_and_f_from_E_chunk(E,b,fmin,fmax,f0,mu,nu,mlp,lpar,nmax,dl,dC,ADC0)
% function [lpar,lperp,f,nit] = compute_lambda_perp_and_f_from_E_chunk(E,b,fmin,fmax,f0,mu,nu,mlp,lpar,nmax,dl,dC,ADC0)
%
%  Computes the components of the micro-structural model given the
%  attenuation signal. This is a helper function to atti2micro. This use
%  case is meant for the case where only two shells are available, so that
%  only two parameters can be inferred: lpar will be kept constant, and
%  the optimization will work over just lperp and f.
%
%    E: N x NB, the signal in the natural domain.
%    b: NB x 1, the set of b-values for each channel in the last dimension
%       of E
%    fmin: X x 1, the minimum admissible value of the partial volume 
%       fraction of free water
%    fmax: N x 1, the minimum admissible value of the partial volume 
%       fraction of free water
%    f0: N x 1, the initial temptative value of the partial volume fraction
%       of free water
%    mu: 1 x 1, penalty weight to avoid lperp being extremely small
%    nu: 1 x 1, penalty weight to avoid lperp being extremely close to lpar
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
%    ADC0: 1x1, free water diffusivity at human body (default: 3.0e-3)
%    mlp: 1x1, minimum allowed value for lperp (default: 0.01e-3)
%
%    lpar: N1 x N2 x N3 x ... x Nd, the diffusivity in the parallel
%       direction. NOTE: this is the same value passed to the function,
%       since no optimization is performed over lpar.
%    lperp: N1 x N2 x N3 x ... x Nd, the diffusivity in the perpendicular
%       direction.
%    f: N x 1
%    nit: N x 1, the number of iterations needed to
%       obtain the solution
%

%%% =======================================================================
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
%%% =======================================================================
% Check sizes:
M  = size(E,1);
%%% =======================================================================
% Find a first iteration
f     = min(max(f0,fmin),fmax);
Ec    = correct_E_with_f(E,f,b,ADC0);
[lpar,lperp,nit] = ...
    compute_lambda_perp_from_E(Ec,b,mu,nu,mlp,lpar,nmax,dl,dC);
%%% =======================================================================
% Normalize b to avoid numerical issues (this scales diffusivities 1000
% times)
b  = b(:)/1000; % NB x 1
%%% =======================================================================
% Normalize diffusivity-related parameters as well:
mlp   = 1000*mlp;
ADC0  = 1000*ADC0;
lpar  = 1000*lpar;
lperp = 1000*lperp;
%%% =======================================================================
% Make sure the solution is physically meaningful by just projecting the
% temptative solutions onto the feasible region:
[bnd1,bnd2,bnd4,bnd5,ins,lperp,f] = ...
    project_on_boundaries(lperp,f,mlp,lpar,fmin,fmax);
%%% =======================================================================
% Constrained Newton-Raphson iterations
doIt = true(M,1);
tau  = 0.001*ones(M,1);
n    = 0;
while( any(doIt) && (n<nmax) )
    % ---------------------------------------------------------------------
    % Process only the necessary data:
    E_     = E(doIt,:);    % P x NB
    lperp_ = lperp(doIt);  % P x 1
    lpar_  = lpar(doIt);   % P x 1
    delta_ = lpar_-lperp_; % P x 1
    f_     = f(doIt);      % P x 1
    fmin_  = fmin(doIt);   % P x 1
    fmax_  = fmax(doIt);   % P x 1
    tau_   = tau(doIt);    % P x 1
    bnd1_  = bnd1(doIt);   % P x 1
    bnd2_  = bnd2(doIt);   % P x 1
    bnd4_  = bnd4(doIt);   % P x 1
    bnd5_  = bnd5(doIt);   % P x 1
    ins_   = ins(doIt);    % P x 1
    % ---------------------------------------------------------------------
    % In general (i.e. for points inside the feasible region), we aim at
    % minimizing a cost function given by:
    %   Q = (1/2) sum_i ( El(bi,f) + lp*bi + F(bi*d) )^2 + mu*P(lp,d)
    %     = (1/2) sum_i c_i(f,lp,d)^2 + mu*P(lp,d)
    % in the three variables f, lp (lperp_), and d (delta_=lpar_-lperp_), 
    % where c_i(f,lp,d) = El(bi,f) + lp*bi + F(bi*d) and function F stands
    % for the log-sqrt-erf decay. The minimum is reached
    %  when the gradient equals 0, i.e.:
    %   Q_lp = sum_i c_i(f,lp,d)*bi          + mu*P_x(lp,d) = 0
    %   Q_d  = sum_i c_i(f,lp,d)*F'(bi*d)*bi + mu*P_y(lp,d) = 0
    %   Q_f  = sum_i c_i(f,lp,d)*El_f(bi,f)                 = 0
    % or, in matrix form:
    %   [ El_f(bi,f), bi, F'(d*bi)*bi ]^T*[ c_i(f,lp,d) ]
    %                                + mu*[ 0, P_x(lp,d), P_y(lp,d) ]^T = 0
    % By taking first order Taylor series expansions around the current
    % estimate of lp and d for the terms c_i, P_x, and P_y, we get a LLS
    % problem with closed form solution to obtain a new iteration
    %   [ El_f(bi,f0), bi, F'(d0*bi)*bi ]^T*[ c_i(f0,lp0,d0) ]
    %      + [ El_f(bi,f0), bi, F'(d0*bi)*bi  ]^T 
    %                  * [ El_f(bi,f0), bi, F'(d0*bi)*bi ]*[ u1, u2, u3 ]^T
    %      + mu*[ 0, P_x(lp0,d0), P_y(lp0,d0) ]^T
    %      + mu*[ 0,            0,            0 ]
    %           [ 0, P_xx(lp0,d0), P_xy(lp0,d0) ]
    %           [ 0, P_xy(lp0,d0), P_yy(lp0,d0) ]*[ u1, u2, u3 ]^T = 0
    % for u1 = (f-f0), u2 = (lp-lp0), and u3 = (d-d0). Hence:
    %   ( Jc^t*Jc + mu*Hp )*[ u1, u2 ]^T = - ( Jc^T*[ ci ] + mu*Jp^T )
    % When mu=0 (no penalty applies) this completely equivalent to make
    % Newton-Raphson iterarions for the non-linear problem:
    %   [ c_i ] = [ 0 ]. 
    % To avoid numerical issues, any matrix inversion is performed over a
    % regularized version: a term of the form tau*I_N is added, so that
    % if the iterations fail tau is increased and the method becomes alike
    % to gradient-descent. If the iterations succeed, tau is decreased so
    % that the method becomes alike to pure Newton-Raphson iterations and
    % its convergence is faster.
    % ---------------------------------------------------------------------
    % Compute the current cost function:
    Ec_    = correct_E_with_f(E_,f_,b,ADC0); % P x NB
    El_    = log(Ec_);                       % P x NB
    argF   = delta_*(b');                    % P x NB
    func_  = sqrt_over_erfL(argF);           % P x NB, the log-sqrt-erf decay
    cost_  = lperp_*(b') + func_ + El_;      % P x NB
    % ---------------------------------------------------------------------
    % Compute the current jacobian:
    d1_ = correct_E_with_f_d1(E_,f_,b,ADC0);   % P x NB
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
    [PV_,Px_,Pxx_] = compute_penalty_term(lperp_,delta_,mu,nu);
    % ---------------------------------------------------------------------
    % Compute the jacobian:
    J11 = sum(d1_.*d1_,2);        % P x 1
    J22 = sum(d2_.*d2_,2) + Pxx_; % P x 1
    LM  = tau_.*(J11+J22)/2;      % P x 1
    J11 = J11 + LM;               % P x 1
    J22 = J22 + LM;               % P x 1
    J12 = d1_*b;                  % P x 1
    J12 = J12 - sum(d1_.*d3_,2);  % P x 1
    % ---------------------------------------------------------------------
    % Invert:
    Jdet = J11.*J22 - J12.*J12; % P x 1
    Ji11 = J22./Jdet;           % P x 1
    Ji22 = J11./Jdet;           % P x 1
    Ji12 = -J12./Jdet;          % P x 1
    % ---------------------------------------------------------------------
    % Compute the step:
    f1 = sum(d1_.*cost_,2);      % P x 1
    f2 = cost_*b + Px_;          % P x 1
    f2 = f2 - sum(d3_.*cost_,2); % P x 1 
    u1 = -Ji11.*f1 - Ji12.*f2;   % P x 1
    u2 = -Ji12.*f1 - Ji22.*f2;   % P x 1
    % ---------------------------------------------------------------------
    % Manage boundary points. There are four boundaries for f<fmin, f>fmax,
    % lperp<mlp, and lperp>lpar. At each of them:
    %   - If the step points inside the feasible region, we do just
    %     nothing, and flag the voxel as "inside" since the step is taking
    %     it out from the boundary.
    %   - If the step points outside the feasible region, we design a new
    %     step with a Newton-Raphson iteration constrained to the
    %     corresponding boundary.
    % ------------ First boundary: lperp = mlp
    if( any(bnd1_) )
        % ---------
        pi_ = ( bnd1_ & (u2>0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd1_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd1_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        % Optimize the problem in one variable by fixing lp = mlp
        if(any(po_))
            J11     = sum(d1_(po_,:).*d1_(po_,:),2);   % Pb1o x 1
            f1      = sum(d1_(po_,:).*cost_(po_,:),2); % Pb1o x 1
            u1(po_) = -f1./J11./(1+tau(po_)/2);        % Pb2o x 1
            u2(po_) = 0;                               % Pb1o x 1
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb1 = po_;               % P x 1
        % ---------
    else
        forceb1 = false(size(bnd1_)); % P x 1
    end
    % ------------ Second boundary: lp = lpar
    if( any(bnd2_) )
        % ---------
        pi_ = ( bnd2_ & (u2<0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd2_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd2_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        % Optimize the problem in one variable by fixing lperp = lpar
        if(any(po_))
            J11     = sum(d1_(po_,:).*d1_(po_,:),2);   % Pb2o x 1
            f1      = sum(d1_(po_,:).*cost_(po_,:),2); % Pb2o x 1
            u1(po_) = -f1./J11./(1+tau(po_)/2);        % Pb2o x 1
            u2(po_) = 0;
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb2 = po_;                % P x 1
        % ---------
    else
        forceb2 = false(size(bnd2_)); % P x 1
    end
    % ------------ Fourth boundary: f < fmin
    if( any(bnd4_) )
        % ---------
        pi_ = ( bnd4_ & (u1>0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd4_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd4_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        if(any(po_))
            % ----------------
            % In this case, we just fix the value of f and call the
            % optimization method for fixed f!
            Ec4 = correct_E_with_f(E_(po_,:),f_(po_),b,ADC0); % P x NB
            [~,lpeb4,~] = ...
                compute_lambda_perp_from_E( Ec4,b*1000,mu,nu,mlp, ...
                lpar_(po_,:)/1000, ...
                nmax,dl,dC );
            % ----------------
            % But re-normalization is required since lpar and lperp are
            % returned in their natural units
            lpeb4 = 1000*lpeb4;
            % ----------------
            % Assign step:
            u1(po_) = 0;
            u2(po_) = lpeb4  - lperp_(po_);
            % ----------------
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb4 = po_;                % P x 1
        % ---------
    else
        forceb4 = false(size(bnd4_)); % P x 1
    end
    % ------------ Fifth boundary: f > fmax
    if( any(bnd5_) )
        % ---------
        pi_ = ( bnd5_ & (u1<0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd5_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd5_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        if(any(po_))
            % ----------------
            % In this case, we just fix the value of f and call the
            % optimization method for fixed f!
            Ec5 = correct_E_with_f(E_(po_,:),f_(po_),b,ADC0); % P x NB
            [~,lpeb5,~] = ...
                compute_lambda_perp_from_E( Ec5,b*1000,mu,nu,mlp, ...
                lpar_(po_,:)/1000, ...
                nmax,dl,dC );
            % ----------------
            % But re-normalization is required since lpar and lperp are
            % returned in their natural units
            lpeb5 = 1000*lpeb5;
            % ----------------
            % Assign step:
            u1(po_) = 0;
            u2(po_) = lpeb5  - lperp_(po_);
            % ----------------
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb5 = po_;                % P x 1
        % ---------
    else
        forceb5 = false(size(bnd5_)); % P x 1
    end
    % ---------------------------------------------------------------------
    % Compute the new parameters:
    fN     = f_     + u1;    % P x 1
    lperpN = lperp_ + u2;    % P x 1
    % ---------------------------------------------------------------------
    % Recheck if the solutions are in bounds or otherwise they are upon the
    % frontiers of the feasible region:
    [bnd1N,bnd2N,bnd4N,bnd5N,~,lperpN,fN] = ...
        project_on_boundaries(lperpN,fN,mlp,lpar_,fmin_,fmax_);
    deltaN = lpar_ - lperpN; % P x 1
    % But we have to make sure that those values processed as boundary
    % points are still boundary points
    % ---
    bnd1N(forceb1) = true;
    % ---
    bnd2N(forceb2) = true;
    bnd1N(bnd2N)   = false;
    % ---
    bnd4N(forceb4) = true;
    bnd1N(bnd4N)   = false;
    bnd2N(bnd4N)   = false;
    % ---
    bnd5N(forceb5) = true;
    bnd1N(bnd5N)   = false;
    bnd2N(bnd5N)   = false;
    bnd4N(bnd5N)   = false;
    % ---
    insN = ~(bnd1N|bnd2N|bnd4N|bnd5N);
    % ---------------------------------------------------------------------
    % Compute the new cost...
    Ec_    = correct_E_with_f(E_,fN,b,ADC0); % P x NB
    El_    = log(Ec_);                       % P x NB
    costN  = lperpN*(b') + sqrt_over_erfL(deltaN*(b')) + El_; % P x NB
    % ... and the new penalty:
    PVN    = compute_penalty_term(lperpN,deltaN,mu,nu);       % P x 1
    % ---------------------------------------------------------------------
    % Check if the iteration actually succeeded:
    cost_         = sum(cost_.*cost_,2)/2 + PV_; % P x 1
    costN         = sum(costN.*costN,2)/2 + PVN; % P x 1
    fails         = costN>cost_;
    tau_(fails)   = 10*tau_(fails);
    tau_(~fails)  = tau_(~fails)/10;
    fN(fails)     = f_(fails);
    lperpN(fails) = lperp_(fails);
    bnd1N(fails)  = bnd1_(fails);
    bnd2N(fails)  = bnd2_(fails);
    bnd4N(fails)  = bnd4_(fails);
    bnd5N(fails)  = bnd5_(fails);
    insN(fails)   = ins_(fails);
    % ---------------------------------------------------------------------
    % Prepare for the new step
    f(doIt)     = fN;
    lperp(doIt) = lperpN;
    bnd1(doIt)  = bnd1N;
    bnd2(doIt)  = bnd2N;
    bnd4(doIt)  = bnd4N;
    bnd5(doIt)  = bnd5N;
    ins(doIt)   = insN;
    tau(doIt)   = tau_;
    nit(doIt)   = nit(doIt) + 1;
    doIt(doIt)  = (abs(cost_-costN)>(dC*dC) ) | ...
        (abs(f_-fN)>dl*333) | ...
        (abs(lperp_-lperpN)>dl*1000) | fails;
    n           = n+1;
    % ---------------------------------------------------------------------
end
%%% =======================================================================
% Scale back:
lpar  = lpar/1000;
lperp = lperp/1000;

%%% =======================================================================
% Compute the corrected logarithmic signal for a given value of the partial
% volume fraction of confined water. It assumes f is within the allowed
% range, check it from the outside:
function E = correct_E_with_f(E,f,bs,ADC0)
is_broadcast_available = is_broadcast_available_test;
th = sqrt(eps);
S0 = exp(-(bs')*ADC0); % 1 x NB
E  = E - (1-f)*S0;     % M x NB
if(is_broadcast_available)
    E = E./f;                     % M x NB
else
    E = bsxfun(@(x,y)(x./y),E,f); % M x NB
end
E(isinf(E)) = 1;
E(E>1)      = 1;
E(E<th)     = th;

%%% =======================================================================
% Compute the derivative of the corrected logarithmic signal for a given 
% value of the partial volume fraction of confined water with respect to 
% such partial volume fraction. It assumes f is within the allowed range,
% check it from the outside:
function d1 = correct_E_with_f_d1(E,f,bs,ADC0)
is_broadcast_available = is_broadcast_available_test;
th = sqrt(eps);
S0 = exp(-(bs')*ADC0); % 1 x NB
E  = E - (1-f)*S0;     % M x NB
E(E<th) = th;
f(f<th) = th;
d1 = 1./E;             % M x NB
if(is_broadcast_available)
    d1 = d1.*S0;                      % M x NB
    d1 = d1 - 1./f;                   % M x NB
else
    d1 = bsxfun(@(x,y)(x.*y),d1,S0);  % M x NB
    d1 = bsxfun(@(x,y)(x-y),d1,1./f); % M x NB
end

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
function [bnd1,bnd2,bnd4,bnd5,ins,lperp,f] = ...
    project_on_boundaries(lperp,f,mlp,lpar,fmin,fmax)
% -------
bnd1        = (lperp<mlp);        % M x 1, first frontier
lperp(bnd1) = mlp;
% -------
bnd2        = (lperp>lpar);       % M x 1, second frontier
lperp(bnd2) = lpar(bnd2);
bnd1(bnd2)  = false;
% -------
bnd4        = (f<fmin);
f(bnd4)     = fmin(bnd4);
bnd1(bnd4)  = false;
bnd2(bnd4)  = false;
% -------
bnd5        = (f>fmax);
f(bnd5)     = fmax(bnd5);
bnd1(bnd5)  = false;
bnd2(bnd5)  = false;
bnd4(bnd5)  = false;
% -------
ins         = ~(bnd1|bnd2|bnd4|bnd5);
% -------

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
