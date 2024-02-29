function [lpar,lperp,f,nit] = compute_lambda_par_and_f_from_E_chunk(E,b,fmin,fmax,f0,lperp,nmax,dl,dC,ADC0)
% function [lpar,lperp,f,nit] = compute_lambda_par_and_f_from_E_chunk(E,b,fmin,fmax,f0,lperp,nmax,dl,dC,ADC0)
%
%  Computes the components of the micro-structural model given the
%  attenuation signal. This is a helper function to atti2micro. This use
%  case is meant for the case where only two shells are available, so that
%  only two parameters can be inferred: lperp will be kept constant, and
%  the optimization will work over just lpar and f.
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
%    mlp: 1x1, minimum allowed value for lperp (default: 0.01e-3)
%
%    lpar: N x 1, the diffusivity in the parallel
%       direction.
%    lperp: N1 x N2 x N3 x ... x Nd, the diffusivity in the perpendicular
%       direction. NOTE: this is the same value passed to the function,
%       since no optimization is performed over lperp.
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
    compute_lambda_par_from_E(Ec,b,lperp,nmax,dl,dC,ADC0);
%%% =======================================================================
% Normalize b to avoid numerical issues (this scales diffusivities 1000
% times)
b  = b(:)/1000; % NB x 1
%%% =======================================================================
% Normalize diffusivity-related parameters as well:
ADC0  = 1000*ADC0;
lpar  = 1000*lpar;
lperp = 1000*lperp;
delta = lpar-lperp;
%%% =======================================================================
% Make sure the solution is physically meaningful by just projecting the
% temptative solutions onto the feasible region:
[bnd2,bnd3,bnd4,bnd5,ins,lperp,delta,f] = ...
    project_on_boundaries(lperp,delta,f,ADC0,fmin,fmax);
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
    delta_ = delta(doIt);  % P x 1
    f_     = f(doIt);      % P x 1
    fmin_  = fmin(doIt);   % P x 1
    fmax_  = fmax(doIt);   % P x 1
    tau_   = tau(doIt);    % P x 1
    bnd2_  = bnd2(doIt);   % P x 1
    bnd3_  = bnd3(doIt);   % P x 1
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
    d3_ = sqrt_over_erfL_d1(argF);             % P x NB
    if(is_broadcast_available)
        d3_ = d3_.*(b');                       % P x NB
    else
        d3_ = bsxfun( @(x,y)(x.*y), d3_, b' ); % P x NB
    end
    % ---------------------------------------------------------------------
    % Compute the jacobian:
    J11 = sum(d1_.*d1_,2);    % P x 1
    J33 = sum(d3_.*d3_,2);    % P x 1
    LM  = tau_.*(J11+J33)/2;  % P x 1
    J11 = J11 + LM;           % P x 1
    J33 = J33 + LM;           % P x 1
    J13 = sum(d1_.*d3_,2);    % P x 1
    % ---------------------------------------------------------------------
    % --- Invert:
    Jdet = J11.*J33 - J13.*J13; % P x 1
    Ji11 = J33./Jdet;           % P x 1
    Ji33 = J11./Jdet;           % P x 1
    Ji13 = -J13./Jdet;          % P x 1
    % ---------------------------------------------------------------------
    % Compute the step:
    f1 = sum(d1_.*cost_,2);    % P x 1
    f3 = sum(d3_.*cost_,2);    % P x 1
    u1 = -Ji11.*f1 - Ji13.*f3; % P x 1
    u3 = -Ji13.*f1 - Ji33.*f3; % P x 1
    % ---------------------------------------------------------------------
    % Manage boundary points. There are five boundaries for f<fmin, f>fmax,
    % delta<0, and lperp+delta>ADC0. At each of them:
    %   - If the step points inside the feasible region, we do just
    %     nothing, and flag the voxel as "inside" since the step is taking
    %     it out from the boundary.
    %   - If the step points outside the feasible region, we design a new
    %     step with a Newton-Raphson iteration constrained to the
    %     corresponding boundary.
    % ------------ Second boundary: delta = 0
    if( any(bnd2_) )
        % ---------
        pi_ = ( bnd2_ & (u3>0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd2_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd2_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        % Optimize the problem in one variable by fixing d = 0, i.e. find
        % the minimum in d of:
        %   Q = (1/2) sum_i ( El(bi,f) + bi*lp + F(0) )^2 + mu*P(lp,0),
        % then:
        %   Q_f  = sum_i c_i(f,lp,0)*El_f(bi,f)        = 0
        %   Q_lp = sum_i c_i(f,lp,0)*bi + mu*P_x(lp,0) = 0
        % or, in matrix form:
        %   [ El_f(bi,f), bi ]^T*[ c_i(lp,0) ] + [ 0, mu*P_x(lp,0) ]^T = 0
        % with the series expansion at lp = lp0:
        %     [ El_f(bi,f0), bi ]^T*[ c_i(f0,lp0,0) ] 
        %   + [ 0, mu*P_x(lp0,0) ]^T
        %   + [ El_f(bi,f0), bi ]^T*[ El_f(bi,f0), bi ] * [ u1, u2 ]^T 
        %   + [ 0, mu*P_xx(lp0,0)*u2 ]^T = 0
        if(any(po_))
            J11     = sum(d1_(po_,:).*d1_(po_,:),2);   % Pb2o x 1
            f1      = sum(d1_(po_,:).*cost_(po_,:),2); % Pb2o x 1
            u1(po_) = -f1./J11./(1+tau(po_)/2);        % Pb2o x 1
            u3(po_) = 0;
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb2 = po_;                % P x 1
        % ---------
    else
        forceb2 = false(size(bnd2_)); % P x 1
    end
    % ------------ Third boundary: lperp + delta = ADC0
    if( any(bnd3_) )
        % ---------
        pi_ = ( bnd3_ & (u3<0) );    % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd3_ & (~pi_)  );   % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd3_(pi_) = false;          % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;           % P x 1
        % ---------
        % Optimize the problem in two variables by fixing lp + d = ADC0:
        %   f  = f;
        %   d  = t;
        %   lp = ADC0 - t;
        % i.e. find the minimum in d of:
        %   Q = (1/2) sum_i ( El(bi,f) + bi*(ADC0-t) + F(bi*t) )^2
        %                            + mu*P(ADC0-t,t),
        % then:
        %   dQ/df = sum_i c_i(f,ADC0-t,t)*El_f(bi,f)
        %   dQ/dt = sum_i c_i(f,ADC0-t,t)*( F'(bi*t)*bi - bi ) 
        %             - mu*P_x(ADC0-t,t) + mu*P_y(ADC0-t,t) = 0
        % or, in matrix form:
        %   [ El_f(bi,f), F'(bi*t)*bi-bi ]^T*[ c_i(f,ADC0-t,t) ]
        %             - [ 0, mu*P_x(ADC0-t,t) + mu*P_y(ADC0-t,t) ]^T = 0
        % with the series expansion at f=f0, t = t0:
        %     [ El_f(bi,f0), F'(bi*t0)*bi-bi ]^T*[ c_i(f0,ADC0-t0,t0) ]
        %   + [ 0, - mu*P_x(ADC0-t0,t0) + mu*P_y(ADC0-t0,t0) ]^T
        %   + [ El_f(bi,f0), F'(bi*t0)*bi-bi ]^T
        %                   *[ El_f(bi,f0), F'(bi*t0)*bi-bi ]*[u1,(t-t0)]^T
        %   - [ 0, mu*( -P_xx(ADC0-t0,t0) + P_xy(ADC0-t0,t0) )(t-t0) ]^T
        %   + [ 0, mu*( -P_yx(ADC0-t0,t0) + P_yy(ADC0-t0,t0) )(t-t0) ]^T
        if(any(po_))
            J11     = sum(d1_(po_,:).*d1_(po_,:),2);   % Pb2o x 1
            f1      = sum(d1_(po_,:).*cost_(po_,:),2); % Pb2o x 1
            u1(po_) = -f1./J11./(1+tau(po_)/2);        % Pb2o x 1
            u3(po_) = 0;
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb3 = po_;                % P x 1
        % ---------
    else
        forceb3 = false(size(bnd3_)); % P x 1
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
            [lpab4,~,~] = ...
                compute_lambda_par_from_E(Ec4,b*1000,lperp_(po_,:)/1000,nmax/4,dl*10,dC*10,ADC0/1000);
            % ----------------
            % But re-normalization is required since lpar and lperp are
            % returned in their natural units
            lpab4 = 1000*lpab4;
            del4  = lpab4-lperp_(po_,:);
            % ----------------
            % Assign step:
            u1(po_) = 0;
            u3(po_) = del4  - delta_(po_);
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
            [lpab5,~,~] = ...
                compute_lambda_par_from_E(Ec5,b*1000,lperp_(po_,:)/1000,nmax/4,dl*10,dC*10,ADC0/1000);
            % ----------------
            % But re-normalization is required since lpar and lperp are
            % returned in their natural units
            lpab5 = 1000*lpab5;
            del5  = lpab5-lperp_(po_,:);
            % ----------------
            % Assign step:
            u1(po_) = 0;
            u3(po_) = del5  - delta_(po_);
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
    fN     = f_     + u1; % P x 1
    lperpN = lperp_;      % P x 1
    deltaN = delta_ + u3; % P x 1
    % ---------------------------------------------------------------------
    % Recheck if the solutions are in bounds or otherwise they are upon the
    % frontiers of the feasible region:
    [bnd2N,bnd3N,bnd4N,bnd5N,~,lperpN,deltaN,fN] = ...
        project_on_boundaries(lperpN,deltaN,fN,ADC0,fmin_,fmax_);
    % But we have to make sure that those values processed as boundary
    % points are still boundary points
    % ---
    bnd2N(forceb2) = true;
    % ---
    bnd3N(forceb3) = true;
    bnd2N(bnd3N)   = false;
    % ---
    bnd4N(forceb4) = true;
    bnd2N(bnd4N)   = false;
    bnd3N(bnd4N)   = false;
    % ---
    bnd5N(forceb5) = true;
    bnd2N(bnd5N)   = false;
    bnd3N(bnd5N)   = false;
    bnd4N(bnd5N)   = false;
    % ---
    insN = ~(bnd2N|bnd3N|bnd4N|bnd5N);
    % ---------------------------------------------------------------------
    % Compute the new cost...
    Ec_    = correct_E_with_f(E_,fN,b,ADC0); % P x NB
    El_    = log(Ec_);                       % P x NB
    costN  = lperpN*(b') + sqrt_over_erfL(deltaN*(b')) + El_; % P x NB
    % ---------------------------------------------------------------------
    % Check if the iteration actually succeeded:
    cost_         = sum(cost_.*cost_,2)/2; % P x 1
    costN         = sum(costN.*costN,2)/2; % P x 1
    fails         = costN>cost_;
    tau_(fails)   = 10*tau_(fails);
    tau_(~fails)  = tau_(~fails)/10;
    fN(fails)     = f_(fails);
    deltaN(fails) = delta_(fails);
    bnd2N(fails)  = bnd2_(fails);
    bnd3N(fails)  = bnd3_(fails);
    bnd4N(fails)  = bnd4_(fails);
    bnd5N(fails)  = bnd5_(fails);
    insN(fails)   = ins_(fails);
    % ---------------------------------------------------------------------
    % Prepare for the new step
    f(doIt)     = fN;
    delta(doIt) = deltaN;
    bnd2(doIt)  = bnd2N;
    bnd3(doIt)  = bnd3N;
    bnd4(doIt)  = bnd4N;
    bnd5(doIt)  = bnd5N;
    ins(doIt)   = insN;
    tau(doIt)   = tau_;
    nit(doIt)   = nit(doIt) + 1;
    doIt(doIt)  = (abs(cost_-costN)>(dC*dC) ) | ...
        (abs(f_-fN)>dl*333) | ...
        (abs(delta_-deltaN)>dl*1000) | fails;
    n           = n+1;
    % ---------------------------------------------------------------------
end
lpar = lperp + delta; % M x 1
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
function [bnd2,bnd3,bnd4,bnd5,ins,lperp,delta,f] = ...
    project_on_boundaries(lperp,delta,f,ADC0,fmin,fmax)
bnd2        = (delta<0);          % M x 1, second frontier
delta(bnd2) = 0;
bnd3        = (lperp+delta>ADC0); % M x 1, third frontier
if(any(bnd3))
   delta(bnd3) = max(ADC0-lperp(bnd3),0);
end
bnd2(bnd3) = false;
bnd4       = (f<fmin);
f(bnd4)    = fmin(bnd4);
bnd2(bnd4) = false;
bnd3(bnd4) = false;
bnd5       = (f>fmax);
f(bnd5)    = fmax(bnd5);
bnd2(bnd5) = false;
bnd3(bnd5) = false;
bnd4(bnd5) = false;
ins        = ~(bnd2|bnd3|bnd4|bnd5);
