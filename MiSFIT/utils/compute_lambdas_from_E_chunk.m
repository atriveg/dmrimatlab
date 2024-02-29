function [lpar,lperp,nit] = compute_lambdas_from_E_chunk(E,b,mu,nu,mlp,nmax,dl,dC,ADC0)
% function [lpar,lperp,nit] = compute_lambdas_from_E_chunk(E,b,mu,nu,mlp,nmax,dl,dC,ADC0)
%
%  Computes the two components of the micro-structural model given the
%  attenuation signal. This is a helper function to atti2micro.
%
%    E: N x NB, the signal in the natural domain.
%    b: NB x 1, the set of b-values for each channel in the last dimension
%       of E
%    mu: 1 x 1, penalty weight to avoid lperp being extremely small
%    nu: 1 x 1, penalty weight to avoid lperp being extremely close to lpar
%    mlp: 1 x 1, minimum allowed value for lperp
%    nmax: 1 x 1, maximum number of Newton-Raphson iterations
%    dl: maximum allowed chage in lambda from one iteration to the next
%    dC: maximum allowed change in the (squared) norm of the vector
%       objective funcion
%    ADC0: 1x1, free water diffusivity at human body temperature
%    mlp: 1x1, minimum allowed value for lperp
%
%    lpar: N x 1, the diffusivity in the parallel
%       direction.
%    lperp: N x 1, the diffusivity in the perpendicular
%       direction.
%    nit: N x 1, the number of iterations needed to
%       obtain the solution
%

%%% =======================================================================
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
%%% =======================================================================
% Normalize b to avoid numerical issues (this scales diffusivities 1000
% times)
b  = b(:)/1000; % NB x 1
b2 = b.*b;      % NB x 1
b3 = b2.*b;     % NB x 1
b4 = b2.*b2;    % NB x 1
% Compute the Least Squares matrix for the first iteration (2nd order
% polynomial):
LS   = [ sum(b2), sum(b3); sum(b3), sum(b4) ];  % 2 x 2
LS   = LS\eye(2);                               % 2 x 2
%%% =======================================================================
% Normalize diffusivity-related parameters as well:
mlp  = 1000*mlp;
ADC0 = 1000*ADC0;
%%% =======================================================================
% Check sizes and compute the log-signal:
M  = size(E,1);
El = log(E); % M x NB
%%% =======================================================================
% Find a first iteration
ND    = 2;                   % 1 x 1
NI    = 3;                   % 1 x 1
d0v   = linspace(ADC0,0,ND); % 1 x ND
cost  = zeros(M,ND);         % M x ND
lpar  = zeros(M,ND);         % M x ND
lperp = zeros(M,ND);         % M x ND
for nd=1:ND
    % ---------------------------------------------------------------------
    % The log-sqrt-erf decay is modeled as a second order polynomial:
    %
    %  El(bi) = -bi*lperp - (bi*(d0/3) - bi^2*(2*d0^2/45) +  resid(bi*d0))
    %
    % So that we solve a quadratic regression problem in "b" given a
    % temptative solution d0:
    %
    %  [El(bi)-resid(bi*d0_temp)] = [ bi,   bi^2 ] * [ x1,   x2 ]^T,
    %    with:
    %              x1 = -lperp-d0/3;
    %              x2 = 2*d0^2/45.
    %
    % We can solve for a new estimate of lperp and d0 as:
    %
    %       d0 = 3*sqrt(5*x2/2);
    %    lperp = -x1 - d0/3 = -x1 - sqrt(5*x2/2);
    %     lpar = lperp + d0 = -x1 + sqrt(10*x2),
    %
    % update d0temp=d0 and iterate.
    % ---------------------------------------------------------------------
    % Initial temptative value of d0:
    d0_nd = ones(M,1)*d0v(nd); % M x 1
    % ---------------------------------------------------------------------
    for ni=1:NI
        % -----------------------------------------------------------------
        % Compute the target depending on the current estimate of d0
        tgt = El + sqrt_over_erfL_correction(d0_nd*(b')); % M x NB
        % -----------------------------------------------------------------
        % Solve the least squares problem using a precomputed matrix:
        if(is_broadcast_available)
            y1 = sum(tgt.*(b'),2);    % M x 1
            y2 = sum(tgt.*(b2'),2);   % M x 1
        else
            y1 = sum( bsxfun( @(x,y)(x.*y), tgt, b' ),  2 ); % M x 1 
            y2 = sum( bsxfun( @(x,y)(x.*y), tgt, b2' ), 2 ); % M x 1
        end
        x1 = LS(1,1)*y1 + LS(1,2)*y2; % M x 1
        x2 = LS(2,1)*y1 + LS(2,2)*y2; % M x 1
        % -----------------------------------------------------------------
        % Solve for the parameters of interest depending on the
        % characteristics of x1 and x2
        lpe_nd = zeros(M,1);           % M x 1
        % ------------------- Case study:
        pp1 = (x1<=0);
        pp2 = (x2>=0);
        % ------------------- Regular situation:
        pp         = (pp1 & pp2);
        d0_nd(pp)  = 3*sqrt((5/2)*x2(pp));
        lpe_nd(pp) = max( -x1(pp) - d0_nd(pp)/3, 0 );
        % ------------------- Weird value of x1
        pp         = ( (~pp1) & pp2 );
        d0_nd(pp)  = 3*sqrt((5/2)*x2(pp));
        lpe_nd(pp) = mlp;
        % ------------------- Weird value of x2
        pp         = ( pp1 & (~pp2) );
        d0_nd(pp)  = 0;
        lpe_nd(pp) = -x1(pp);
        % ------------------- Weird value of both x1 and x2
        pp         = ( (~pp1) & (~pp2) );
        d0_nd(pp)  = min(mlp+ADC0/10,ADC0-mlp);
        lpe_nd(pp) = mlp;
        % -----------------------------------------------------------------
        % Now, d0_nd has been properly updated and a new estimate of the
        % residual is available
        % -----------------------------------------------------------------
    end
    % ---------------------------------------------------------------------
    % By now the fixed-point-like iterations have finished, and we can
    % compute the final cost achieved for the present initial temptative
    % value of d0:
    C           = El + lpe_nd*(b') + sqrt_over_erfL(d0_nd*(b')); % M x NB
    cost(:,nd)  = sum(C.*C,2);                                   % M x ND
    lpar(:,nd)  = d0_nd + lpe_nd;                                % M x ND
    lperp(:,nd) = lpe_nd;                                        % M x ND
    % ---------------------------------------------------------------------
end
% We have tried all possible initial d0. Keep the solution with minimum
% error:
if(ND>1)
    [~,idx] = min( cost, [], 2 );
    lpar    = lpar(  sub2ind([M,ND],(1:M)',idx) );
    lperp   = lperp( sub2ind([M,ND],(1:M)',idx) );
    delta   = lpar - lperp;
end
%%% =======================================================================
% Make sure the solution is physically meaningful by just projecting the
% temptative solutions onto the feasible region:
[bnd1,bnd2,bnd3,ins,lperp,delta] = ...
    project_on_boundaries(lperp,delta,ADC0,mlp);
%%% =======================================================================
% Constrained Newton-Raphson iterations
doIt = true(M,1);
tau  = 0.001*ones(M,1);
nit  = zeros(M,1);
n    = 0;
while( any(doIt) && (n<nmax) )
    % ---------------------------------------------------------------------
    % Process only the necessary data:
    El_    = El(doIt,:);   % P x NB
    lperp_ = lperp(doIt);  % P x 1
    delta_ = delta(doIt);  % P x 1
    tau_   = tau(doIt);    % P x 1
    bnd1_  = bnd1(doIt);   % P x 1
    bnd2_  = bnd2(doIt);   % P x 1
    bnd3_  = bnd3(doIt);   % P x 1
    ins_   = ins(doIt);    % P x 1
    % ---------------------------------------------------------------------
    % In general (i.e. for points inside the feasible region), we aim at
    % minimizing a cost function given by:
    %   Q = (1/2) sum_i ( El(bi) + lp*bi + f(bi*d) )^2 + mu*P(lp,d)
    %     = (1/2) sum_i c_i(lp,d)^2 + mu*P(lp,d)
    % in the two variables lp (lperp_) and d (delta_=lpar_-lperp_), where 
    % c_i(lp,d) = El(bi) + lp*bi + f(bi*d)and function f stands for the 
    % log-sqrt-erf decay. The minimum is reached
    %  when the gradient equals 0, i.e.:
    %   Q_lp = sum_i c_i(lp,d)*bi          + mu*P_x(lp,d) = 0
    %   Q_d  = sum_i c_i(lp,d)*f'(bi*d)*bi + mu*P_y(lp,d) = 0
    % or, in matrix form:
    %   [ bi, f'(d*bi)*bi ]^T*[ c_i(lp,d) ]
    %                                  + mu*[ P_x(lp,d), P_y(lp,d) ]^T = 0
    % By taking first order Taylor series expansions around the current
    % estimate of lp and d for the terms c_i, P_x, and P_y, we get a LLS
    % problem with closed form solution to obtain a new iteration
    %   [ bi, f'(d0*bi)*bi ]^T*[ c_i(lp0,d0) ]
    %      + [ bi, f'(d0*bi)*bi ]^T * [ bi, f'(d0*bi)*bi ]*[ u1, u2 ]^T
    %      + mu*[ P_x(lp0,d0), P_y(lp0,d0) ]^T
    %      + mu*[ P_xx(lp0,d0), P_xy(lp0,d0) ]
    %           [ P_xx(lp0,d0), P_xy(lp0,d0) ]*[ u1, u2 ]^T = 0
    % for u1 = (lp-lp0) and u2 = (d-d0). Hence:
    %   ( Jc^t*Jc + mu*Hp )*[ u1, u2 ]^T = - ( Jc^T*[ ci ] + mu*Jp^T )
    % When mu=0 (no penalty applies) this completely equivalent to make
    % Newton-Raphson iterarions for the non-linear problem:
    %   [ c_i ] = [ 0 ]. 
    % To avoid numerical issues, any matrix inversion is performed over a
    % regularized version: a term of the form tau*I_N is added, so that
    % if the iterations fail tau is increased and the method becomes alike
    % to gradient-descent. If the iterations succeed, tau is decreased so
    % that the method becomes a like to pure Newton-Raphson iterations and
    % its convergence is faster.
    % ---------------------------------------------------------------------
    % Compute the current cost function:
    argF   = delta_*(b');               % P x NB
    func_  = sqrt_over_erfL(argF);      % P x NB, the log-sqrt-erf decay
    cost_  = lperp_*(b') + func_ + El_; % P x NB
    % ---------------------------------------------------------------------
    % Compute the current jacobian:
    d1_ = b';                       % 1 x NB
    d2_ = sqrt_over_erfL_d1(argF);  % P x NB
    if(is_broadcast_available)
        d2_ = d2_.*(b');                       % P x NB
    else
        d2_ = bsxfun( @(x,y)(x.*y), d2_, b' ); % P x NB
    end
    % ---------------------------------------------------------------------
    % Compute the current penalty term:
    [PV_,Px_,Py_,Pxx_,Pxy_,Pyy_] = compute_penalty_term(lperp_,delta_,mu,nu);
    % ---------------------------------------------------------------------
    % Pseudo-invert the jacobian:
    J11 = sum(d1_.*d1_,2) + Pxx_;    % P x 1
    J22 = sum(d2_.*d2_,2) + Pyy_;    % P x 1
    LM  = tau_.*(J11+J22)/2;         % P x 1
    J11 = J11 + LM;                  % P x 1
    J22 = J22 + LM;                  % P x 1
    J12 = d2_*b + Pxy_;              % P x 1
    % --- Invert:
    Jdet = J11.*J22 - J12.*J12; % P x 1
    Ji11 = J22./Jdet;           % P x 1
    Ji22 = J11./Jdet;           % P x 1
    Ji12 = -J12./Jdet;          % P x 1
    Ji21 = -J12./Jdet;          % P x 1
    % ---------------------------------------------------------------------
    % Compute the step:
    f1 = cost_*b + Px_;            % P x 1
    f2 = sum(d2_.*cost_,2) + Py_;  % P x 1
    u1 = -Ji11.*f1 - Ji12.*f2;     % P x 1
    u2 = -Ji21.*f1 - Ji22.*f2;     % P x 1
    % ---------------------------------------------------------------------
    % Manage boundary points. There are three boundaries for lperp<mlp,
    % delta<0, and lperp+delta>ADC0. At each of them:
    %   - If the step points inside the feasible region, we do just
    %     nothing, and flag the voxel as "inside" since the step is taking
    %     it out from the boundary.
    %   - If the step points outside the feasible region, we design a new
    %     step with a Newton-Raphson iteration constrained to the
    %     corresponding boundary.
    % ------------ First boundary: lperp = mlp
    if( any(bnd1_) )
        % ---------
        pi_ = ( bnd1_ & (u1>0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd1_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd1_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        % Optimize the problem in one variable by fixing lp = mlp, i.e.
        % find the minimum in d of:
        %   Q = (1/2) sum_i ( El(bi) + bi*mlp + f(bi*d) )^2 + mu*P(mlp,d),
        % then:
        %   Q_d = sum_i c_i(mlp,d)*f'(bi*d)*bi + mu*P_y(mlp,d) = 0
        % or, in matrix form:
        %   [ f'(d0*bi)*bi ]^T*[ c_i(mlp,d) ] + mu*P_y(mlp,d) = 0
        % with the series expansion at d = d0:
        %     [ f'(d0*bi)*bi ]^T*[ c_i(mlp,d0) ] + mu*P_y(mlp,d0)
        %   + [ f'(d0*bi)*bi ]^T*[ f'(d0*bi)*bi ]*u2
        %   + mu*P_yy(mlp,d0)*u2 = 0
        if(any(po_))
            M22     = sum(d2_(po_,:).*d2_(po_,:),2) + Pyy_(po_);  % Pb1o x 1
            y2      = sum(d2_(po_,:).*cost_(po_,:),2) + Py_(po_); % Pb1o x 1
            u2(po_) = -y2./M22./(1+tau(po_)/2);                   % Pb1o x 1
            u1(po_) = 0;                                          % Pb1o x 1
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb1 = po_;               % P x 1
        % ---------
    else
        forceb1 = false(size(bnd1_)); % P x 1
    end
    % ------------ Second boundary: delta = 0
    if( any(bnd2_) )
        % ---------
        pi_ = ( bnd2_ & (u2>0) );  % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd2_ & (~pi_)  ); % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd2_(pi_) = false;        % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;         % P x 1
        % ---------
        % Optimize the problem in one variable by fixing d = 0, i.e. find
        % the minimum in d of:
        %   Q = (1/2) sum_i ( El(bi) + bi*lp + f(0) )^2 + mu*P(lp,0),
        % then:
        %   Q_lp = sum_i c_i(lp,0)*bi + mu*P_x(lp,0) = 0
        % or, in matrix form:
        %   [ bi ]^T*[ c_i(lp,0) ] + mu*P_x(lp,0) = 0
        % with the series expansion at lp = lp0:
        %     [ bi ]^T*[ c_i(lp,0) ] + mu*P_x(lp,0)
        %   + [ bi ]^T*[ bi ]*u1 + mu*P_xx(lp,0)*u1 = 0
        if(any(po_))
            M11     = sum(b2) + Pxx_(po_);       % Pb2o x 1
            y1      = cost_(po_,:)*b + Px_(po_); % Pb2o x 1
            u1(po_) = -y1./M11./(1+tau(po_)/2);  % Pb2o x 1
            u2(po_) = 0;                         % Pb2o x 1
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
        pi_ = ( bnd3_ & (u1+u2<0) ); % P x 1, this takes the voxel into the feasible region
        po_ = ( bnd3_ & (~pi_)  );   % P x 1, this takes the voxel outside the feasible region
        % ---------
        bnd3_(pi_) = false;          % P x 1, this is a regular Newton-Raphson step
        ins_(pi_)  = true;           % P x 1
        % ---------
        % Optimize the problem in one variable by fixing lp + d = ADC0:
        %   d  = t;
        %   lp = ADC0 - t;
        % i.e. find the minimum in d of:
        %   Q = (1/2) sum_i ( El(bi) + bi*(ADC0-t) + f(bi*t) )^2
        %                            + mu*P(ADC0-t,t),
        % then:
        %   dQ/dt = sum_i c_i(ADC0-t,t)*( f'(bi*t)*bi - bi ) 
        %             - mu*P_x(ADC0-t,t) + mu*P_y(ADC0-t,t) = 0
        % or, in matrix form:
        %   [ f'(bi*t)*bi-bi ]^T*[ c_i(ADC0-t,t) ]
        %             - mu*P_x(ADC0-t,t) + mu*P_y(ADC0-t,t) = 0
        % with the series expansion at t = t0:
        %     [ f'(bi*t0)*bi-bi ]^T*[ c_i(ADC0-t0,t0) ]
        %   - mu*P_x(ADC0-t0,t0) + mu*P_y(ADC0-t0,t0)
        %   + [ f'(bi*t0)*bi-bi ]^T*[ f'(bi*t0)*bi-bi ](t-t0)
        %   - mu*( -P_xx(ADC0-t0,t0) + P_xy(ADC0-t0,t0) )(t-t0)
        %   + mu*( -P_yx(ADC0-t0,t0) + P_yy(ADC0-t0,t0) )(t-t0)
        if(any(po_))
            M12 = d2_(po_,:);                     % Pb3o x NB
            if(is_broadcast_available)
                M12 = M12 - b';                   % Pb3o x NB
            else
                M12 = bsxfun(@(x,y)(x-y),M12,b'); % Pb3o x NB
            end
            y12 = sum( cost_(po_,:).*M12, 2 ) - Px_(po_) + Py_(po_); % Pb3o x 1
            M12 = sum( M12.*M12, 2 ) + ...
                ( Pxx_(po_) + Pyy_(po_) - 2*Pxy_(po_) );             % Pb3o x 1
            u2(po_) = -y12./M12./(1+tau(po_)/2);                     % Pb3o x 1
            u1(po_) = -u2(po_);                                      % Pb3o x 1
        end
        % ---------
        % We will have to make sure that those values moved within the
        % boundary will remain flagged as such
        forceb3 = po_;                % P x 1
        % ---------
    else
        forceb3 = false(size(bnd3_)); % P x 1
    end
    % ---------------------------------------------------------------------
    % Compute the new parameters:
    lperpN = lperp_ + u1;                    % P x 1
    deltaN = delta_ + u2;                    % P x 1
    % ---------------------------------------------------------------------
    % Recheck if the solutions are in bounds or otherwise they are upon the
    % frontiers of the feasible region:
    [bnd1N,bnd2N,bnd3N,~,lperpN,deltaN] = ...
        project_on_boundaries(lperpN,deltaN,ADC0,mlp);
    % But we have to make sure that those values processed as boundary
    % points are still boundary points
    bnd1N(forceb1) = true;
    bnd2N(forceb2) = true;
    bnd1N(bnd2N)   = false;
    bnd3N(forceb3) = true;
    bnd1N(bnd3N)   = false;
    bnd2N(bnd3N)   = false;
    insN = ~(bnd1N|bnd2N|bnd3N);
    % ---------------------------------------------------------------------
    % Compute the new cost...
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
    lperpN(fails) = lperp_(fails);
    deltaN(fails) = delta_(fails);
    bnd1N(fails)  = bnd1_(fails);
    bnd2N(fails)  = bnd2_(fails);
    bnd3N(fails)  = bnd3_(fails);
    insN(fails)   = ins_(fails);
    % ---------------------------------------------------------------------
    % Prepare for the new step
    lperp(doIt) = lperpN;
    delta(doIt) = deltaN;
    bnd1(doIt)  = bnd1N;
    bnd2(doIt)  = bnd2N;
    bnd3(doIt)  = bnd3N;
    ins(doIt)   = insN;
    tau(doIt)   = tau_;
    nit(doIt)   = nit(doIt) + 1;
    doIt(doIt)  = (abs(cost_-costN)>(dC*dC) ) | ...
        (abs(delta_-deltaN)>dl*1000) | ...
        (abs(lperp_-lperpN)>dl*1000) | fails;
    n           = n+1;
    % ---------------------------------------------------------------------
end
lpar = lperp + delta; % M x 1
%%% =======================================================================
% Scale back:
lpar  = lpar/1000;
lperp = lperp/1000;


%%% =======================================================================
% Compute the error between the actual (logarithmic) decay due to the
% parallel component and its second-order series expansion at x=0
function resid = sqrt_over_erfL_correction(x)
resid = sqrt_over_erfL(x) - (x/3-x.*x/(45/2));

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
function [bnd1,bnd2,bnd3,ins,lperp,delta] = project_on_boundaries(lperp,delta,ADC0,mlp)
bnd1        = (lperp<mlp);        % M x 1, first frontier
lperp(bnd1) = mlp;
bnd2        = (delta<0);          % M x 1, second frontier
delta(bnd2) = 0;
bnd1(bnd2)  = false;
bnd3        = (lperp+delta>ADC0); % M x 1, third frontier
if(any(bnd3))
    % -------------------------------------------------
    scl = ADC0./(lperp(bnd3)+delta(bnd3)); % P x 1
    lpe = lperp(bnd3);                     % P x 1
    del = delta(bnd3);                     % P x 1
    % -------------------------------------------------
    bnd31      = (scl.*lpe>=mlp);        % P x 1, we can scale without going beyond f1
    lpe(bnd31) = lpe(bnd31).*scl(bnd31); % P x 1
    del(bnd31) = del(bnd31).*scl(bnd31); % P x 1
    % -------------------------------------------------
    bnd31      = ~bnd31;       % P x 1
    lpe(bnd31) = mlp;          % P x 1
    del(bnd31) = ADC0 - mlp;   % P x 1
    % -------------------------------------------------
    lperp(bnd3) = lpe;
    delta(bnd3) = del;
    % -------------------------------------------------
end
bnd1(bnd3) = false;
bnd2(bnd3) = false;
ins        = ~(bnd1|bnd2|bnd3);

%%% =======================================================================
function [PV,Px,Py,Pxx,Pxy,Pyy] = compute_penalty_term(lperp,delta,mu,nu)
% ---------------------------
pp        = (lperp<100*sqrt(eps));
lperp(pp) = 100*sqrt(eps);
pp        = (delta<100*sqrt(eps));
delta(pp) = 100*sqrt(eps);
% ---------------------------
PV  =  nu*(lperp./delta) + mu*(delta./lperp); % P x 1, the penalty itself
% ---------------------------
Px  = -mu*(delta./(lperp.*lperp)) + nu./delta;
Py  =  mu./lperp - nu*(lperp./(delta.*delta));
% ---------------------------
Pxx =  2*mu*(delta./(lperp.*lperp.*lperp));
Pyy =  2*nu*(lperp./(delta.*delta.*delta));
Pxy = -mu./(lperp.*lperp)-nu./(delta.*delta);
