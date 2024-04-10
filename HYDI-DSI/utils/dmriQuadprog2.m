function [x,res,iters,y] = dmriQuadprog2(Q,f,A,a,B,b,lb,ub,miters,otol,stol)
% [x,res,iters,y] = dmriQuadprog2( Q, f, A, a, B, b, lb, ub, 
%              [miters, [otol, [stol, [ctol, [initl]]]]] )
%
%    Solves a very specific quadratic program typically appearing in
%    diffusion MRI:
%
%         min_x   0.5 · x^T·Q^T·x + f^T·x
%         s.t.:   A*x <= a
%                 B*x  = b
%                 x   >= lb
%                 x   <= ub
%
%    which typically fits some positive function to noisy measurements by
%    means of least squares. The equality constraint often answers to the
%    physical constraint of the searched function to have mass 1.
%
%    This function is designed because Matlab's quadprog seems to have a
%    poor performance for large scale problems.
%
%  INPUTS:
%
%    Q:  NxN matrix of doubles, symmetric and positive definite
%    f:  a Nx1 vector of doubles
%    A:  a M x N matrix
%    a:  a N x 1 vector
%    B:  a P x N matrix
%    b:  a P x 1 vector
%    lb: a N x 1 vector
%    ub: a N x 1 vector
%
%  OPTIONAL POSITIONAL INPUT ARGUMENTS:
%
%    miters: 1x1, maximum number of iterations (default: 100)
%    otol: 1x1, tolerance in the optimality condition: if the relative
%        change in the cost function w.r.t. that in the previous iteration
%        becomes smaller than otol, exit (default: 1.0e-6)
%    stol: 1x1, tolerance in the optimizations step: if the optimization
%        step becomes smaller than stol times the Newton-Raphson step, exit
%        (default: 1.0e-6)
%
%       NOTE: the sizes and types of these arguments,as well as their
%             consistency, are not checked for the sake of performance.
%
%  OUTPUTS:
%
%    x:     Nx1, the estimated positive function minimizing the cost function
%    res:   1x1, the residual 0.5 · x^T·A·x + b^T·x
%    iters: 1x1, number of iterations spent
%    y:     Cx1, Lagrange multipliers of the problem
if(nargin<9),   miters = 10000;   end
if(nargin<10),  otol   = 1.0e-16; end
if(nargin<11),  stol   = 1.0e-12; end

if(miters<1)
    error('at least 1 iteration required');
end

C = [B;A;-eye(size(lb,1));eye(size(ub,1))];
u = [b;a;-lb;ub];

% Normalize the rows of C and u (i.e. normalize the constraints) to
% attain a uniform RMS value.
rmsv = sqrt(sum(C.*C,2));
C    = C./rmsv;
u    = u./rmsv;

pe = false(size(u));
pe(1:size(b,1)) = true;
pn = ~pe;

Qi = Q\eye(size(Q,1));
g  = C*Qi*f + u;
H  = C*Qi*(C');

y0 = zeros(size(u));

y    = y0;
bnd  = ( (y0<0.0) & pn );
mu   = norm(g)/trace(Qi)/1000;
fct  = 1.0;
res0 = 0.5*(y0'*H*y0) + g'*y0;

streak = 3;
ngoods = 0;

for n=1:miters
    dy   = -mu*fct*(H*y0+g);
    bndn = bnd & (dy<0);
    dy(bndn) = 0;
    y    = y0 + dy;
    bndn = ( (y<0.0) & pn );
    y(bndn) = 0;
    res = 0.5*(y'*H*y)+g'*y;
    if(res<res0)
        % Success! Update x and increase fct
        ngoods = min(ngoods+1,streak);
        y0  = y;
        bnd = bndn;
        if(ngoods>streak-1)
            fct  = 2*fct;
        end
        % Also, if the change in the residual is small enough, we can exit:
        if( abs((res0-res)/res0) < otol )
            res0 = res;
            break;
        else
            res0 = res;
        end
    else
        % Wrong step. Try decreasing mu
        ngoods = 0;
        fct = fct/2;
    end
    % ---------------------------------------------------------------------
    % If the optimization step became too small without improving the
    % objective function, we should exit
    if(fct<stol)
        break;
    end
end
iters = n;
res   = res0;
x     = -Qi*(f+(C')*y);
end