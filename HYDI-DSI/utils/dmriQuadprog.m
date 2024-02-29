function [x,res,iters] = dmriQuadprog(A,b,f,x0,miters,otol,stol,ctol)
% [x,res,iters] = dmriQuadprog( H, b, f, x0,
%              [miters, [otol, [stol, [ctol]]]] )
%
%    Solves a very specific quadratic program typically appearing in
%    diffusion MRI:
%
%         min_x   0.5 · x^T·A·x + b^T·x
%         s.t.:   x>=0
%                 f^T·x = 1
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
%    A:  NxN matrix of doubles, symmetric and positive definite
%    b:  Nx1 vector of doubles
%    f:  Nx1 vector of doubles such that f(i)>0 for all i=1..N
%    x0: Nx1 vector of doubles, the initial estimate of the solution. NOTE:
%        it is assumed x0 fulfills all constraints. It WILL NOT be checked
%        internally.
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
%    ctol: 1x1, tolerance in the constraints (maximum distance to the
%        boundaries before the solution is considered to be over them) 
%        (default: 1.0e-8)
%
%       NOTE: the sizes and types of these arguments,as well as their
%             consistency, are not checked for the sake of performance.
%
%  OUTPUTS:
%
%    x:   Nx1, the estimated positive function minimizing the cost function
%    res: 1x1, the residual 0.5 · x^T·A·x + b^T·x
if(nargin<5), miters = 100;    end
if(nargin<6), otol   = 1.0e-6; end
if(nargin<7), stol   = 1.0e-6; end
if(nargin<8), ctol   = 1.0e-8; end

if(miters<1)
    x = x0;
    res = nan;
    iters = 0;
    return;
end

% Get rid of the equality constraint by redefining the variables of the
% problem:
H00 = A(1,1);
H10 = A(2:end,1);
H11 = A(2:end,2:end);
b0  = b(1,1);
b1  = b(2:end,1);
f0  = f(1,1);
f1  = f(2:end,1);
% ---
H   = H11 + (H00/(f0*f0))*(f1*f1') - (f1*H10'+H10*f1')/f0;
g   = b1 + H10/f0 - (H00/(f0*f0)+b0/f0)*f1;
dir = (H\g); % This is probably the slowest part!
% Find the "active boundaries", i.e. those where the components of x are
% null or where f1'*x=1:
x0   = x0(2:end);
bnd1 = (x0<ctol); % (N-1)x1
bnd2 = (f1'*x0>1-ctol);
% Start iterations (gradient projection)
mu      = 1/10;                 % This controls the step size
res0    = 0.5*(x0'*H*x0)+g'*x0; % This is the original residual
for n=1:miters
    % ---------------------------------------------------------------------
    % Move towards the Newton-Raphson calculated direction:
    dx = -mu*(x0+dir);
    % Check boundaries: if we are at a boundary point and the movement
    % points outside, we have to project:
    bnd1n     = bnd1 & (dx<0);
    dx(bnd1n) = 0;
    if(bnd2)
        if(f1'*dx>0)
            dx(~bnd1n) = dx(~bnd1n) - ...
                f1(~bnd1n)*(f1(~bnd1n)'*dx((~bnd1n)))/...
                (f1(~bnd1n)'*f1(~bnd1n));
        end
    end
    % Update:
    x = x0 + dx;
    % ---------------------------------------------------------------------
    % Make sure the solution is still in bounds, and update the active
    % boundaries:
    bnd1n    = (x<ctol); % (N-1)x1
    x(bnd1n) = 0;
    bnd2n    = (f1'*x>1-ctol);
    if(bnd2n)
        x = x./(f1'*x);
    end
    % ---------------------------------------------------------------------
    % Compute the new residual:
    res = 0.5*(x'*H*x)+g'*x;
    % Check if the residual is now smaller, so that the iterations
    % succeeded:
    if(res<res0)
        % Success! Update x and increase mu
        x0   = x;
        bnd1 = bnd1n;
        bnd2 = bnd2n;
        mu   = 2*mu;
        % Also, if the change in the residual is small enough, we can exit:
        if( abs((res0-res)/res0) < otol )
            res0 = res;
            break;
        else
            res0 = res;
        end
    else
        % Wrong step. Try decreasing mu
        mu = mu/2;
    end
    % ---------------------------------------------------------------------
    % If the optimization step became too small without improving the
    % objective function, we should exit
    if(mu<stol)
        break;
    end
    % ---------------------------------------------------------------------
end
iters = n;
% Finally, undo variable stretching:
res = res0;
x   = [ (1-f1'*x0)/f0; x0 ];
