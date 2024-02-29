function [G2_def,ptr_def] = decimateGradients(G,N2,varargin)
% function [G2,n2] = decimateGradients(G,N2,'Option1', 
%                                               value1,'Option2',value2...)
%
%    Decimates a set G of gradients (whose size is Nx3, each row being a
%    unit norm gradient) to find a set G2 of N2<N gradients taken from G
%    that are roughly evenly spaced:
%
%        G2 = G(n2,:);
%
%    Here, "roughly evenly spaced" means the output set of gradients, G2, 
%    try to minimize a figure of merit inspired by electrostatics, i.e., 
%    minimize the electrostatic potentials:
%
%    Q(n2) = sum_{i=1}^N2 sum_{j=1,j~=i}^N2 1/acos(|g_{n2(i)}^t g_{n2(j)}|)
%
%    where n2 is an 1xN2 vector with the selected indices among the N
%    gradients in the original collection G.
%
%    The optimization is merely heuristic, so no warranty of global
%    optimality is provided:
%
%      0- Set Q0 = Inf.
%      1- Repeat 'trials' times:
%        2- Randomly initialize n2 with randperm(N,N2).
%        3- Repeat at most 'maxiters' times:
%          4- For each position i=1..N2 in n2, subsequentially:
%             + Remove position i, m = setdiff(n2,n2(i)).
%             + Among all gradients j in the pool not in n2,
%               j = setdiff(1:N,n2), find the one such that [m,j] has
%               minimum energy, j0 = arg min_j{Q([m,j])}
%             + If Q([m,j0]) becomes smaller than Q(n2), update n2(i)=j0
%          5- If no indices in n2 have actually changed in step 4, stop 
%             repeating step 3 and go to 6
%        6- If the final cost Q obtained in step 3 is smaller than the 
%           current Q0, udpate Q0=Q and keep the current n2.
%
%    Options:
%
%        seed: seed to initalize the random number generator. Pass a
%           positive integer to obtain a predictable behavior each time the
%           function is invoked (default: [], meaning no re-seeding is
%           performed).
%        trials: the number of random initializations of the pool of
%           gradients chosen. The largest this number, the more likely a
%           globally optimal solution is obtained (default: 100).
%        maxiters: for each trial, the maximum number of optimization
%           iterations (i.e. the maximum number of loops along n2)
%           (default: 1000).
%        verbose: wether (true) or not (false) output additional
%           information about the optimization procedure (default: false)
%        plot: wether (true) or not (false) make a 3-D plot of the
%           resulting gradients (regd together with the original gradients
%           (green) and their respective specular reflections over the unit
%           sphere to check the results (default: false)
assert(ismatrix(G),'G must be a matrix');
assert(size(G,2)==3,'G must be Gx3');
N = size(G,1);
assert(N2<N,'N2 must be smaller than size(G,1)');
% Parse the optional input arguments:
opt.seed = [];       optchk.seed = [false,false];
opt.trials = 100;    optchk.trials = [true,true];   % always 1x1 double
opt.maxiters = 1000; optchk.maxiters = [true,true]; % always 1x1 double
opt.verbose = false; optchk.verbose = [true,true];  % always 1x1 boolean
opt.plot = false;    optchk.plot = [true,true];     % always 1x1 boolean
opt = custom_parse_inputs(opt,optchk,varargin{:});

% Normalize the gradient directions, just in case:
mods = sqrt(sum(G.*G,2));
G    = G./[mods,mods,mods];
% Compute all distances between each gradient and the remaining ones:
distances = compute_all_distances(G);
if(opt.verbose)
    Q = compute_Q_subset(distances,1:N);
    fprintf(1,'FOM in the orginal subset (the smaller the better): %1.5f\n',Q);
end
% Randomly intitalize the subset
if(~isempty(opt.seed))
    rng(opt.seed);
end
FOM = Inf;
for trial=1:opt.trials
    ptr = sort(randperm(N,N2));
    Q0  = compute_Q_subset(distances,ptr);
    Q   = Q0;
    if(opt.verbose)
        fprintf(1,'[%d] FOM in the initial random subset (the smaller the better): %1.5f\n',trial,Q);
    end
    for it=1:opt.maxiters
        count = 0;
        order = randperm(N2);
        for k=1:N2
            n    = order(k);
            Q0n  = Q;
            ptr2 = ptr([1:n-1,n+1:N2]);
            [Q,optimal] = find_optimal_entry(distances,ptr2);
            if(Q0n-Q>1.0e-3*Q0n)
                count = count+1;
                ptr(n) = optimal;
            end
        end
        if(opt.verbose)
            fprintf('[%d] FOM evolved from Q0=%1.5f to Q=%1.5f\n',trial,Q0,Q);
            fprintf('[%d] %d positions changed\n',trial,count);
        end
        if(count<0.99)
            break;
        else
            Q0 = Q;
        end
    end
    if(Q<FOM)
        FOM     = Q;
        G2_def  = G(ptr,:);
        ptr_def = ptr;
    end
end
if(opt.verbose)
    Q = compute_Q_subset(distances,ptr_def);
    fprintf(1,'FOM in the final subset (the smaller the better): %1.5f\n',Q);
end
if(opt.plot)
    hf = figure(3013);
    close(hf);
    figure(3013);
    subplot(1,2,1);
    plot3( [G(:,1);-G(:,1)], [G(:,2);-G(:,2)], [G(:,3);-G(:,3)], ...
        'Color', [0.0,0.5,0.0], 'LineStyle', 'none', 'Marker', '.', ...
        'MarkerSize', 16 );
    axis('equal');
    grid('on');
    rotate3d('on');
    subplot(1,2,2);
    plot3( [G2_def(:,1);-G2_def(:,1)], [G2_def(:,2);-G2_def(:,2)], ...
        [G2_def(:,3);-G2_def(:,3)], ...
        'Color', [0.5,0.0,0.0], 'LineStyle', 'none', 'Marker', '.', ...
        'MarkerSize', 16 );
    axis('equal');
    grid('on');
    rotate3d('on');
end

function distances = compute_all_distances(G)
[G1x,G1y]  = meshgrid(G(:,1),G(:,1));
[G2x,G2y]  = meshgrid(G(:,2),G(:,2));
[G3x,G3y]  = meshgrid(G(:,3),G(:,3));
css        = abs(G1x.*G1y+G2x.*G2y+G3x.*G3y);
css(css>1) = 1;
distances  = acos(css);
distances  = 1./distances;
distances(logical(eye(size(distances,1)))) = 0;

function Q = compute_Q_subset(distances,ptr)
distances = distances(ptr,ptr);                 % N2xN2
Q         = sum(distances(:));

function [Q,j0] = find_optimal_entry(distances,ptr)
N  = size(distances,1);          % Total number of gradients
ptr2       = setdiff(1:N,ptr);   % N3x1
distances2 = distances(ptr,ptr); % (N2-1)x(N2-1)
%%%
distances3 = distances(ptr,ptr2); % (N2-1)xN3
Qd         = 2*sum(distances3,1); % 1xN3
%%%
[Q,j0] = min(Qd);
j0     = ptr2(j0(1));
Q      = Q(1) + sum(distances2(:));
