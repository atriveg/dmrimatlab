function [grads,linit,lend] = icosamplesSphere(L,varargin)
% function [grads,linit,lend] = icosamplesSphere(L,'opt1',val1,'opt2',val2...)
%
%    Creates (roughly) uniform samplings of the unit sphere based on
%    succesive refinements of an icosahedron. The sampling is
%    hierarchically organized in L succesive levels (each level provides
%    more samples than the pevious one) that are roughly interleaved (the
%    samples in one level try to fill the gaps left by the samples in the
%    previous levels), so that appropriate multi-shell samplings can be
%    obtained by assigning a diffeent b-value to each refinement level. In
%    brief:
%
%      - We start with an icosahedon, i.e. 12 vertices and 20 facets. In
%        the first level, l=1, the gradients will be taken from the 12
%        vertices. However, since the icosahedron is symetric and two
%        opposite gradients represent the same diffusion measurement, only
%        6 gradients are returned.
%
%      - For level l=2, the gradient directions corresponding to the 20
%        facets are chosen. However, the ctual amount of diffusion
%        gradients is 10 since we have to remove opposites.
%
%      - If L>2, the icosahedron is refiend by splitting each facet in 4
%        triangles, inserting new vertexes at the mid points of the
%        original edges (see refinePolyhedron.m). A gran total of 42
%        vertices is obtained this way, among which 30 of them are new
%        points (the remaining 12 where already in the original
%        icosahedron). From these new 30 points, we keep 15 of them and
%        remove their opposites.
%
%      - If L>3, the center points of the newly introduced facets are used
%        as gradients: the original 20 facets are refined into 80 facets,
%        so that 80-20=60 new center points are obtained. Only 30 up of
%        these 60 points are actually used since the  opposites must be
%        removed.
%
%      - This process is recursively repeated if L>4 to obtain the
%        following hierarchy of interleaved gradient directions:
%
%      L:           1   2   3   4   5   6   7    ...   L
%      grads_L:     6   10  15  30  60  120 240  ...   2^(L-3)*15
%      total_grads: 6   16  31  61  121 241 481  ...   16+15*(2^(L-2)-1)
%
%    L: a 1x1 integer, the number of levels in the hierarchy
%    grads: a Mx3 double, with M=6 for L=1 or M=16+15*(2^(L-2)-1) for L>1.
%      The complete set of gradients in the hierarchy. All of them are
%      normalized to modulus 1
%    linit, lend: both Lx1 integers describing the hierarchy. For l=1..L,
%      the gradients at level l in the hierarchy are recovered as:
%           grads_l = grads(linit(l):lend(l),:);
%
%    Optional arguments as pair/value:
%
%       O1: wether (true) or not (false) using optimizeGradients.m to
%           further improve the distribution of gradients at each level by
%           increasing the homogenity of their mutual distances (default:
%           false)
%
%    In case O1 is set true:
%
%       p: 1x1 double, it defines the energy term in [1] (default: 1)
%       iters: 1x1 double, maximum number of iterations for the
%          optimization (default: 100)
%       lambda: 1x1 double, Levenberg-Marquardt regularization term in the 
%          Newton-Raphson step (default: 1.0)
%       delta: maximum allowed change in the parameters before the
%          iterations stop (default: 1.0e-5)
%       verbose: 1x1 logical, print additional information (default: false)

% Parse the optional input arguments:
opt.O1 = false;      optchk.O1 = [true,true];      % always 1x1 boolean
opt.p = 1;           optchk.p = [true,true];       % always 1x1 double
opt.iters = 100;     optchk.iters = [true,true];   % always 1x1 double
opt.lambda = 1.0;    optchk.lambda = [true,true];  % always 1x1 double
opt.delta = 1.0e-5;  optchk.delta = [true,true];   % always 1x1 double
opt.verbose = false; optchk.verbose = [true,true]; % always 1x1 boolean
opt = custom_parse_inputs(opt,optchk,varargin{:});

assert(numel(L)==1,'Input L must be scalar');
assert(L==floor(L),'Input L must be an integer');
assert(L>=1,'Input L must be greater than or equal to 1');

if(L==1)
    m = 6;
elseif(L==2)
    m = [6;10];
else
    m = [6;10;15*2.^(0:L-3)'];
end

M     = cumsum(m); % 6, 16, 31, 61, 121...
lend  = M;
linit = M-m+1;

grads = zeros(M(end),3);

% Start by creating an icoshedron and its successive refinements
nref       = floor((L-1)/2); % Each refinement provides two levels
polyhedron = cell(1,nref+1);
polyhedron{1} = createIcosahedron;
for n=1:nref
    polyhedron{n+1} = refinePolyhedron(polyhedron{n});
    % If necessary, optimize the gradient directions to obtain a more
    % homogeneous distribution:
    if(opt.O1)
        cverts = polyhedron{n+1}.vertices;
        [cverts2,idx1,idx2] = remove_duplicated_gradients(cverts);
        cverts2 = optimizeGradients( cverts2, 'p', opt.p, ...
            'iters', opt.iters, 'lambda', opt.lambda, ...
            'delta', opt.delta, 'verbose', opt.verbose );
        cverts(idx1,:) = cverts2;   % Optimized gradients
        cverts(idx2,:) = -cverts2;  % Their opposites
        % Now, update the vertices in the whole hierarchy. Forntunately, nw
        % vertices are added at the end, so this task is easy:
        for il=1:n+1
            iEnd = size( polyhedron{il}.vertices, 1 );
            polyhedron{il}.vertices = cverts(1:iEnd,:);
        end
    end
end

% Now, pick up gradients from either the vertices or the facets:
for l=1:L
    pInd = floor((l-1)/2)+1;
    if(rem(l,2)==1) % Take gradients from vertices
        % refinePolyhedron insert new vertices in the last positions, so it
        % is very easy to pick up the appropriate ones knowing how many
        % vertices we must take in this level:
        newGrads = polyhedron{pInd}.vertices(end-(2*m(l))+1:end,:); % Pick up the last 2*m(l) vertices
    else % Take gradients from facets
        % refinePolyhedron substitutes each old facet with 4 new facets. In
        % case l=1, we must choose all facets. In case l>1, we must discard
        % those inserted in the central parts since they ae duplicates:
        facets = polyhedron{pInd}.facets;
        if(l>2)
            idx = 4:4:size(facets,1);
            idx = setdiff(1:size(facets,1),idx);
            facets = facets(idx,:);
        end
        % Now we have the proper facets. Find their baricenters:
        vertices = polyhedron{pInd}.vertices;
        newGrads = vertices(facets(:,1),:) + ...
            vertices(facets(:,2),:) + vertices(facets(:,3),:);
        moduli   = sqrt(sum(newGrads.*newGrads,2));
        newGrads = newGrads./repmat(moduli,[1,3]);
    end
    % Either way, we have twice the gradients we need because there are
    % duplicates. Remove them and make sure all gradients lay within the
    % sub-space where 0<=phi<pi:
    [newGrads,~,~] = remove_duplicated_gradients(newGrads);
    % Place gradients in the north hemisphere:
    newGrads( newGrads(:,3)<0, : ) = -newGrads( newGrads(:,3)<0, : );
    % Now there are no duplicates. Put the new gradients in place:
    grads(linit(l):lend(l),:) = newGrads;
    % If needed, improve gradient positions coming from facets:
    if( (l>1) && (opt.O1) )
        % Find gradients corresponding to vertices of the polyhedra:
        idx = linit(1):lend(1);
        for il=3:2:L
            idx = [idx,(linit(il):lend(il))]; %#ok<AGROW> % Odd levels provide vertex-based gradients
        end
        % Move only facet-based gradients:
        grads(1:lend(l),:) = optimizeGradients( grads(1:lend(l),:), ...
            'p', opt.p, 'iters', opt.iters, 'lambda', opt.lambda, ...
            'delta', opt.delta, 'verbose', opt.verbose, ...
            'exclude', idx );
    end
end
