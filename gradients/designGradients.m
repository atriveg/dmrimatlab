function [G,sets] = designGradients(N,varargin)
% function [G,sets] = designGradients(N,'Option1',value1,'Option2',value2...)
%
%    Designs single-shell and multi-shell samplings evenly spaced over
%    (half) the surface of the unit sphere. N can be either a scalar 
%    integer (for single-shell samplings) or a vector of integers (for
%    multi-shell samplings), for which each entry represents the number of
%    gradients per shell. The algorithm proceeds in two steps:
%
%      1- Using spiralPhyllotaxis, it finds a roughly evenly spaced
%         distribution fulfilling the requirements on the number of 
%         gradients per shell.
%      2- Using either optimizeGradients (for single-shell samplings) or
%         optimizeGrradientSets (for multi-shell samplings), it refines the
%         positions of the gradients previously found in step 1.
%
%    The second output sets has size sum(N)x1, and values in the range
%    1,2,...,length(N). Its entries tell which shell the corresponging
%    gradient (row within G) belongs to.
%
%    OPTIONAL ARGUMENTS:
%
%       plot: whether (true) or not (false) plot the results of the
%         optimization (default: false)
%
%    These arguments are directly passed to either optimizeGradients or 
%    optimizeGradientSets, and their meaning can be checked in the help
%    comments therein:
%
%       p: 1x1 double, (default: 1)
%       iters: 1x1 double, (default: 100)
%       lambda: 1x1 double, (default: 1.0)
%       delta: 1x1 double, (default: 1.0e-5)
%       verbose: 1x1 logical, (default: false)
%
%    For multi-shell samplings only:
%
%       weights: 1xlength(N), (default: [1;0;0;...;1])
%       nocomb: 1x1 logical, (default: false)

% Parse the optional input arguments:
%%%
opt.plot = false;    optchk.plot = [true,true];    % always 1x1 boolean
opt.p = 1;           optchk.p = [true,true];       % always 1x1 double
opt.iters = 100;     optchk.iters = [true,true];   % always 1x1 double
opt.lambda = 1.0;    optchk.lambda = [true,true];  % always 1x1 double
opt.delta = 1.0e-5;  optchk.delta = [true,true];   % always 1x1 double
opt.verbose = false; optchk.verbose = [true,true]; % always 1x1 boolean
opt.weights      = zeros(1,length(N));
opt.weights(1)   = 1;
opt.weights(end) = 1;
optchk.weights = [true,true];
opt.nocomb = false;  optchk.nocomb = [true,true];  % always 1x1 boolean
%%%
opt = custom_parse_inputs(opt,optchk,varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[G,sets] = spiralPhyllotaxis(N);
if(length(N)==1) % Single-shell
    G = optimizeGradients( G, ...
        'p', opt.p, ...
        'iters', opt.iters, ...
        'lambda', opt.lambda, ...
        'delta', opt.delta, ...
        'verbose', opt.verbose );
else % Multi-shell
    G = optimizeGradientSets( G, sets, ...
        'p', opt.p, ...
        'iters', opt.iters, ...
        'lambda', opt.lambda, ...
        'delta', opt.delta, ...
        'verbose', opt.verbose, ...
        'weights', opt.weights, ...
        'nocomb', opt.nocomb );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(opt.plot)
    polyhedron = createIcosahedron;
    for n=1:3
        polyhedron = refinePolyhedron(polyhedron);
    end
    close(figure(2023));
    hf = figure(2023);
    % ------------------------------------
    subplot(1,2,1); hold('on'); rotate3d('on');
    patch( 'Vertices', polyhedron.vertices, ...
        'Faces', polyhedron.facets, 'FaceColor', [0.2,0.2,0.4], ...
        'FaceAlpha', 0.9, 'EdgeAlpha', 0 );
    light;
    lighting('phong');
    axis('equal');
    axis('off');
    title('All shells','FontSize',32);
    % ------------------------------------
    subplot(1,2,2); hold('on'); rotate3d('on');
    patch( 'Vertices', polyhedron.vertices, ...
        'Faces', polyhedron.facets, 'FaceColor', [0.2,0.2,0.4], ...
        'FaceAlpha', 0.9, 'EdgeAlpha', 0 );
    light;
    lighting('phong');
    axis('equal');
    axis('off');
    title('Current shell','FontSize',32);
    % ---------------------------------------------------------------------
    ng = 1;
    ni = 1;
    colors = [ 1,0,0; 0,1,0; 0,0,1; 1,1,0; 1,0,1; 0,1,1 ];
    for s=1:length(N)
        % ---------------------------------------------------------
        ne = ni + N(s) - 1;
        % ---------------------------------------------------------
        color = colors( rem(s-1,size(colors,1))+1, : );
        hp = zeros(1,N(s));
        for n=1:N(s)
            % -------------------------------------------
            subplot(1,2,1);
            plot3( ...
                [G(ng,1),-G(ng,1)], ...
                [G(ng,2),-G(ng,2)], ...
                [G(ng,3),-G(ng,3)], ...
                'LineStyle', 'none', ...
                'Marker', '.', ...
                'MarkerSize', 60, ...
                'Color', color );
            axis('equal');
            % -------------------------------------------
            subplot(1,2,2);
            hp(n) = plot3( ...
                [G(ng,1),-G(ng,1)], ...
                [G(ng,2),-G(ng,2)], ...
                [G(ng,3),-G(ng,3)], ...
                'LineStyle', 'none', ...
                'Marker', '.', ...
                'MarkerSize', 60, ...
                'Color', color );
            axis('equal');
            % -------------------------------------------
            pause(0.1);
            % -------------------------------------------
            ng = ng+1;
        end
        % ---------------------------------------------------------
        subplot(1,2,2);
        title('Current shell [HIT ENTER FOR THE NEXT ONE]','FontSize',32);
        pause;
        delete(hp);
        title('Current shell','FontSize',32);
        % ---------------------------------------------------------
        ni = ne + 1;
    end
    close(hf);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
