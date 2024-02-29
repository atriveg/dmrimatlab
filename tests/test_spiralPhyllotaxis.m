% test_spiralPhyllotaxis.m
Ns     = [8,16,32,32,64,64];
[G,ps] = spiralPhyllotaxis(Ns);
G      = optimizeGradientSets( G, ps, ...
    'iters', 1000, ...
    'delta', 1.0e-9, ... 
    'lambda', 10, ...
    'p', 1, ...
    'weights', [1,1,0,0,0,1], ...
    'verbose', true );
% -------------------------------------------------------------------------
polyhedron = createIcosahedron;
for n=1:3
    polyhedron = refinePolyhedron(polyhedron);
end
% -------------------------------------------------------------------------
close(figure(1003));
hf = figure(1003);
% ------------------------------------
subplot(1,2,1);
hold('on');
rotate3d('on');
patch( 'Vertices', polyhedron.vertices, ...
    'Faces', polyhedron.facets, 'FaceColor', [0.2,0.2,0.4], ...
    'FaceAlpha', 0.9, 'EdgeAlpha', 0 );
light;
lighting('phong');
axis('equal');
axis('off');
title('All shells','FontSize',32);
% ------------------------------------
subplot(1,2,2);
hold('on');
rotate3d('on');
patch( 'Vertices', polyhedron.vertices, ...
    'Faces', polyhedron.facets, 'FaceColor', [0.2,0.2,0.4], ...
    'FaceAlpha', 0.9, 'EdgeAlpha', 0 );
light;
lighting('phong');
axis('equal');
axis('off');
title('Current shell','FontSize',32);
% -------------------------------------------------------------------------
ng = 1;
ni = 1;
colors = [ 1,0,0; 0,1,0; 0,0,1; 1,1,0; 1,0,1; 0,1,1 ];
for s=1:length(Ns)
    % ---------------------------------------------------------
    ne = ni + Ns(s) - 1;
    % ---------------------------------------------------------
    color = colors(s,:);
    hp = zeros(1,Ns(s));
    for n=1:Ns(s)
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
% -------------------------------------------------------------------------
