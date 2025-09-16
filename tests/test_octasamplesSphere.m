% test_octasamplesSphere.m
L = 7;
[grads,linit,lend] = octasamplesSphere(L,'O1',true,'verbose',true);

pol = createOctahedron;
for n=1:4
    pol = refinePolyhedron(pol);
end

hf = figure(1001);
close(hf);
hf = figure(1001);
hold('on');
grid('on');
rotate3d('on');

patch( 'Vertices', pol.vertices, ...
    'Faces', pol.facets, 'FaceColor', [0.2,0.2,0.2], ...
    'FaceAlpha', 0.5, 'EdgeAlpha', 0 );
light;
lighting('gouraud');

colors = [ ...
    0.5, 0.0, 0.0;
    0.0, 0.5, 0.0;
    0.0, 0.0, 0.5;
    0.3, 0.3, 0.0;
    0.3, 0.0, 0.3;
    0.0, 0.3, 0.3;
    0.4, 0.2, 0.0;
    0.2, 0.0, 0.4;
    0.0, 0.4, 0.2];
for l=1:L
    plot3( grads(linit(l):lend(l),1), ...
        grads(linit(l):lend(l),2), ...
        grads(linit(l):lend(l),3), ...
        'LineStyle', 'none', 'Marker', '.', ...
        'MarkerSize', 40, 'Color', colors(l,:) );
    plot3( -grads(linit(l):lend(l),1), ...
        -grads(linit(l):lend(l),2), ...
        -grads(linit(l):lend(l),3), ...
        'LineStyle', 'none', 'Marker', '.', ...
        'MarkerSize', 40, 'Color', colors(l,:) );
    axis('equal');
    pause;
    title('PRESS ENTER to continue');
end
axis('equal');

