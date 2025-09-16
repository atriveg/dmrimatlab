% test_optimizeGradients
N = 4;
type = 'i';

for n=1:N
    if(n==1)
        switch(type)
            case 'i'
                pol = createIcosahedron;
            case 'o'
                pol = createOctahedron;
        end
    else
        pol = refinePolyhedron(pol);
    end
end

hf = figure(1001);
close(hf);
hf = figure(1001);

subplot(1,2,1);
patch( 'Vertices', pol.vertices, ...
    'Faces', pol.facets, 'FaceColor', [0,0.5,0], ...
    'FaceAlpha', 0.8 );
axis('equal');
light;
rotate3d('on');
lighting('gouraud');

%%% -------
for n=1:N
    if(n==1)
        switch(type)
            case 'i'
                pol = createIcosahedron;
            case 'o'
                pol = createOctahedron;
        end
    else
        pol = refinePolyhedron(pol);
    end
    vts  = pol.vertices;
    G    = size(vts,1)/2;
    idx  = 1:2*G;
    idx2 = 1:G;
    dp   = vts(:,1)*vts(:,1)'+vts(:,2)*vts(:,2)'+vts(:,3)*vts(:,3)';
    dp   = abs(dp);
    dp(logical(eye(2*G))) = 0;
    for g=1:G
        [~,pos] = max(dp(idx(g),:));
        idx2(g) = pos;
        idx = setdiff(idx,pos);
    end
    vts = vts(idx,:);
    vts = optimizeGradients(vts,'verbose',true);
    vts2 = pol.vertices;
    vts2(idx,:)  = vts;
    vts2(idx2,:) = -vts;
    pol.vertices = vts2;
end
%%% -------
subplot(1,2,2);
patch( 'Vertices', pol.vertices, ...
    'Faces', pol.facets, 'FaceColor', [0,0.5,0], ...
    'FaceAlpha', 0.8 );
axis('equal');
light;
rotate3d('on');
lighting('gouraud');
