function tetrahedron = createTetrahedron
% tetrahedron = createTetrahedron
%
%   Creates a structure defining an tetrahedron with vertices over the unit
%   sphered and centered at zero:
%
%     tetrahedron.vertices: 4x3, each row is a vertex.
%     tetrahedron.facets: 4x3, each row is a vector with 3 indices in the
%          range [1,4] pointing to the vertices defining the facet.
%     tetrahedron.edges: 6x2, each row is a vector with 2 indices in the
%          range [1,4] pointing to the vertices defining the edge.
%
%   You can check the result by simply plotting the icosahedron:
%
%     >> polyhedron = createTetrahedron
%     >> patch( 'Vertices', polyhedron.vertices, ...
%               'Faces', polyhedron.facets, 'FaceColor', [0,0.5,0], ...
%               'FaceAlpha', 0.8 );
%     >> axis('equal');
%     >> light;
%     >> lighting('phong');
vertices = [    0,0,0;
                1,0,0;
                1/2,sqrt(3)/2,0
                1/2,sqrt(3)/6,sqrt(6)/3];
C        = mean( vertices );
for k=1:4
    vertices(k,:) = vertices(k,:) - C;
    vertices(k,:) = vertices(k,:)./norm(vertices(k,:));
end

facets   = [ [1,2,3]; [1,2,4]; [2,3,4]; [1,3,4] ];
edges    = [1,2;2,3;3,1;1,4;2,4;3,4];

tetrahedron.vertices = vertices;
tetrahedron.facets   = facets;
tetrahedron.edges    = edges;