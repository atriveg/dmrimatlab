function octahedron = createOctahedron
% function octahedron = createOctahedron
%
%   Creates a structure defining an octahedron with vertices over the unit
%   sphered and centered at zero:
%
%     octahedron.vertices: 6x3, each row is a vertex.
%     octahedron.facets: 8x3, each row is a vector with 3 indices in the
%          range [1,6] pointing to the vertices defining the facet.
%     octahedron.edges: 12x2, each row is a vector with 2 indices in the
%          range [1,6] pointing to the vertices defining the edge.
%
%   You can check the result by simply plotting the icosahedron:
%
%     >> polyhedron = createOctahedron
%     >> patch( 'Vertices', polyhedron.vertices, ...
%               'Faces', polyhedron.facets, 'FaceColor', [0,0.5,0], ...
%               'FaceAlpha', 0.8 );
%     >> axis('equal');
%     >> light;
%     >> lighting('phong');

vertices = [
    -sqrt(2)/2,-sqrt(2)/2,0;
    sqrt(2)/2,-sqrt(2)/2,0;
    sqrt(2)/2,sqrt(2)/2,0;
    -sqrt(2)/2,sqrt(2)/2,0;
    0,0,1;
    0,0,-1];

facets = [ [1,2,5]; [2,3,5]; [3,4,5]; [1,4,5]; [1,2,6]; [2,3,6]; [3,4,6]; [1,4,6] ];
edges  = [1,2;1,4;1,5;1,6;3,2;3,4;3,5;3,6;2,5;4,5;4,6;2,6];

octahedron.vertices = vertices;
octahedron.facets   = facets;
octahedron.edges    = edges;