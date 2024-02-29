function icosahedron = createIcosahedron
% function icosahedron = createIcosahedron
%
%   Creates a structure defining an icosahedron with vertices over the unit
%   sphered and centered at zero:
%
%     icosahedron.vertices: 12x3, each row is a vertex.
%     icosahedron.facets: 20x3, each row is a vector with 3 indices in the
%          range [1,12] pointing to the vertices defining the facet.
%     icosahedron.edges: 30x2, each row is a vector with 2 indices in the
%          range [1,12] pointing to the vertices defining the edge.
%
%   You can check the result by simply plotting the icosahedron:
%
%     >> polyhedron = createIcosahedron;
%     >> patch( 'Vertices', polyhedron.vertices, ...
%               'Faces', polyhedron.facets, 'FaceColor', [0,0.5,0], ...
%               'FaceAlpha', 0.8 );
%     >> axis('equal');
%     >> light;
%     >> lighting('phong');

phi = (1+sqrt(5))/2;
vertices = [
     [0,1, phi];
     [0,1,-phi];
     [1, phi,0];
     [1,-phi,0];
     [ phi,0,1];
     [-phi,0,1];
    -[0,1, phi];
    -[0,1,-phi];
    -[1, phi,0];
    -[1,-phi,0];
    -[ phi,0,1];
    -[-phi,0,1];
];

C = mean(vertices);
for n=1:size(vertices,1)
    vertices(n,:) = vertices(n,:)-C;
    vertices(n,:) = vertices(n,:)./norm(vertices(n,:));
end

edges = [
    1,3;
    1,5;
    1,8;
    1,6;
    1,10;
    7,2;
    7,4;
    7,9;
    7,11;
    7,12;
    5,8;
    8,6;
    6,10;
    10,3;
    3,5;
    4,12;
    12,2;
    2,11;
    11,9;
    9,4;
    3,2;
    2,10;
    10,11;
    11,6;
    6,9;
    9,8;
    8,4;
    4,5;
    5,12;
    12,3];

facets = [ ...
    [1,3,5];[1,3,10];[1,10,6];[1,6,8];[1,8,5];...
        [7,4,9];[7,9,11];[7,11,2];[7,2,12];[7,12,4];...
        [4,9,8];[4,8,5];[4,5,12];[12,5,3];[12,2,3];...
        [2,3,10];[2,10,11];[6,10,11];[11,6,9];[6,9,8] ...
        ];
    

icosahedron.vertices = vertices;
icosahedron.facets   = facets;
icosahedron.edges    = edges;