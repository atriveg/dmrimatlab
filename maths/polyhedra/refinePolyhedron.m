function polyhedron = refinePolyhedron(polyhedron)
% function polyhedron = refinePolyhedron(polyhedron)
%
%   This function takes a convex polyhedron made up of triangular facets 
%   whose vertices lay within the unit sphere centered at zero and
%   defined by the structure:
%
%          % NVx3, the set of vertices in the unit sphere, each row
%          represents a vertex:
%            >> polyhedron.vertices;
%          % NFx3, the set of triangular facets, each row is a
%          3-element vector of indices in the range [1,NV] pointing to the
%          three vertices that define the facet:
%           >> polyhedron.facets
%          % NEx2, the set of edges between pairs of vertices, each row is
%          a two-element vector of indices in the range [1,NV] pointing to
%          the two vertices that define the edge:
%           >> polyhedron.edges
%
%   and outputs a nw polyhedron defined by a similar structure and obtained
%   by dividing each triangular facet of the original polyhedon into 4
%   different facets by inserting new vertices at the mid-point of the
%   edges of the former polyhedron and projecting them onto the unit
%   sphere.
%
%   For both polyhedrons, since the facets are all triangles and each edge
%   is shared by two facets, it should hold:
%
%          3*NF = 2*NE
%
%   Since both polyhedrons are convex, it should hold also:
%
%          NF+NV = NE+2
%
%   So that:
%
%          NV = (NF+4)/2
%          NE = 3*NF/2
%          NF_out = 4*NF_in
%
%   You can check the result by simply plotting the polyhedron:
%
%     >> polyhedron = createIcosahedron;
%     >> patch( 'Vertices', polyhedron.vertices, ...
%               'Faces', polyhedron.facets, 'FaceColor', [0,0.5,0], ...
%               'FaceAlpha', 0.8 );
%     >> polyhedron = refinePolyhedron(polyhedron);
%     >> patch( 'Vertices', polyhedron.vertices, ...
%               'Faces', polyhedron.facets, 'FaceColor', [0,0,0.5], ...
%               'FaceAlpha', 0.5 );
%     >> axis('equal');
%     >> light;
%     >> lighting('phong');


%--------------------------------------------------------
% Original number of elements:
NF = size( polyhedron.facets, 1 );
NV = size( polyhedron.vertices, 1 );
NE = size( polyhedron.edges, 1 );
assert(size(polyhedron.facets,2)==3,'Facets must be triangles');
assert(size(polyhedron.vertices,2)==3,'Vertices must be 3-D vectors');
assert(size(polyhedron.edges,2)==2,'Edges are defined by two vertices');
assert(NV==(NF+4)/2,'Weird sizes! Check NV=(NF+4)/2');
assert(NE==3*NF/2,'Weird sizes! Check NE=3*NF/2');
%--------------------------------------------------------
% The final number of vertices is NV+NE:
NV2      = NV+NE;
vertices = zeros( NV2, 3 );
vertices(1:NV,:) = polyhedron.vertices;
%--------------------------------------------------------
% Determine the final number of edges:
NE2      = 2*NE+3*NF;
edges    = zeros( NE2, 2 );
%--------------------------------------------------------
% Determine the final number of facets:
NF2      = 4*NF;
facets   = zeros( NF2, 3 );
%--------------------------------------------------------
% Create a map of middle-points between pairs of vertices:
MAP      = -1.*ones(NV,NV);
for n=1:NE
    % Find the middle point of each edge, and create a new vertex in that
    % position:
    %   The position of the new vertex in the list is:
    newPos = NV+n;
    %   The two old vertices describing this edge have positions:
    vtx1p  = polyhedron.edges(n,1);
    vtx2p  = polyhedron.edges(n,2);
    %   So the new entry in the map is:
    MAP(vtx1p,vtx2p) = newPos;
    MAP(vtx2p,vtx1p) = newPos;
    %   The value of the new vertex is:
    newVtx = ( polyhedron.vertices(vtx1p,:) + polyhedron.vertices(vtx2p,:) )./2;
    %   Project onto the unit sphere:
    newVtx = newVtx./norm(newVtx);
    %   Store the new vertex:
    vertices(newPos,:) = newVtx;
end
% AT THIS POINT, WE ALREADY HAVE THE NEW LIST OF VERTICES
%--------------------------------------------------------
% Iterate through facets to compute the new facets; compute also the
% corresponding new edges.
for n=1:NF
    % Get the pointers to the new facet:
    facet = polyhedron.facets(n,:);
    % This implementation works only for triangled polyhedra, so we keep
    % only the first three vertices of the facet regardless on this is a
    % closed polyline.
    p1    = facet(1);
    p2    = facet(2);
    p3    = facet(3);
    % FIRST of the new 4 facets:
    facet1 = [ p1, MAP(p1,p2), MAP(p1,p3) ];
    % SECOND of the new 4 facets:
    facet2 = [ MAP(p1,p2), p2, MAP(p2,p3) ];
    % THIRD of the new 4 facets:
    facet3 = [ MAP(p2,p3), p3, MAP(p1,p3) ];
    % FOURTH of the new 4 facets:
    facet4 = [ MAP(p1,p2), MAP(p2,p3), MAP(p1,p3) ];
    % Store all the facets in the list:
    facets(4*(n-1)+1,:) = facet1;
    facets(4*(n-1)+2,:) = facet2;
    facets(4*(n-1)+3,:) = facet3;
    facets(4*(n-1)+4,:) = facet4;
    % Store the new edges:
    edges(2*NE+3*(n-1)+1,:) = [MAP(p1,p2),MAP(p2,p3)];
    edges(2*NE+3*(n-1)+2,:) = [MAP(p2,p3),MAP(p1,p3)];
    edges(2*NE+3*(n-1)+3,:) = [MAP(p1,p2),MAP(p1,p3)];
end
%--------------------------------------------------------
% Iterate through edges to compute the splitted edges:
for n=1:NE
    p1 = polyhedron.edges(n,1);
    p2 = polyhedron.edges(n,2);
    edges(2*(n-1)+1,:) = [ p1, MAP(p1,p2) ];
    edges(2*(n-1)+2,:) = [ MAP(p1,p2), p2 ];
end
%--------------------------------------------------------
% Rebuild the output structure:
polyhedron.facets   = facets;
polyhedron.edges    = edges;
polyhedron.vertices = vertices;

%--------------------------------------------------------
% Make sure all facets are properly oriented:
polyhedron = reorientFacets(polyhedron);


% -------------------------------------------------------------------------
function polyhedron = reorientFacets(polyhedron)
facets = polyhedron.facets;
for f=1:size(facets,1)
    p1 = polyhedron.vertices(facets(f,1),:);
    p2 = polyhedron.vertices(facets(f,2),:);
    p3 = polyhedron.vertices(facets(f,3),:);
    pc = (p1+p2+p3);
    pc = pc./norm(pc);
    n  = cross(p2-p1,p3-p1);
    n  = n./norm(n);
    if( pc*n'<0 )
        facets(f,:) = [ facets(f,1), facets(f,3), facets(f,2)];
    end
end
polyhedron.facets = facets;
