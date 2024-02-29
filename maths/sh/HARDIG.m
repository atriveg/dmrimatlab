function  [G,theta,phi] = HARDIG( G )
% function  [G,theta,phi] = HARDIG( G )
%
% Compute spherical coordinates from Cartesian components. G is assumed to
% be a Nx3 matrix, each row representing a point on the Carteisan space.
%
% All the points are put in the semisphere 0<=phi<pi, reflecting them with
% respect to the origin if needed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize all vectors to achieve norm 1:
for k=1:size(G,1)
    G(k,:) = G(k,:)./norm( G(k,:) );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the spherical coordinates:
tol  = eps;
% Initialize the output:
theta = zeros( size(G,1), 1 );
phi   = zeros( size(G,1), 1 );
for k=1:size(G,1) % For each point
    % Get the vector:
    g        = G(k,:);
    % Compute theta as the inverse cosine:
    theta(k) = acos( g(3) ); % Always in the range [0,pi]
    % Use atan2 to compute the phi (sin(theta(k)) is always >=0)
    if( (abs(g(1))<tol) && (abs(g(2))<tol) )
        % The atan2 is not defined in this case, so the best we can do
        % is placing a 0
        phi(k) = 0;
    else
        % atan2() outputs a value in the range [-pi,pi).
        phi(k) = atan2(g(2),g(1));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make sure all gradient directions are always in the y>=0 hemisphere:
for k=1:size(G,1)
    if( (pi-theta(k)) < tol ) % -z axis
        theta(k) = 0;
        phi(k)   = 0;
    elseif( phi(k)<0 ) % Wrong hemisphere. Use the antipode
        theta(k) = pi-theta(k);
        phi(k)   = phi(k) + pi;
    elseif( (pi-phi(k))<tol ) % '-x'--'z' plane. Use the antipode
        theta(k) = pi-theta(k);
        phi(k)   = phi(k) - pi;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the transformed matrix G whose samples has been transformed to
% the semisphere 0<=phi<pi
G = [ sin(theta).*cos(phi), sin(theta).*sin(phi), cos(theta) ];
