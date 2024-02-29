function B = GenerateSHMatrix( L, G )
% function B = GenerateSHMatrix( L, G )
% L: The maximum order of the SH expansion (even)
% G: May be:
%            a) A matrix with each gradient direction [ g_xi, g_yi, g_zi ]
%               with g_(x,y,z)i column vectors
%            b) A matrix with the angular coordinates: [ theta_i, phi_i ];
%               physics convention, i.e. 0 <= theta_i < pi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global parameters:
if( L/2 ~= floor(L/2) )
    warning(['The index L has been changed from ',num2str(L),' to ',num2str(L-1)]);
    L = L-1;
end
N = size(G,1);
R = ( L/2 + 1 )*( L + 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the spherical coordinates:
if( size(G,2)>2 )
    [~,theta,phi] = HARDIG( G );
else
    theta     = G(:,1);
    phi       = G(:,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the matrix
B = zeros( N, R );
% The 0-th order SH is trivial:
B(:,1) = 1/2/sqrt(pi);
% For each order of the SH expansion:
for l=2:2:L
    AUX = zeros( N, 2*l+1 );
    % Positive m:
    AUX( :, l+1:end ) = legendre( l, cos(theta) )';
    % Negative m:
    for m=-1:-1:-l
        AUX( :, m+l+1 ) = ( (-1)^m * factorial(l+m)/factorial(l-m) ).*AUX( :, l-m+1 );
    end
    % Compute the part that depends on phi (on the modified, real basis) and normalize:
    for m=-l:-1
        AUX( :, m+l+1 ) = sqrt( (2*l+1)/(2*pi) * factorial(l-m)/factorial(l+m) ).*AUX( :, m+l+1 ).*cos( m.*phi );
    end
    AUX( :, l+1 ) = sqrt( (2*l+1)/(4*pi) ).*AUX( :, l+1 );
    for m=1:l
        AUX( :, m+l+1 ) = sqrt( (2*l+1)/(2*pi) * factorial(l-m)/factorial(l+m) ).*AUX( :, m+l+1 ).*sin( m.*phi );
    end
    B( :, l/2*(l-1)+1:l/2*(l+3)+1 ) = AUX;  
end

% Done
