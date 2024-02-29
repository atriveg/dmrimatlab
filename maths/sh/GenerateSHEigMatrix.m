function out = GenerateSHEigMatrix( L )
% function out = GenerateSHEigMatrix( L )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global parameters:
if( L/2 ~= floor(L/2) )
    warning(['The index L has been changed from ',num2str(L),' to ',num2str(L-1)]);
    L = L-1;
end
R   = ( L/2 + 1 )*( L + 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out = zeros( 1, R );
for l=2:2:L
    out( l/2*(l-1)+1:l/2*(l+3)+1 ) = -l*(l+1);
end
out = diag( out );