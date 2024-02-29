function FRT = GenerateFRTMatrix( L )
% L: The order of the SH decomposition (even)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global parameters:
if( L/2 ~= floor(L/2) )
    warning(['The index L has been changed from ',num2str(L),' to ',num2str(L-1)]);
    L = L-1;
end
R   = ( L/2 + 1 )*( L + 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FRT = ones( 1, R );

for l=2:2:L
    FRT( l/2*(l-1)+1:l/2*(l+3)+1 ) = (-1)^(l/2)*prod(1:2:l-1)/prod(2:2:l);
end
FRT = 2*pi*FRT;
FRT = diag( FRT );