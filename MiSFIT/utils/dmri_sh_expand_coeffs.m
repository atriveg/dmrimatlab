function ptr  = dmri_sh_expand_coeffs( L )
% function ptr  = dmri_sh_expand_coeffs( L )
%
%   This is a simple helper function to expand SH coefficients that are
%   defined only for {l,m=0}, i.e., C_0, C_2, C_4 to full SH coefficients
%   for all the range: C_0, C_2, C_2, C_2, C_2, C_2, C_4, C_4, C_4, C_4,
%    C_4, C_4, C_4, C_4, C_4...
ptr = ones(1,(L+1)*(L+2)/2);
pos = 1;
for l=2:2:L
    nl  = 2*l+1;
    ptr(1,pos+1:pos+nl) = l/2+1;
    pos = pos +nl;
end
