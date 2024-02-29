function SHp = dipy2dmrilab_SH_coeffs(SH)
% function SHp = dipy2dmrilab_SH_coeffs(SH)
%
%   Converts SH coefficients computed with Python's dipy package to SH
%   coefficients computed with this package. The conversion works back and 
%   forth since it is the same in both directions (a change in trhe sign of
%   odd, negative values of index m).
%
% INPUTS:
%
%   SH: an N-D array containing the SH coefficients to convert. The actual
%      indexing of the SH coefficients is assumed to correspond to the last
%      dimension of the array, that should therefore have size 6, 15, 28, 
%      45...
%
% OUTPUTS:
%
%   SHp: The converted coefficients with the new convention

% Check the mandatory input argments:
if(nargin<1)
    error('At least the SH coefficients volume must be supplied');
end

D  = ndims(SH);
SZ = size(SH);
LL = SZ(D);
L = (sqrt(8*(LL-1)+9)-3)/2; % SH order
if( abs(L-round(L))>1.0e-9 || L<2 )
    error('Weird size of the SH volume. Its fourth dimension should have size 6, 15, 28, 45, ..., (L+1)(L+2)/2, with L=2,4,6,...');
end

% Reshape the volume to work with arbitrary directions:
sz  = [ prod(SZ(1:D-1)), LL ];
SHp = reshape(SH,sz);

% Recover the m indices:
mm = zeros(1,(L+1)*(L+2)/2);
pp = 2;
for l=2:2:L
    nl = 2*l+1;
    mm(1,pp:pp+nl-1) = -l:l;
    pp = pp + nl;
end

% Which m indices have to be inverted?
pt = (mm<0) & ( abs((-mm/2)-round(-mm/2))>0.25 );

% Invert them:
SHp(:,pt) = -SHp(:,pt);

% Reshape back to the original size of the array:
SHp = reshape(SHp,SZ);
