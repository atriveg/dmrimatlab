function w = carlsonRd( x, y, z, maxthreads )
% function w = carlsonRd( x, y, z[, maxthreads] )
%
%   Computes the Carlson elliptic integral of the second kind,
%   Rd, at each element of x, y and z, which must be either
%   scalars or the same size as the others. The integral is
%   defined as:
%
%     Rd(x,y,z) = (3/2)*int_0^{infty} ( (t+x)(t+y)(t+z)^3 )^{-1/2} dt
%
%   This is a re-implementation of Milan Batista's rd function in:
%
%      https://es.mathworks.com/matlabcentral/fileexchange/65915-elfun18/files/elFun18v1_5/elfun18/rd.m
%
%   that has been coded as a parallel mex file in C++. The fourth
%   (optional) argument is the maximum number of threads to be
%   used for computations (though this value can be reduced if it
%   exceeds the number of logical cores in the machine).
%
%   NOTE: each of x, y and z must be real-valued. They all must be
%   the same size, unless one (or several) of them is scalar.
%   Besides, corresponding entries must fulfill:
%
%     x,y: must be greater or equal than 0, but only one of them
%          can exactly equal 0.
%       z: must be strictly greater than 0.
%
%   References (taken from the original implementation):
%
%   [1] B.C.Carlson,"Numerical computation of real or complex elliptic
%       integrals", Numerical Algorithms 10, 1995, pp 13-26.
%   [2] B.C.Carlson, "Elliptic Integrals" In: F.W.J.Olver et al., NIST
%       Handbook of Mathematical Functions, Chap.19, Cambridge, 2010.
%       Handbook of Mathematical Functions, Cambridge, 2010.
%   [3] W.H.Press et al., "Numerical Recipes in Fortran 90", 2nd ed.,
%       Cambridge, 1996
%
%   This is a helper function to dkt2kurtosis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check inputs:
if(nargin<4)
    maxthreads = 1e6;
else
    assert( isscalar(maxthreads), 'maxthreads must be scalar');
    maxthreads = round(maxthreads);
    assert( maxthreads>0, 'maxthreads must be a positive integer');
end
% --------
if( isempty(x) || isempty(y) || isempty(z) )
    w = [];
    return;
end
% --------
sz = [1,1];
if( ~isscalar(x) )
    sz = size(x);
end
if( ~isscalar(y) )
    if(isscalar(sz))
        sz = size(y);
    else
        assert( isequal(size(y),sz), 'x, y, z must be either scalars or have the same size' );
    end
end
if( ~isscalar(z) )
    if(isscalar(sz))
        sz = size(z);
    else
        assert( isequal(size(z),sz), 'x, y, z must be either scalars or have the same size' );
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process:
x = x(:);
y = y(:);
z = z(:);
% --------
[w,~] = carlsonR_( double(x), double(y), double(z), ...
    'd', maxthreads );
% --------
w = reshape(w,sz);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
