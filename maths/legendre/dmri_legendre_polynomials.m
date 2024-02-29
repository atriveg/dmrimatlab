function p = dmri_legendre_polynomials(L)
% function p = dmri_legendre_polynomials(L)
%
%    Returns the coefficients of the first L even-order non-associated
%    Legendre polynomials. Therefore, p contains only even-power
%    coefficients. NOTE: no sanity check is performed over L, so it is the
%    responsibility of the user to input an even integer g.e.t. zero.
%
%      L: 1x1 even integer
%      p: (L/2+1)x(L/2+1) double:
%          - each row is a Legendre polynomial, from l=0,2,...,L
%          - each column is an even power, from x^0,x^2,...,x^L
p = zeros(L/2+1,L/2+1);
if(L<1)
    p = 1;
else
    % Recursively compute by using:
    %   n*P_{n}(x) = (2n-1)*x*P_{n-1}(x) - (n-1)*P_{n-2}(x)
    p0     = [1, zeros(1,L/2)]; % Initalize: P_0(x) = 1
    p1     = [1, zeros(1,L/2)]; % Initalize: P_1(x) = x
    p(1,:) = p0;                % First even polynomial
    n      = 2;                 % Order of the polynomial to compute
    for l=1:L/2
        % Compute the next even polynomial:
        p0 = ( (2*n-1)*[0,p1(1:end-1)]- (n-1)*p0 )/n;
        n  = n+1; % Next polynomial is odd
        % Compute the next odd polynomial:
        p1 = ( (2*n-1)*p0- (n-1)*p1 )/n; % Last one becomes cropped, but it won't be used anyway
        n  = n+1; % Next polynomial is even
        % Assign the next even polynomial
        p(l+1,:) = p0;
    end
end
