function [y0,y1,y2] = laguerreL_1_2(x)
% function [y0,y1,y2] = laguerreL_1_2(x)
%
%    Computes the Laguerre function associated to Rician noise
%    distributions and its first two derivatives:
%
%       y0(x) = sqrt(pi/2)*laguerrel(1/2,-x^2/2);
%       y1(x) = d y0(x) / dx
%       y2(x) = d^2 y0(x) / dx^2
%
%    It makes use of several equalities and assymptotic equations to avoid
%    the long computation times and stability problems of Matlabs's
%    built-in laguerreL

% --------------------
x  = double(x);   % Bessel function implementations seem to work bad for singles
pb = (abs(x)<50); % Beyond x=50, bessel functions become NaN.
pa = (~pb);
% --------------------
y0 = zeros(size(x));
y1 = zeros(size(x));
y2 = zeros(size(x));
if(any(pb(:)))
    % --------------------
    x2 = x(pb).*x(pb);
    i0 = besseli(0,x2/4);
    i1 = besseli(1,x2/4);
    ef = sqrt(pi/2)*exp(-x2/4)/2;
    % --------------------
    y0(pb) = ef.*(x2.*(i1+i0) + 2*i0);
    % --------------------
    y1(pb) = ef.*x(pb).*(i0+i1);
    % --------------------
    y2(pb) = ef.*(i0-i1);
    % --------------------
end
if(any(pa(:)))
    % --------------------
    x   = x(pa);
    xi1 = 1./x;
    xi2 = xi1.*xi1;
    xi3 = xi1.*xi2;
    xi4 = xi2.*xi2;
    xi5 = xi3.*xi2;
    xi6 = xi3.*xi3;
    xi7 = xi3.*xi4;
    % --------------------
    y0(pa) = x + xi1/2 + xi3/8 + 3*xi5/16;
    % --------------------
    y1(pa) = 1 - xi2/2 - 3*xi4/8 - 15*xi6/16;
    % --------------------
    y2(pa) = xi3 + 3*xi5/2 + 45*xi7/8;
    % --------------------
end
