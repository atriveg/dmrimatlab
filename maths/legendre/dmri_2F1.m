function i0 = dmri_2F1(g,x)
% function i0 = dmri_2F1(g,x)
%
%    This is an auxiliar function to compute the gaussian hypergeometric
%    function 2F1 in the very special case:
%
%        2F1([1/2,g/2],[3/2],x), with:
%                 0 < g < 2,
%         -infinity < x < 0, (for x = -1/rho,  rho in (-infty,infty) )
%                 0 < x < 1, (for x = 1/kappa, kappa in [1,infty)    )
%
% This implementation is basically a particular case of Andrew D.
% Horchler's hypergeom2F1(a,b,c,z,tol) function:
%
%   Andrew Horchler (2020). 
%   hypergeomq (https://www.github.com/horchler/hypergeomq), 
%   GitHub. Retrieved March 18, 2020.
%   Revision: 1.0, 5-19-14
%
% for a = 1/2, b = g/2, c = 3/2, when abs(x)<1. For abs(x)>1 (basically,
% for 2F1(...,-1/rho)) we use a different approach.
%
% NOTE: this function has been replaced with a much faster, multi-threaded
% mex  function
i0 = dmri_2F1_(g,x);

% i0 = zeros(size(x));
% pp = ( (x<1) & (x>-1) );
% z  = x(pp);
% if(~isempty(z))
%     % -----------------------------------------------------
%     tol     = sqrt(eps)/10000;
%     itermax = 100000;
%     % -----------------------------------------------------
%     h = zeros(size(z));
%     for j = 1:numel(z)
%         Z  = z(j);
%         yi = (1/2)*(g/2)*Z/(3/2);
%         y  = 1+yi;
%         for i = 1:itermax
%             yi = yi*(1/2+i)*(g/2+i)*Z/((i+1)*(3/2+i));
%             y  = y + yi;
%             if(abs(yi)<tol)
%                 break;
%             end
%         end
%         h(j) = y;
%     end
%     % -----------------------------------------------------
%     i0(pp) = h;
% end
% % ------------------------------------------------------------------------
% i0(x>=1) = NaN;
% % ------------------------------------------------------------------------
% x  = -x;
% % ------------------------------------------------------------------------
% pp = ( (x>=1) & (x<=10) );
% if(any(pp))
%     t  = linspace(1,10,100);
%     v  = 1./hypergeom([1/2,g/2],3/2,-t);
%     i0(pp) = 1./interp1(-t,v,-x(pp));
% end
% % ------------------------------------------------------------------------
% pp = ( (x>10) & (x<=100) );
% if(any(pp))
%     t  = linspace(10,100,50);
%     v  = 1./hypergeom([1/2,g/2],3/2,-t);
%     i0(pp) = 1./interp1(-t,v,-x(pp));
% end
% % ------------------------------------------------------------------------
% pp = (x>100);
% if(any(pp))
%     if(abs(g-1)<100*sqrt(eps))
%         arg = sqrt(x(pp));
%         ser = log(4*x(pp))/2./arg;
%         arg = arg.*x(pp);
%         ser = ser + 1/4./arg;
%         arg = arg.*x(pp);
%         ser = ser - 3/32./arg;
%         arg = arg.*x(pp);
%         ser = ser + 5/96./arg;
%         arg = arg.*x(pp);
%         ser = ser - 35/1024./arg;
%         arg = arg.*x(pp);
%         ser = ser + 63/2560./arg;
%     else
%         arg = 1;
%         ser = 1/(1-g); arg = arg.*x(pp);
%         ser = ser + (g/2/(g+1))./arg;  arg = arg.*x(pp);
%         ser = ser - (g*(g+2)/8/(g+3))./arg;  arg = arg.*x(pp);
%         ser = ser + (g*(g+2)*(g+4)/48/(g+5))./arg;  arg = arg.*x(pp);
%         ser = ser - (g*(g+2)*(g+4)*(g+6)/384/(g+7))./arg; arg = arg.*x(pp);
%         ser = ser + (g*(g+2)*(g+4)*(g+6)*(g+8)/3840/(g+9))./arg;
%         ser = ser.*(x(pp).^(-g/2));
%         ser = ser + sqrt(pi)/2*gamma((g-1)/2)/gamma(g/2)./sqrt(x(pp));
%     end
%     i0(pp) = ser;
% end
% % ------------------------------------------------------------------------
% i0(isinf(x)&(x>0)) = 0;
% i0(isinf(x)&(x<0)) = NaN;
% % ------------------------------------------------------------------------
