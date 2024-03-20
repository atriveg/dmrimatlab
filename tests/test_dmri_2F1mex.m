% test_dmri_2F1mex
clear;
x     = -[-1,-0.99999999,-0.9999999,-0.999999,-0.99999,-0.9999,-0.999,-0.9,-0.5,0.12,0.77,0.999,0.999999,0.999999999,1.0,1.000000001,1.000001,1.001,1.13,5.2,13.4,24.7,51.7,99.2,200.3,400.5,800.3,1200.7,2000.9];
% Should not use interger values of gamma 2 or greater, since otherwise the
% (necessary) computation of Gamma(3/2-gamma/2) and Gamma(1-gamma/2) will 
% fail because the argument is a negative integer/zero
%
% Note this is an issue within the toolbox, since these gamma values are
% are treated without calling the dmri_2F1 function, in a recursive
% fashion.
%
% What might be an issue (check in the plots) is using g>=2 for x->1.
% However, x->1 should be a degeneracy, and any INTEGER value of g>=2
% will be correctly managed with recursive rules.
gamma = [-1.0,-0.5,0.11,0.33,0.67,0.99,1.0,1.001,1.33,1.66,1.78,3.1,4.9,6.1];

close(figure(1001));
figure(1001);
xlabel('x');
ylabel('2F1([1/2,\gamma/2],3/2,x)');
hold('on');
cols = jet( length(gamma) );
grid('on');
hl = zeros(1,length(gamma));
hltext = cell(1,length(gamma));
for n=1:length(gamma)
    f1 = hypergeom([1/2,gamma(n)/2],3/2,x);
    if(length(f1)~=length(x))
        % hypergeom fails at x=1 for certain values of gamma, and it only
        % returns a single "inf" value crashing the test. Correct this
        % missbehavior:
        f1 = [ nan, hypergeom([1/2,gamma(n)/2],3/2,x(2:end)) ];
    end
    f2 = dmri_2F1mex( ...
        gamma(n), ...
        x );
    hl(n) = plot(x,f1,'Marker','o','MarkerSize',20,'LineStyle','none','Color',cols(n,:));
    hltext{n} = ['gamma = ',num2str(gamma(n))];
    plot(x,f2,'Marker','*','MarkerSize',20,'LineStyle','none','Color',cols(n,:));
end
legend(hl,hltext{:});
