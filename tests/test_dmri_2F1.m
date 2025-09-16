% test_dmri_2F1
sf = check_software_platform;
if(sf==2)
    pkg load gsl;
end


x     = -[ -0.999,-0.9,-0.5, [1.0,1.13,5.2,13.4,24.7,51.7,99.2,200.3,400.5,800.3]/1000 ];
gamma = [0.11,0.33,0.67,1.0,1.33,1.66,1.78,5];

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
    if(sf==1)
        f1 = hypergeom([1/2,gamma(n)/2],3/2,x);
        if(length(f1)~=length(x))
            % hypergeom fails at x=1 for certain values of gamma, and it only
            % returns a single "inf" value crashing the test. Correct this
            % missbehavior:
            f1 = [ nan, hypergeom([1/2,gamma(n)/2],3/2,x(2:end)) ];
        end
    else
        f1 = gsl_sf_hyperg_2F1(1/2,gamma(n)/2,3/2,x);
    end
    f2 = dmri_2F1(gamma(n),x);
    hl(n) = plot(x,f1,'Marker','o','MarkerSize',20,'LineStyle','none','Color',cols(n,:));
    hltext{n} = ['gamma = ',num2str(gamma(n))];
    plot(x,f2,'Marker','*','MarkerSize',20,'LineStyle','none','Color',cols(n,:));
end
legend(hl,hltext{:});
