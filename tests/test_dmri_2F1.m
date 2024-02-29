% test_dmri_2F1

x     = -[-0.999,-0.9,-0.5,0.12,0.77,1.0,1.13,5.2,13.4,24.7,51.7,99.2,200.3,400.5,800.3];
gamma = [0.11,0.33,0.67,1.33,1.66,1.78,5];

close(figure(1001));
figure(1001);
xlabel('x');
ylabel('2F1([1/2,\gamma/2],3/2,x)');
hold('on');
cols = [1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;1,0,1;0.6,0.6,0.3];
grid('on');
hl = zeros(1,length(gamma));
hltext = cell(1,length(gamma));
for n=1:length(gamma)
    f1 = hypergeom([1/2,gamma(n)/2],3/2,x);
    f2 = dmri_2F1(gamma(n),x);
    hl(n) = plot(x,f1,'Marker','o','MarkerSize',20,'LineStyle','none','Color',cols(n,:));
    hltext{n} = ['gamma = ',num2str(gamma(n))];
    plot(x,f2,'Marker','*','MarkerSize',20,'LineStyle','none','Color',cols(n,:));
end
legend(hl,hltext{:});