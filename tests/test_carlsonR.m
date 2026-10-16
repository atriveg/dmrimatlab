%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_carlsonR.m
rng(13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z = 1;
y = [1/16,1/4,1,4,16,64,256,1024];
x = 10.^(linspace(-1,3,1000));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lg = {};
for n=1:length(y)
    lg{n} = sprintf('y = %1.2f',y(n));
end
colors = rand(length(y),3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(1));
hf = figure(1);
hold('on');
grid('on');
hl = zeros(1,length(y));
for n=1:length(y)
    hl(n) = plot( x, carlsonRd(x,y(n),z), 'Color', colors(n,:), 'LineWidth', 2 );
end
xlabel('x','FontSize',20);
ylabel('R_d(x,y,1)','FontSize',20);
legend(hl,lg{:},'Location','NorthEast');
set( get(hf,'CurrentAxes'), 'XScale', 'log' );
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(2));
hf = figure(2);
hold('on');
grid('on');
hl = zeros(1,length(y));
for n=1:length(y)
    hl(n) = plot( x, carlsonRf(x,y(n),z), 'Color', colors(n,:), 'LineWidth', 2 );
end
xlabel('x','FontSize',20);
ylabel('R_f(x,y,1)','FontSize',20);
legend(hl,lg{:},'Location','NorthEast');
set( get(hf,'CurrentAxes'), 'XScale', 'log' );
drawnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
