%%% test_shodf2samples.m
function test_shodf2samples

sf = check_software_platform;
if(sf==2)
    pkg load statistics;
end

clear;
close('all');
%clc;
N  = 10e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SH = [1,0.14,-0.25,0.13,-0.31,0.22];
SH = SH./sqrt(sum(SH.*SH));
SH = sh2squaredsh(reshape(SH,[1,1,1,numel(SH)]));
SH = SH(:)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = size(SH,2);
L = ( sqrt(9-4*(2-2*K)) - 3 )/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c1 = shodf2samples( SH, N );
bds = bound_sh(L,10000);
c2 = shodf2samples( SH, N, bds, [], 1e6 );
fprintf(1,'[Internally computed, brute force] c: %1.9f\n',c1);
fprintf(1,'      [Computed from Y_l^m bounds] c: %1.9f\n',c2);
tic;
[theta,phi,c] = shodf2samples( SH, N, [], [], 1e6 );
dirs = [theta,phi];
fprintf(1,'%1.2g samples generated in %1.2f seconds (c: %1.9f)\n',N,toc,c);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M  = 100;
tt = linspace(0,pi,M);
pp = linspace(0,2*pi,M);
[TT,PP] = meshgrid(tt,pp);
% -----
B   = GenerateSHMatrix(L,[TT(:),PP(:)]);
ODF = B*(SH').*sin(TT(:));
ODF = reshape(ODF,[M,M]);
% -----
pTheta = sum(ODF,1)*(pp(2)-pp(1));
pPhi   = sum(ODF,2)*(tt(2)-tt(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(1));
figure(1);
hold('on');
grid('on');
histogram(dirs(:,1),40,'Normalization','pdf');
plot(tt,pTheta,'LineStyle','--','Color',[.5,.0,.0],'LineWidth',2);
xlabel('theta');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(2));
figure(2);
hold('on');
grid('on');
histogram(dirs(:,2),40,'Normalization','pdf');
plot(pp,pPhi,'LineStyle','--','Color',[.5,.0,.0],'LineWidth',2);
xlabel('phi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close(figure(3));
hf3 = figure(3);
hold('on');
grid('on');
hs = surf(TT,PP,ODF,ones(size(TT)),'EdgeColor','none','FaceAlpha',0.5);
set(hs,'FaceColor','flat');
%hs.CDataMode = 'manual';
set(hs,'CData',ones(size(TT)));
colormap([.5,.0,.0]);
histogram2(dirs(:,1),dirs(:,2),'Normalization','pdf');
xlabel('theta');
ylabel('phi');
rotate3d('on');
title('rotate me');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% -------------------------------------------------------------------------
function bds = bound_sh(L,N)
if(nargin<2), N=1000; end
t   = linspace(-1,1,N);
bds = [];
for l    = 0:2:L
    P    = legendre(l,t);
    bdsl = zeros(1,2*l+1);
    for m=0:1:l
        values = sqrt( (2*l+1)/(2*pi) * factorial(l-m) / factorial(l+m) ) * ...
            P(m+1,:);
        bdsl(m+l+1) = max(abs(values));
    end
    bdsl(l+1) = bdsl(l+1)/sqrt(2);
    if(l>0)
        bdsl(1:l) = bdsl(2*l+1:-1:l+2);
    end
    bds = [bds,bdsl]; %#ok<AGROW>
end
end
