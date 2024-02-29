%%% -----------------------------------------------------------------------
function [bs,ps,Ns] = auto_detect_shells(bi,bth)
% bi:  Gx1, acquired shells "as they are"
% bth: 1x1, threshold to tell if bi is equal or not to bj
% bs:  Nsx1, the (averaged) b-value assigned to each shell
% ps:  Gx1, tells which shell (from 1 to Ns) each bi belongs to
% Ns:  1x1, the total number of shells

[bi2,pt] = sort(bi,'ascend');      % Gx1, bi2 = bi(pt)
bid      = [bth+1;diff(bi2)];      % Gx1
change   = (bid>bth+100*eps);      % Gx1
change   = cumsum(double(change)); % Gx1, equivalent bs labelled together
pti      = (1:length(bi))';        % Gx1 
pti      = pti(pt);                % Gx1
ps       = zeros(size(bi));        % Gx1
ps(pti)  = change;                 % Gx1, shell order of each bi
Ns       = max(change);            % 1x1, number of shells
bs       = zeros(Ns,1);            % Gx1, effective b-value for each shell
for n=1:Ns
    bs(n) = mean( bi(abs(ps-n)<1/2) );
end