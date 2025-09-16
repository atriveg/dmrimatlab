%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test_mexAllSquaredSHFactors.m
% Change Lmax to change how far to go (may become slow!)
% Prints only the number of non-null factors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_mexAllSquaredSHFactors
clear;
close('all');
clc;
Lmax = 20;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hmsg1 = 'L';        L1 = length(hmsg1);
hmsg2 = 'KxKxKp=N'; L2 = length(hmsg2);
hmsg3 = 'P';        L3 = length(hmsg3);
hmsg4 = 'sparsity'; L4 = length(hmsg4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg1 = cell(1,Lmax/2+1);
msg2 = cell(1,Lmax/2+1);
msg3 = cell(1,Lmax/2+1);
msg4 = cell(1,Lmax/2+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for L=0:2:Lmax
    % ---
    factors = mexAllSquaredSHFactors(L);
    % ---
    P  = length(factors);
    K  = (L+1)*(L+2)/2;
    Kp = (2*L+1)*(2*L+2)/2;
    N  = K*K*Kp;
    % ---
    msg1{L/2+1} = sprintf('%i',L);                 L1 = max(L1,length(msg1{L/2+1}));
    msg2{L/2+1} = sprintf('%ix%ix%i=%i',K,K,Kp,N); L2 = max(L2,length(msg2{L/2+1}));
    msg3{L/2+1} = sprintf('%i',P);                 L3 = max(L3,length(msg3{L/2+1}));
    msg4{L/2+1} = sprintf('%1.2f',100*P/N);        L4 = max(L4,length(msg4{L/2+1}));
    % ---
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pad_n_print(hmsg1,L1);
pad_n_print(hmsg2,L2);
pad_n_print(hmsg3,L3);
pad_n_print(hmsg4,L4);
fprintf(1,'\n\n');
for L=0:2:Lmax
    pad_n_print(msg1{L/2+1},L1);
    pad_n_print(msg2{L/2+1},L2);
    pad_n_print(msg3{L/2+1},L3);
    pad_n_print(msg4{L/2+1},L4,true);
    fprintf(1,'\n');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pad_n_print(msg,L,percent)
if(nargin<3)
    percent = false;
end
pad = 3;
if(percent)
    fprintf(1, '%s%s%%', repmat(' ',[1,L+pad-length(msg)]), msg );
else
    fprintf(1, '%s%s', repmat(' ',[1,L+pad+1-length(msg)]), msg );
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%














