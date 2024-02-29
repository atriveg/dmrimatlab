function lambda = dmri_gcv(Phi,L,y,lmax,lmin,Nl)
% Phi: M x N
% L:   N x N
% y:   M x 1
% lmax: scalar
% lmin: scalar
% Nl: scalar

lstep = exp(log(lmin/lmax)/Nl);
[M,N] = size(Phi);
l     = exp( log(lmax) + (0:Nl)*log(lstep) );

Q0  = inf;
IM  = eye(M);
IN  = eye(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Nl+1
    Al = (Phi')*Phi + l(n)*((L')*L);
    Sl = Phi*(Al\IN)*Phi';
    den = M - trace(Sl);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    el = ((Sl-IM)*y);
    el = el/den;
    Q  = (el')*el;
    if(Q>Q0)
        break;
    else
        Q0 = Q;
    end
end

if(n==1)
    lambda = nan;
else
    if(n<Nl)
        lambda = l(n-1);
    else
        lambda = l(n);
    end
end

end



































