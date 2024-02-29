function iota = dmri_compute_Pkernel_integrals(lpar,lperp,gamma,N)
% function iota = dmri_compute_Pkernel_integrals(lpar,lperp,gamma,N)
%
%     Computes the value of the integral
%
%          int_{-1}^{1}x^(2n)/(1-x^2/kappa)^(gamma/2);
%          kappa = lpar/(lpar-lperp)
%
%     for n=0,1,...,N and for arbitrary real gamma>0. Such integrals are
%     used to define the SH coefficient weights to compute the convolutions
%     defining arbitrary order moments.
%
%       INPUTS:
%
%       + lpar: Mx1, parallel diffusivity obtained via the atti2micro.m
%         function.
%       + lperp: Mx1, perpendicular diffusivity obtained via the 
%         atti2micro.m function.
%       These two parameters should fulfill lpar(i)>lperp(i)>0.
%       + gamma: 1x1, real, gamma>0. The parameter of the integral
%       + N: 1x1, non-begative integer. The maximum value of n to be
%         computed, so that the maximum power will be 2N.
%
%       OUPUTS:
%
%       + iota: Mx(N+1), real. Each column corresponds to the integral for
%         each n=0,1,...,N. Each column corresponds to all degrees
%         integrals at each combination of {lpar, lperp}.
%
%     NOTE: no sanity checks are perfomed over any one of the input
%     arguments, so it is the user's responsibility to input appropriate
%     sizes and types.

%%% -----------------------------------------------------------------------
% Check mode depending on gamma:
if(abs(round(gamma)-gamma)>1000*eps)
    gamma_is_integer = false;
    gamma_is_odd     = false;
else
    gamma_is_integer = true;
    if(abs(gamma/2-round(gamma/2))>0.1)
        gamma_is_odd = true;
    else
        gamma_is_odd = false;
    end
end
%%% -----------------------------------------------------------------------
% Initalize the output:
M    = size(lpar,1);
iota = zeros(M,N+1);
%%% -----------------------------------------------------------------------
% Compute kappa:
kappa  = lpar./(lpar-lperp);
ql0    = lpar./lperp;
kappas = sqrt(kappa);
%%% -----------------------------------------------------------------------
% Find the value for n=0:
if(gamma_is_integer)
    if(gamma_is_odd)
        iota(:,1) = 2*kappas.*acsc(kappas);
        i01       = iota(:,1); % For future use
        g0        = 1;
    else
        if(gamma>=1)
            iota(:,1) = 2*kappas.*acoth(kappas);
            g0        = 2;
        else
            iota(:,1) = 2.*ones(size(kappa));
            g0        = 0;
        end
    end
    ql = ql0.^(g0/2-1);
    for g = g0+2:2:gamma
        ql        = ql.*ql0;
        iota(:,1) = (2/(g-2))*ql - ((3-g)/(g-2)).*iota(:,1);
    end
else
    iota(:,1) = 2*dmri_2F1(gamma,1./kappa);
    ql        = ql0.^(gamma/2-1);
end
%%% -----------------------------------------------------------------------
if(N>0)
    if( gamma_is_odd && (gamma>2) )
        %%% ---------- n < gamma/2 - 1/2
        for n=1:round(gamma/2-1/2-1)
            iota(:,n+1) = (-kappa/(n+1/2-gamma/2)).* ...
                ( ql - (n-1/2).*iota(:,n) );
        end
        %%% ---------- n = gamma/2 - 1/2
        kt  = 1./(kappa-1);
        kt0 = 1./sqrt(kt);
        for k=3:2:gamma
            kappas = kappas.*kappa;
            kt0    = kt0.*kt;
            i01    = (2*kappas/(k-2)).*kt0 - kappa.*i01;
        end
        iota( :, round(gamma/2-1/2)+1 ) = i01;
        %%% ---------- n > gamma/2 - 1/2
        for n=round(gamma/2-1/2+1):N
            iota(:,n+1) = (-kappa/(n+1/2-gamma/2)).* ...
                ( ql - (n-1/2).*iota(:,n) );
        end
    else
        for n=1:N
            iota(:,n+1) = (-kappa/(n+1/2-gamma/2)).* ...
                ( ql - (n-1/2).*iota(:,n) );
        end
    end
end
%%% -----------------------------------------------------------------------
% Fix the very special case when lpar=lperp and kappa -> inf, then the
% integral reduces to:
%       int_{-1}^{1}x^(2n)/(1-x^2/kappa)^(gamma/2);
%    -> int_{-1}^{1}x^(2n) = 1/(n+1/2)
for n=0:N
    pp           = isinf(kappa); % M x 1
    iota(pp,n+1) = 1/(n+1/2);    % M x (N+1)
end
%%% -----------------------------------------------------------------------
