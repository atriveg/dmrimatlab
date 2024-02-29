function shc = compute_shodf_from_micro(atti,lpar,lperp,bs,gi,p,L,lambda,optimal,chunksz)
% function shc = compute_shodf_from_micro(atti,lpar,lperp,bs,gi,p,L,lambda,optimal,chunksz)
%
%   Computes the SH coefficients of the ODF provided a micro-structural
%   model. This is a helper function to micro2shodf.
%
%     atti:    Q x G_n
%     lpar:    Q x 1
%     lperp:   Q x 1
%     bs:      Ns x 1
%     gi:      G x 3
%     p:       1 x NS cell
%     L:       1 x 1
%     lambda:  1 x 1
%     optimal: 1 x 1 boolean

% -------------------------------------------------------------------------
is_broadcast_available = is_broadcast_available_test;
use_parallel           = use_parallel_test;
% -------------------------------------------------------------------------
% Compute the convolution coefficients for the Funk-Hecke theorem:
F   = dmri_compute_convolution_weights_ODF(bs,lpar,lperp,L); % Qx(L/2+1)xNs
ptr = dmri_sh_expand_coeffs(L); % 1 x (L+1)(L+2)/2, will be used to expand F
ptr = ptr(2:end);               % 1 x (L+1)(L+2)/2-1, skip DC component
Ns  = size(bs,1);
% -------------------------------------------------------------------------
% Remove the DC component to the atti at each shell:
for n=1:Ns
    if(is_broadcast_available)
        atti(:,p{n}) = atti(:,p{n}) - F(:,1,n)/(4*pi);                       % QxG_n - Qx1 -> QxG_n
    else
        atti(:,p{n}) = bsxfun( @(x,y)(x-y), atti(:,p{n}), F(:,1,n)/(4*pi) ); % QxG_n - Qx1 -> QxG_n
    end
end
% -------------------------------------------------------------------------
% Initialize SH maths for each shell
B   = cell(1,Ns);
LR  = GenerateSHEigMatrix( L ); % KxK, where K=(L+1)(L+2)/2
LR  = LR(2:end,2:end);          % (K-1)x(K-1)
LR  = LR.*LR;                   % (K-1)x(K-1)
% ----------
% Cell arrray WLS as precomputed here is only used if the 'suboptimal'
% flag is active, otherwise it will be overwritten. In case the
% estimation is subotimal, besides, the HA expansion concerns E(q) and
% not the ODF, so that some extra regularization is needed so that the
% high-pass filter difision by F doesn't blow the noise. The penalty
% has to be higher for smaller values of b, since the lower b ther
% stronger the high-pass filtering
WLS = cell(1,Ns);   % Used only if optimal==false
l1t = 2.7e-3;       % Typical value for lpar
l2t = 0.2e-3;       % Typical valie for lperp
Ft  = dmri_compute_convolution_weights_ODF(bs,l1t,l2t,L); % 1 x (L/2+1) x Ns
Ft  = 1./Ft(1,ptr,:);                                     % 1 x ((L+1)(L+2)/2-1) x Ns
Ft  = Ft.*Ft;                                             % 1 x ((L+1)(L+2)/2-1) x Ns
% ----------
for n=1:Ns
    Bn     = GenerateSHMatrix( L, gi(p{n},:) ); % G_nxK, where K=(L+1)(L+2)/2
    Bn     = Bn(:,2:end);                       % G_nx(K-1), remove DC component
    B{n}   = Bn;                                % G_nx(K-1)
    % ----------------------
    LR2 = LR.*diag(Ft(1,:,n));
    WLS{n} = (Bn'*Bn+0.25*lambda.*LR2)\(Bn');   % ((K-1)x(K-1))^(-1)*((K-1)xG_n) -> (K-1)xG_n
    WLS{n} = WLS{n}';                           % G_nx(K-1), for future convenience
end
% -------------------------------------------------------------------------
% We can now proceed with the estimation:
Q    = size(atti,1);
K    = (L+1)*(L+2)/2;
shc  = zeros(Q,K-1);
if(optimal)
    % No other way than for-looping through the data. This is the slowest 
    % part of the algorithm, so we will use parallel matlab workers 
    % if available:
    if(use_parallel)
        parfor q=1:Q % At each imaged voxel...
            % Compute the matrixes involved in the LS problem:
            [WLS,PRE] = LS_compute_matrixes( F(q,:,:), atti(q,:), B, p, LR, lambda, ptr, K );
            % Invert:
            shc(q,:) = (WLS\PRE)'; % 1x(K-1)
        end
    else % No parallel workers
        for q=1:Q % At each imaged voxel...
            % Compute the matrixes involved in the LS problem:
            [WLS,PRE] = LS_compute_matrixes( F(q,:,:), atti(q,:), B, p, LR, lambda, ptr, K );
            % Invert:
            shc(q,:) = (WLS\PRE)'; % 1x(K-1)
        end
    end
else
    % Work chunk-by-chunk
    if(use_parallel)
        shcc  = cell(1,ceil(Q/chunksz));
        Fc    = cell(1,ceil(Q/chunksz));
        attic = cell(1,ceil(Q/chunksz));
        for ck=1:ceil(Q/chunksz)
            idi       = (ck-1)*chunksz+1;
            idf       = min(ck*chunksz,Q);
            Fc{ck}    = F(idi:idf,:,:);
            attic{ck} = atti(idi:idf,:);
        end
        parfor ck=1:ceil(Q/chunksz)
            % ------------------------------------------
            % Characterize the present chunk:
            idi = (ck-1)*chunksz+1;
            idf = min(ck*chunksz,Q);
            Qc  = idf-idi+1;
            % ------------------------------------------
            sh  = compute_odf_suboptimal( Fc{ck}, attic{ck}, ...
                B, p, WLS, ptr, Qc, K, Ns );
            % ------------------------------------------
            shcc{ck} = sh;
        end
        for ck=1:ceil(Q/chunksz)
            idi = (ck-1)*chunksz+1;
            idf = min(ck*chunksz,Q);
            shc(idi:idf,:) = shcc{ck};
        end
    else
        for ck=1:ceil(Q/chunksz)
            % ------------------------------------------
            % Characterize the present chunk:
            idi = (ck-1)*chunksz+1;
            idf = min(ck*chunksz,Q);
            Qc  = idf-idi+1;
            % ------------------------------------------
            sh  = compute_odf_suboptimal( F(idi:idf,:,:), atti(idi:idf,:), ...
                B, p, WLS, ptr, Qc, K, Ns );
            % ------------------------------------------
            shc(idi:idf,:) = sh;
        end
    end
end

% -------------------------------------------------------------------------
function [WLS,PRE] = LS_compute_matrixes( F, atti, B, p, LR, lambda, ptr, K )
% F:      1 x (L/2+1) x Ns
% atti:   1 x G
% B:      1 x Ns cell, G_n x (K-1) each
% p:      1 x Ns cell, 1 x G_n each
% LR:     (K-1) x (K-1)
% lambda: 1x1
% ptr:    1 x K-1
% K:      1 x 1
% ------------------------------------------
% Intialize sum:
WLS = LR*lambda;       % (K-1)x(K-1)
PRE = zeros(K-1,1);    % (K-1)x1
% ------------------------------------------
% Sum...
for n=1:size(F,3)
    Fn  = F(1,ptr,n);                  % 1x(K-1)
    Bn  = (B{n})*diag(Fn);             % G_nx(K-1)
    WLS = WLS + (Bn')*Bn;              % (K-1)x(K-1)
    PRE = PRE + (Bn')*(atti(1,p{n})'); % (K-1)x1
end

% -------------------------------------------------------------------------
function shc = compute_odf_suboptimal( F, atti, B, p, WLS, ptr, Qc, K, Ns )
is_broadcast_available = is_broadcast_available_test;
% ------------------------------------------
sh  = zeros(Qc,K-1,Ns); % Qc x K-1 x Ns, SH coefficients at each shell
% Compute the ODF suggested by each shell:
for n=1:Ns
    sh(:,:,n) = atti(:,p{n})*WLS{n};   % Qc x (K-1), sh of E(q)
    sh(:,:,n) = sh(:,:,n)./F(:,ptr,n); % Qc x (K-1), sh of Phi(u)
end
% ------------------------------------------
% Compute the global error committed by the estimation of Phi(u) at
% each shell
wgs = zeros(Qc,Ns);  % Qc x Ns
for n=1:Ns
    err = zeros(Qc,1);
    for m=1:Ns
        rec = sh(:,:,n).*F(:,ptr,m); % Qc x (K-1)
        rec = rec*(B{m}');           % ( Qc x (K-1) ) x ( G_n x (K-1) )' -> Qc x G_n
        rec = rec - atti(:,p{m});    % Qc x G_n
        err = err + sum(rec.*rec,2); % Qc x 1
    end
    wgs(:,n) = err;
end
% ------------------------------------------
% Design weights depending on the errors committed:
% wgs has size Qc x Ns
wgs = 1./(wgs+1000*sqrt(eps));
wgn = sum(wgs,2);
% Normalize weights so they sum up to 1
if(is_broadcast_available)
    wgs = wgs./wgn;                     % Qc x Ns
else
    wgs = bsxfun(@(x,y)(x./y),wgs,wgn); % Qc x Ns
end
% ------------------------------------------
% Weighted average of the SH coefficients obtained at each shell
% for the ODF depending on the errors they produce:
shc = zeros(Qc,K-1);
for n=1:Ns
    if(is_broadcast_available)
        shc = shc + wgs(:,n).*sh(:,:,n);                           % Qc x (K-1)
    else
        shc = shc + bsxfun( @(x,y)(x.*y), wgs(:,n), sh(:,:,n)   ); % Qc x (K-1)
    end
end
