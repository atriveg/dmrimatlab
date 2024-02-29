function [lpar,lperp,nit] = compute_lambdas_from_E(E,b,mu,nu,mlp,nmax,dl,dC,ADC0,chunksz)
% function [lpar,lperp,nit] = compute_lambdas_from_E(E,b,mu,nu,mlp,nmax,dl,dC,ADC0,chunksz)
%
%  Computes the two components of the micro-structural model given the
%  attenuation signal. This is a helper function to atti2micro.
%
%    E: N1 x N2 x N3 x ... x Nd x NB, the signal in the natural domain.
%    b: NB x 1, the set of b-values for each channel in the last dimension
%       of E
%    mu: 1 x 1, penalty weight to avoid lperp being extremely small
%      (default 0.001)
%    nu: 1 x 1, penalty weight to avoid lperp being extremely close to lpar
%      (default 0.001)
%    mlp: 1 x 1, minimum allowed value for lperp (default: 0.01e-3)
%    nmax: 1 x 1, maximum number of Newton-Raphson iterations 
%       (default: 100).
%    dl: maximum allowed chage in lambda from one iteration to the next
%       (default: 1.0e-6).
%    dC: maximum allowed change in the (squared) norm of the vector
%       objective funcion (default: 1.0e-3)
%    ADC0: 1x1, free water diffusivity at human body (default: 3.0e-3)
%    chunksz: 1x1, data are processed chunk-by-chunk (default: 100)
%
%    lpar: N1 x N2 x N3 x ... x Nd, the diffusivity in the parallel
%       direction.
%    lperp: N1 x N2 x N3 x ... x Nd, the diffusivity in the perpendicular
%       direction.
%    nit: N1 x N2 x N3 x ... x Nd, the number of iterations needed to
%       obtain the solution
%

%%% =======================================================================
use_parallel = use_parallel_test;
%%% =======================================================================
if( nargin<3 )
    mu = 0.001;
end
if( nargin<4 )
    nu = 0.001;
end
if( nargin<5 )
    mlp = 0.01;
end
if( nargin<6 )
    nmax = 100;
end
if( nargin<7 )
    dl = 1.0e-6;
end
if( nargin<8 )
    dC = 1.0e-3;
end
if( nargin<9 )
    ADC0 = 3.0e-3;
end
if( nargin<10 )
    chunksz = 10000;
end
%%% =======================================================================
% Reshape to work with arbitrary dimensions:
sz = size(E);
NB   = sz(end);
sz2  = sz(1:end-1);
M    = prod(sz2);
E    = reshape(E,[M,NB]); % M x NB
%%% =======================================================================
if(chunksz<M)
    % -----------
    lpar  = zeros(M,1);
    lperp = zeros(M,1);
    nit   = zeros(M,1);
    % -----------
    if(use_parallel)
        pool    = gcp;
        CHUNKSZ = 5*(pool.NumWorkers)*chunksz;
    else
        CHUNKSZ = 10*chunksz;
    end
    for CK=1:ceil(M/CHUNKSZ)
        % -----------
        IDI = (CK-1)*CHUNKSZ+1;
        IDF = min(CK*CHUNKSZ,M);
        % -----------
        % Piece of data to process in this meta-chunk:
        Ec   = E(IDI:IDF,:); % Meta chunk
        QM   = size(Ec,1);   % Size of this meta-chunk
        lpac = zeros(QM,1);
        lppc = zeros(QM,1);
        nitc = zeros(QM,1);
        % -----------
        if(use_parallel)
            % -- Split the meta-chunk in chunks ad store in cell arrays:
            Ecc  = cell(1,ceil(QM/chunksz));
            for ck=1:ceil(QM/chunksz)
                idi       = (ck-1)*chunksz+1;
                idf       = min(ck*chunksz,QM);
                Ecc{ck}   = Ec(idi:idf,:);
            end
            % -- Process each chunk and store in cell arrays:
            lpacc = cell(1,ceil(QM/chunksz));
            lppcc = cell(1,ceil(QM/chunksz));
            nitcc = cell(1,ceil(QM/chunksz));
            parfor ck=1:ceil(QM/chunksz)
                [lpacc{ck},lppcc{ck},nitcc{ck}] = ...
                    compute_lambdas_from_E_chunk(Ecc{ck},b,mu,nu,mlp,nmax,dl,dC,ADC0);
            end
            % -- Fill the output of the meta-chunk:
            for ck=1:ceil(QM/chunksz)
                idi       = (ck-1)*chunksz+1;
                idf       = min(ck*chunksz,QM);
                lpac(idi:idf) = lpacc{ck};
                lppc(idi:idf) = lppcc{ck};
                nitc(idi:idf) = nitcc{ck};
            end
        else
            % -- Process chunk-by-chunk:
            for ck=1:ceil(QM/chunksz)
                % -----------
                idi       = (ck-1)*chunksz+1;
                idf       = min(ck*chunksz,QM);
                % -----------
                [lpacc,lppcc,nitcc] = ...
                    compute_lambdas_from_E_chunk(Ec(idi:idf,:), ...
                    b,mu,nu,mlp,nmax,dl,dC,ADC0);
                % -----------
                lpac(idi:idf) = lpacc;
                lppc(idi:idf) = lppcc;
                nitc(idi:idf) = nitcc;
                % -----------
            end
        end
        % -----------
        % Fill this meta-chunk:
        lpar(IDI:IDF,1)  = lpac;
        lperp(IDI:IDF,1) = lppc;
        nit(IDI:IDF,1)   = nitc;
        % -----------
    end
    % -----------
else
    [lpar,lperp,nit] = compute_lambdas_from_E_chunk(E,b,mu,nu,mlp,nmax,dl,dC,ADC0);
end
%%% =======================================================================
% Reshape and scale back:
if(numel(sz2)>1)
    lpar  = reshape(lpar,sz2);
    lperp = reshape(lperp,sz2);
    nit   = reshape(nit,sz2);
end
%%% =======================================================================
