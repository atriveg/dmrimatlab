function [lpar,lperp,f,nit] = compute_lambda_par_and_f_from_E(E,b,fmin,fmax,f0,lperp,nmax,dl,dC,ADC0,chunksz)
% function [lpar,lperp,f,nit] = compute_lambda_par_and_f_from_E(E,b,fmin,fmax,f0,lperp,nmax,dl,dC,ADC0,chunksz)
%
%  Computes the components of the micro-structural model given the
%  attenuation signal. This is a helper function to atti2micro. This use
%  case is meant for the case where only two shells are available, so that
%  only two parameters can be inferred: lperp will be kept constant, and
%  the optimization will work over just lpar and f.
%
%    E: N1 x N2 x N3 x ... x Nd x NB, the signal in the natural domain.
%    b: NB x 1, the set of b-values for each channel in the last dimension
%       of E
%    fmin: N1 x N2 x N3 x ... x Nd, the minimum admissible value of the 
%       partial volume fraction of free water
%    fmax: N1 x N2 x N3 x ... x Nd, the minimum admissible value of the 
%       partial volume fraction of free water
%    f0: N1 x N2 x N3 x ... x Nd, the initial temptative value of the 
%       partial volume fraction of free water
%    lperp: N1 x N2 x N3 x ... x Nd, non-optimized value of lperp, which
%       can be different at each voxel. NOTE: no sanity checks will be 
%       performed on lperp to ensure it lies within a physically meaningful
%       range (default: all zeros)
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
%       direction. NOTE: this is the same value passed to the function,
%       since no optimization is performed over lperp.
%    f: N1 x N2 x N3 x ... x Nd, the partial volume fraction of non-free
%       water.
%    nit: N1 x N2 x N3 x ... x Nd, the number of iterations needed to
%       obtain the solution
%

%%% =======================================================================
use_parallel = use_parallel_test;
%%% =======================================================================
% Reshape to work with arbitrary dimensions:
sz   = size(E);
NB   = sz(end);
sz2  = sz(1:end-1);
M    = prod(sz2);
E    = reshape(E,[M,NB]);   % M x NB
fmin = reshape(fmin,[M,1]); % M x 1
fmax = reshape(fmax,[M,1]); % M x 1
f0   = reshape(f0,[M,1]);   % M x 1
%%% =======================================================================
if( nargin<6 )
    lperp = zeros(M,1);
else
    if(isscalar(lperp))
        lperp = lperp*ones(M,1);
    end
end
if( nargin<7 )
    nmax = 100;
end
if( nargin<8 )
    dl = 1.0e-6;
end
if( nargin<9 )
    dC = 1.0e-3;
end
if( nargin<10 )
    ADC0 = 3.0e-3;
end
if( nargin<11 )
    chunksz = 10000;
end

if(chunksz<M)
    % -----------
    lpar  = zeros(M,1);
    f     = zeros(M,1);
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
        Ec     = E(IDI:IDF,:); % Meta chunk
        lperpc = lperp(IDI:IDF,:);
        fmc    = fmin(IDI:IDF);
        fMc    = fmax(IDI:IDF);
        f0c    = f0(IDI:IDF);
        QM     = size(Ec,1);   % Size of this meta-chunk
        lpac   = zeros(QM,1);
        lppc   = zeros(QM,1);
        fc     = zeros(QM,1);
        nitc   = zeros(QM,1);
        % -----------
        if(use_parallel)
            % -- Split the meta-chunk in chunks ad store in cell arrays:
            Ecc     = cell(1,ceil(QM/chunksz));
            lperpcc = cell(1,ceil(QM/chunksz));
            fmcc    = cell(1,ceil(QM/chunksz));
            fMcc    = cell(1,ceil(QM/chunksz));
            f0cc    = cell(1,ceil(QM/chunksz));
            for ck=1:ceil(QM/chunksz)
                idi         = (ck-1)*chunksz+1;
                idf         = min(ck*chunksz,QM);
                Ecc{ck}     = Ec(idi:idf,:);
                lperpcc{ck} = lperpc(idi:idf,:);
                fmcc{ck}    = fmc(idi:idf,:);
                fMcc{ck}    = fMc(idi:idf,:);
                f0cc{ck}    = f0c(idi:idf,:);
            end
            % -- Process each chunk and store in cell arrays:
            lpacc = cell(1,ceil(QM/chunksz));
            lppcc = cell(1,ceil(QM/chunksz));
            fcc   = cell(1,ceil(QM/chunksz));
            nitcc = cell(1,ceil(QM/chunksz));
            parfor ck=1:ceil(QM/chunksz)
                [lpacc{ck},lppcc{ck},fcc{ck},nitcc{ck}] = ...
                    compute_lambda_par_and_f_from_E_chunk(...
                    Ecc{ck}, b, fmcc{ck}, fMcc{ck}, f0cc{ck}, ...
                    lperpcc{ck}, nmax, dl, dC, ADC0 );
            end
            % -- Fill the output of the meta-chunk:
            for ck=1:ceil(QM/chunksz)
                idi       = (ck-1)*chunksz+1;
                idf       = min(ck*chunksz,QM);
                lpac(idi:idf) = lpacc{ck};
                lppc(idi:idf) = lppcc{ck};
                fc(idi:idf)   = fcc{ck};
                nitc(idi:idf) = nitcc{ck};
            end
        else
            % -- Process chunk-by-chunk:
            for ck=1:ceil(QM/chunksz)
                % -----------
                idi       = (ck-1)*chunksz+1;
                idf       = min(ck*chunksz,QM);
                % -----------
                [lpacc,lppcc,fcc,nitcc] = ...
                    compute_lambda_par_and_f_from_E_chunk( Ec(idi:idf,:), b, ...
                    fmc(idi:idf), fMc(idi:idf), f0c(idi:idf), ...
                    lperpc(idi:idf,:), nmax, dl, dC, ADC0 );
                % -----------
                lpac(idi:idf) = lpacc;
                lppc(idi:idf) = lppcc;
                fc(idi:idf)   = fcc;
                nitc(idi:idf) = nitcc;
                % -----------
            end
        end
        % -----------
        % Fill this meta-chunk:
        lpar(IDI:IDF,1)  = lpac;
        lperp(IDI:IDF,1) = lppc;
        f(IDI:IDF,1)     = fc;
        nit(IDI:IDF,1)   = nitc;
        % -----------
    end
    % -----------
else
    [lpar,lperp,f,nit] = ...
        compute_lambda_par_and_f_from_E_chunk( E, b, fmin, fmax, f0, ...
        lperp, nmax, dl, dC, ADC0 );
end
%%% =======================================================================
% Reshape and scale back:
if(numel(sz2)>1)
    lpar  = reshape(lpar,sz2);
    lperp = reshape(lperp,sz2);
    f     = reshape(f,sz2);
    nit   = reshape(nit,sz2);
end
%%% =======================================================================
