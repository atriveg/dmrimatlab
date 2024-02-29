function atti = micro2atti(sh,lpar,lperp,f,gi,bi,varargin)
% function atti = micro2atti(sh,lpar,lperp,f,gi,bi, ...
%                                           'opt1',value1,'opt2',value2...)
%
%   Reconstructs the attenuation signal E(q) (without the free water
%   compartment) provided a microstructure model described by lpar and
%   lperp and an ODF described by its SH coefficients. E(q) is evaluated at
%   a multi-shell sampling with directions gi and b-values bi.
%
%       sh: MxNxPx((L+1)(L+2)/2), a double array with the SH coefficients 
%           of the ODF at each voxel up to order L. Should be computed with
%           micro2shodf.
%       lpar: MxNxP, a double array with the parallel diffusivity at each
%           voxel, estimated with atti2micro.
%       lperp: MxNxP, a double array with the traverse diffusivity at each
%           voxel, estimated with atti2micro.
%       f: a MxNxP double array with the partial volume fraction of
%           intra-cellular water (should fulfill 0<=f<=1). If an empty
%           array is passed, then f=1 for all voxels is assumed, so that
%           ones(M,N,P) has the same effect as [].
%       gi: Gx3, the gradients table, each row must be a unit-norm
%           direction.
%       bi: Gx1, the b-values describing a multi-shell sampling. This means
%           bi should have few unique values (typically 2 to 4):
%                length(unique(bi)) << length(bi).
%
%       atti: MxNxPxG, the attenuation signal reconstructed at the desired
%           points.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style. General parameters:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      bth: 1x1, b-values that differ from each other less than this
%         threshold are assumed to belong to the same shell (default: 100).
%      ADC0: estimated diffusivity of free water at body temperature 
%         (Diso). Should use the same as in atti2micro (default: 3.0e-3).
%      chunksz: The processing is done in blocks of size chunksz to avoid
%         blowing up the memory (default: 100).

%%% ----------------
[M,N,P] = size(lpar);
assert(isequal(size(lpar),size(lperp)),'lpar and lperp must be the same size');
if(~isempty(f))
    assert(isequal(size(f),size(lpar)),'lpar and f must be the same size');
end
%%% ----------------
[M2,N2,P2,K] = size(sh);
assert(isequal([M2,N2,P2],[M,N,P]),'The first 3 dimensions of sh must match those of lpar and lperp');
L = (-3+sqrt(1+8*K))/2;
assert( abs(round(L)-L)<1000*eps, 'This is a weird size for the fourth dimension of sh; make sure it is a SH volume' );
assert( L>=2, 'This method makes no sense with trivial SH volumes with L=0' );
L = round(L);
%%% ----------------
assert(ismatrix(gi),'gi must be a matrix');
assert(ismatrix(bi),'bi must be a matrix');
assert(size(gi,2)==3,'gi must be Gx3');
assert(size(bi,2)==1,'bi must be Gx1');
assert(size(gi,1)==size(bi,1),'The first dimensions of gi and bi must be the same size');
% -------------------------------------------------------------------------
opt.mask = true(M,N,P); optchk.mask = [true,true];     % boolean the same size as the image field
opt.bth  = 100;         optchk.bth = [true,true];      % always 1x1 double
opt.chunksz = 100;      optchk.chunksz = [true,true];  % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];     % always 1x1 double
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
use_parallel           = use_parallel_test;
% -------------------------------------------------------------------------
[bs,ps,Ns] = auto_detect_shells(bi,opt.bth);
% Unroll and mask to work comfortably:
mask  = reshape( opt.mask, [M*N*P,1] ); % M*N*P  x 1
% ---
lpar  = reshape( lpar,     [M*N*P,1] ); % M*N*P  x 1
lpar  = lpar(mask,:);                   % Q x 1
lperp = reshape( lperp,    [M*N*P,1] ); % M*N*P  x 1
lperp = lperp(mask,:);                  % Q x 1
% ---
sh    = reshape( sh, [M*N*P,K] );       % (M*N*P) x K
sh    = sh(mask,:);                     % Q x K
% ---
Q     = size(sh,1);
% -------------------------------------------------------------------------
% Compute the attenuation signal at each q-space sample
atti0 = zeros(Q,size(gi,1));
% -------------------------------------------------------------------------
% Work chunk-by-chunk
ptr = dmri_sh_expand_coeffs(L); % 1 x K
Gn  = zeros(1,Ns);
Bn  = cell(1,Ns);
for n=1:Ns % For each shell...
    pn     = (ps==n); % ptr tells which samples belong to the current shell
    Bn{n}  = GenerateSHMatrix( L, gi(pn,:) )'; % K x Gn(n)
    Gn(n)  = size(Bn{n},2);
    Bn{n}  = reshape( Bn{n}, [1,K,Gn(n)] );    % 1 x K x Gn(n)
end
if(use_parallel)
    prp = gcp;
    CHUNKSZ = opt.chunksz*prp.NumWorkers*20;
    for CK=1:ceil(Q/CHUNKSZ)
        % -----------------------------------------------------------------
        IDI = (CK-1)*CHUNKSZ+1;
        IDF = min(CK*CHUNKSZ,Q);
        QQ  = IDF-IDI+1;
        % -----------------------------------------------------------------
        LPAR  = lpar(IDI:IDF);
        LPERP = lperp(IDI:IDF);
        SH    = sh(IDI:IDF,:);
        ATTI  = zeros(QQ,size(gi,1));
        % -----------------------------------------------------------------
        NCK    = ceil(QQ/opt.chunksz);
        LPARC  = cell(1,NCK);
        LPERPC = cell(1,NCK);
        SHC    = cell(1,NCK);
        ATTIC  = cell(1,NCK);
        % -----------------------------------------------------------------
        for ck=1:NCK
            idi = (ck-1)*opt.chunksz+1;
            idf = min(ck*opt.chunksz,QQ);
            LPARC{ck}  = LPAR(idi:idf);
            LPERPC{ck} = LPERP(idi:idf);
            SHC{ck}    = SH(idi:idf,:);
        end
        % -----------------------------------------------------------------
        parfor ck=1:NCK
            % --------------------------------------------
            % Compute the convolutuion weights:
            eODF = dmri_compute_convolution_weights_ODF( ...
                bs, LPARC{ck}, LPERPC{ck}, L ); % QC x (L/2+1) x Ns
            % --------------------------------------------
            % Actually compute the convolution:
            if(is_broadcast_available)
                sh2 = SHC{ck}.*eODF(:,ptr,:);  % QC x K x Ns
            else
                sh2 = bsxfun( @(x,y)(x.*y), ...
                    SHC{ck}, eODF(:,ptr,:) );    % QC x K x Ns
            end
            % --------------------------------------------
            for n=1:Ns % For each shell...
                pn = (ps==n); % ptr tells which samples belong to the current shell
                if(is_broadcast_available)
                    attin = sum( ...
                        sh2(:,:,n).*Bn{n}, ... % QC x K x Gn(n)
                        2 );                   %#ok<PFBNS> % QC x 1 x Gn(n)
                else
                    attin = sum( ...
                        bsxfun(@(x,y)(x.*y),sh2(:,:,n),Bn{n}), ... % QC x K x Gn(n)
                        2 );                                       % QC x 1 x Gn(n)
                end
                ATTIC{ck}(:,pn) = reshape(attin,[size(attin,1),Gn(n)]); %#ok<PFBNS>
            end
        end
        % -----------------------------------------------------------------
        for ck=1:NCK
            idi = (ck-1)*opt.chunksz+1;
            idf = min(ck*opt.chunksz,QQ);
            ATTI(idi:idf,:) = ATTIC{ck};
        end
        % -----------------------------------------------------------------
        atti0(IDI:IDF,:) = ATTI;
        % -----------------------------------------------------------------
    end
else
    for ck=1:ceil(Q/opt.chunksz)
        % --------------------------------------------
        idi = (ck-1)*opt.chunksz+1;
        idf = min(ck*opt.chunksz,Q);
        QC  = idf-idi+1;
        % --------------------------------------------
        % Compute the convolutuion weights:
        eODF = dmri_compute_convolution_weights_ODF( ...
            bs, lpar(idi:idf), lperp(idi:idf), L ); % QC x (L/2+1) x Ns
        % --------------------------------------------
        % Actually compute the convolution:
        if(is_broadcast_available)
            shc = sh(idi:idf,:).*eODF(:,ptr,:);  % QC x K x Ns
        else
            shc = bsxfun( @(x,y)(x.*y), ...
                sh(idi:idf,:), eODF(:,ptr,:) );    % QC x K x Ns
        end
        % --------------------------------------------
        for n=1:Ns % For each shell...
            pn = (ps==n); % ptr tells which samples belong to the current shell
            if(is_broadcast_available)
                attin = sum( ...
                    shc(:,:,n).*Bn{n}, ... % QC x K x Gn(n)
                    2 );                   % QC x 1 x Gn(n)
            else
                attin = sum( ...
                    bsxfun(@(x,y)(x.*y),shc(:,:,n),Bn{n}), ... % QC x K x Gn(n)
                    2 );                                       % QC x 1 x Gn(n)
            end
            atti0(idi:idf,pn) = reshape(attin,[QC,Gn(n)]);
        end
        % --------------------------------------------
    end
end
% -------------------------------------------------------------------------
% In case an actual volume with the free-water compartment was provided,
% use it here:
if(~isempty(f))
    f   = reshape( f, [M*N*P,1] );  % (M*N*P)x1
    f   = f(mask,:);                % Qx1
    E0  = (1-f)*exp(-opt.ADC0*bi'); % QxG
    if(is_broadcast_available)
        atti0 = atti0.*f;                           % QxG
    else
        atti0 = bsxfun( @(x,y)(x.*y), atti0, f );   % QxG
    end
    atti0 = atti0 + E0;
end
% -------------------------------------------------------------------------
% Create the output and assign:
atti = zeros(M*N*P,size(gi,1));
atti(mask,:) = atti0;
atti = reshape(atti,[M,N,P,size(gi,1)]);





























