function Ku = dkt2dirkurtosis( dkt, dti, ui, varargin )
% function Kn = dkt2dirkurtosis( dkt, dti, ni, ...
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Computes the value of the kurtosis at desired directions ui from the
%   4th order kurtosis tensor dkt (W_ijkl) and the 2nd order diffusion 
%   tensor dti (D_ij). According to the signal model:
%
%     log(atti) = log(S) - log(S_0)
%       = - b * sum_{i=1,j=1}^3 ni·nj·D_{ij}
%         + b^2/6*MD^2 * sum_{i=1,j=1,k=1,l=1}^3 ni·nj·nk·nl·W_{ijkl}
% 
%   the kurtosis computed at an arbitrary unit direction
%   n = [n_1,n_2,n_3]^T reads:
%
%     K(n) = MD^2 * (sum_{i=1,j=1}^3 ni·nj·D_{ij})^{-2}
%                            * sum_{i=1,j=1,k=1,l=1}^3 ni·nj·nk·nl·W_{ijkl}
%
%   which is the value computed by this function.
%
%   MANDATORY INPUTS:
%
%      dkt: MxNxPx15, where the last dimension contains the 15 unique
%         components of the 4th order kurtosis tensor, W_{ijkl}. The
%         ordering of these components is consistent with that used in
%         HOT-related functions. It is returned by atti2dkt.
%      dti: MxNxPx6, where the last dimension contains the 6 unique 
%         components of the 2nd order diffusion tensor, D_{ik}. It is
%         returned by atti2dkt (or atti2dti).
%      ni: MxNxPx3, where the last dimension represents a unit vector at
%         each voxel within the FoV. It might be, for example, one of the
%         eigenvectors of dti as they are computed with dti2spectrum. If ni
%         is, in particular, the main eigenvector, then the function will
%         return the axial kurtosis at each voxel.
%
%   OUTPUTS:
%
%      Kn: MxNxP, the directional kurtosis computed at each voxel for the
%         given direction ni(x,y,z,:).
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      chunksz: the computation internally uses matrix products with SH
%         coefficients with chunks of data. This is the size of such chunks
%         (default: 1000).
%      maxthreads: the algorithm internally resorts to hot2sh, which 
%         is run with multiple threads. This is the maximum allowed number 
%         of threads, which can indeed be reduced if it exceeds the number 
%         of logical cores (default: the number of logical cores in the 
%         machine).

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<3)
    error('At lest the dkt, the dti and ni must be supplied');
end
% ---------
[M,N,P,K]        = size(dkt);
[M2,N2,P2,K2,L2] = size(dti);
[M3,N3,P3,K3]    = size(ui);
% ---------
assert(K==15,'The 4-th dimension of a proper dkt volume must have size 15');
if(K2==3)
    if(L2==3)
        dti = cat( 4, dti(:,:,:,1,1), dti(:,:,:,1,2), dti(:,:,:,1,3), ...
            dti(:,:,:,2,2), dti(:,:,:,2,3), dti(:,:,:,3,3) );
    else
        error('The dti volume provided does not look like a diffusion tensor volume');
    end
elseif(K2==6)
    if(L2~=1)
        error('The dti volume provided does not look like a diffusion tensor volume');
    end
else
    error('The dti volume provided does not look like a diffusion tensor volume');
end
assert(K3==3,'The ni volume must have 4th dimension with size 3');
% ---------
assert(isequal([M,N,P],[M2,N2,P2]),'The FoVs of dkt and dti do not match');
assert(isequal([M,N,P],[M3,N3,P3]),'The FoVs of dkt and ni do not match');
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true]; % always 1x1 char
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Convert the dkt volume to a volume of SH so that we can easily evaluate
% at arbitrary directions:
dkt = hot2sh( dkt, 'maxthreads', opt.maxthreads, 'mask', opt.mask );
% Unroll the dkt, the dti and ui to work comfortably and easily apply the mask:
dkt = reshape( dkt, [M*N*P,K] );
dkt = dkt( opt.mask, : );
dti = reshape( dti, [M*N*P,6] );
dti = dti( opt.mask, : );
ui = reshape( ui, [M*N*P,3] );
ui = ui( opt.mask, : );
% Create an empty array for the output and work chunk-by-chunk:
NV  = size(dkt,1);
Ku_ = zeros( NV, 1 );
% Now, chunk-by-chunk, evaluate the kurtosis and diffusion for each of the
% main diffusion directions:
for ck=1:ceil(NV/opt.chunksz) % For each chunk with max size NV
    % The range for this chunk is:
    idi = (ck-1)*opt.chunksz+1;   % from
    idf = min(ck*opt.chunksz,NV); % to
    % Gather the principal diffusion directions for the voxels
    % of this chunk; from them, compute the SH encoding matrix:
    B1 = GenerateSHMatrix( 4, ui(idi:idf,:) ); % chunksz x 15
    % Compute the encoding matrix for the diffusion tensor as well:
    B2 = [ ui(idi:idf,1).*ui(idi:idf,1), 2*ui(idi:idf,1).*ui(idi:idf,2), ...
        2*ui(idi:idf,1).*ui(idi:idf,3), ui(idi:idf,2).*ui(idi:idf,2), ...
        2*ui(idi:idf,2).*ui(idi:idf,3), ui(idi:idf,3).*ui(idi:idf,3) ];
    % Gather all the SH coefficients for the voxels of this chunk:
    x1  = dkt(idi:idf,:); % chunksz x 15
    % Gather all the DT coefficients for the voxels of this chunk:
    x2  = dti(idi:idf,:); % chunksz x 6
    % Compute the value of both tensors at the desired directions:
    kurtosis  = sum(B1.*x1,2); % chunksz x 1
    diffusion = sum(B2.*x2,2); % chunksz x 1
    % Compute the MD at each voxel:
    meandiff  = ( x2(:,1) + x2(:,4) + x2(:,6) )/3;
    % The directional kurtosis is:
    Ku_(idi:idf) = kurtosis.*(meandiff.*meandiff)./(diffusion.*diffusion);
end
% Cast the results to the proper size:
Ku = zeros(M,N,P);
Ku(opt.mask) = Ku_;
end
