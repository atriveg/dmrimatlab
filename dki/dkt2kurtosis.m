function [mk,ak,rk] = dkt2kurtosis( dkt, dti, varargin )
% function [MK,AK,RK] = dkt2kurtosis( dkt, dti, ...
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Computes the mean kurtosis mk, the axial kurtosis ak and the radial
%   kurtosis rk from the 4th order kurtosis tensor dkt and the 2nd order
%   diffusion tensor dti. According to the signal model:
%
%     log(atti) = log(S) - log(S_0)
%       = - b * sum_{i=1,j=1}^3 ni·nj·D_{ij}
%         + b^2/6*MD^2 * sum_{i=1,j=1,k=1,l=1}^3 ni·nj·nk·nl·W_{ijkl}
% 
%   where D_{ij} is the 2nd order diffusion tensor, W_{ijkl} is the 4th 
%   order kurtosis tensor and MD is the mean diffusivity computed from 
%   D_{ij}, the directional kurtosis for each unit direction n=[n1,n2,n3]^T 
%   reads:
%
%     K(n) = MD^2 * (sum_{i=1,j=1}^3 ni·nj·D_{ij})^{-2}
%                                 * sum_{i=1,j=1}^3 ni·nj·D_{ij}
%
%   If u1, u2 and u3 are the three orthogonal eigenvectors of D_{ij}, we
%   calculate:
%
%       MK = (K(u1)+K(u2)+K(u3))/3;
%       AK = K(u1);
%       RK = (K(u2)+K(u3))/2;
%
%   MANDATORY INPUTS:
%
%      dkt: MxNxPx15, where the last dimension contains the 15 unique
%         components of the 4th order Kurtosis tensor, W_{ijkl}. The
%         ordering of these components is consistent with that used in
%         HOT-related functions. It is returned by atti2dkt.
%      dti: MxNxPx6, where the last diemsnion contains the 6 unique 
%         components of the 2nd order diffusion tensor, D_{ik}. It is
%         returned by atti2dkt.
%
%   OUTPUTS:
%
%      MK, AK, RK: MxNxP each, the respective values of the mean kurtosis,
%      the axial kurtosis and the radial kurtosis.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      chunksz: the computation internally uses matrix products with SH
%         coefficients with chunks of data. This is the size of such chunks
%         (default: 1000).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      maxthreads: the algorithm interanlly resorts to dti2spectrum, which 
%         is run with multiple threads. This is the maximum allowed number 
%         of threads, which can indeed be reduced if it exceeds the number 
%         of logical cores (default: the number of logical cores in the 
%         machine).

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<2)
    error('At lest the dkt and the dti must be supplied');
end
% ---------
[M,N,P,K]        = size(dkt);
[M2,N2,P2,K2,L2] = size(dti);
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
% ---------
assert(isequal([M,N,P],[M2,N2,P2]),'The FoVs of dkt and dti do not match');
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true]; % always 1x1 char
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Convert the dkt volume to a volume of SH so that we can easily evaluate
% at random directions:
dkt = hot2sh( dkt, 'maxthreads', opt.maxthreads, 'mask', opt.mask );
% Compute the spectrum of the diffusion tensor:
[u1,u2,u3,l1,l2,l3] = dti2spectrum( dti, 'maxthreads', opt.maxthreads, 'mask', opt.mask );
% Unroll the dkt to work comfortably and easily apply the mask:
dkt = reshape( dkt, [M*N*P,K] );
dkt = dkt( opt.mask, : );
% Do the same for the u's volumes:
u1  = reshape( u1, [M*N*P,3] );
u1  = u1( opt.mask, : );
u2  = reshape( u2, [M*N*P,3] );
u2  = u2( opt.mask, : );
u3  = reshape( u3, [M*N*P,3] );
u3  = u3( opt.mask, : );
% Do the same for the l's volumes:
l1  = reshape( l1, [M*N*P,1] );
l1  = l1( opt.mask, : );
l2  = reshape( l2, [M*N*P,1] );
l2  = l2( opt.mask, : );
l3  = reshape( l3, [M*N*P,1] );
l3  = l3( opt.mask, : );
% Create empty arrays for the outputs and work chunk-by-chunk:
NV  = size(dkt,1);
ak1 = zeros( NV, 1 );
ak2 = zeros( NV, 1 );
ak3 = zeros( NV, 1 );
% Now, chunk-by-chunk, evaluate the kurtosis and diffusion for each of the
% main diffusion directions:
for ck=1:ceil(NV/opt.chunksz) % For each chunk with max size NV
    % The range for this chunk is:
    idi = (ck-1)*opt.chunksz+1;   % from
    idf = min(ck*opt.chunksz,NV); % to
    % Gather the principal diffusion directions for the voxels
    % of this chunk; from them, compute the SH encoding matrix:
    B1 = GenerateSHMatrix( 4, u1(idi:idf,:) ); % chunksz x 15
    B2 = GenerateSHMatrix( 4, u2(idi:idf,:) ); % chunksz x 15
    B3 = GenerateSHMatrix( 4, u3(idi:idf,:) ); % chunksz x 15
    % Gather all the SH coefficients for the voxels of this
    % chunk:
    x  = dkt(idi:idf,:); % chunksz x 15
    % Compute the kurtoses for each main direction:
    ak1(idi:idf,1) = sum(B1.*x,2); % chunksz x 1
    ak2(idi:idf,1) = sum(B2.*x,2); % chunksz x 1
    ak3(idi:idf,1) = sum(B3.*x,2); % chunksz x 1
    % Normalize to compute the actual kurtoses:
    MD = (l1(idi:idf,1)+l2(idi:idf,1)+l3(idi:idf,1))/3;
    ak1(idi:idf,1) = (MD.*MD).*ak1(idi:idf,1)./(l1(idi:idf,1).*l1(idi:idf,1));
    ak2(idi:idf,1) = (MD.*MD).*ak2(idi:idf,1)./(l2(idi:idf,1).*l2(idi:idf,1));
    ak3(idi:idf,1) = (MD.*MD).*ak3(idi:idf,1)./(l3(idi:idf,1).*l3(idi:idf,1));
end
% Cast the results to the proper size:
mk = zeros(M,N,P);
ak = zeros(M,N,P);
rk = zeros(M,N,P);
mk(opt.mask) = (ak1+ak2+ak3)/3;
ak(opt.mask) = ak1;
rk(opt.mask) = (ak2+ak3)/2;

end
