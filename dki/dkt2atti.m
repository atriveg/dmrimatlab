function atti = dkt2atti( dkt, dti, S0, gi, bi, varargin )
% function atti = dkt2atti( dkt, dti, S0, gi, bi, 'opt1', value1, 'opt2', value2, ... )
%
%   Computes the attenuation signal from its diffusion tensor + kurtosis
%   tensor representation (and the correction S0 to the baseline image),
%   all of them returned by atti2dkt, according to the signal model:
%
%     log(atti) = log(S) - log(S_0)
%       = - b * sum_{i=1,j=1}^3 ni·nj·D_{ij}
%         + b^2/6*MD^2 * sum_{i=1,j=1,k=1,l=1}^3 ni·nj·nk·nl·W_{ijkl}
% 
%   where MD is the mean diffusivity at each voxel computed from D_{ij}.
%
%   MANDATORY INPUTS:
%
%      dkt: MxNxPx15, where the last dimension contains the 15 unique
%         components of the 4th order Kurtosis tensor, W_{ijkl}. The
%         ordering of these components is consistent with that used in
%         HOT-related functions. It is returned by atti2dkt.
%      dti: MxNxPx6, where the last diemsnion contains the 6 unique 
%         components of the 2nd order diffusion tensor, D_{ik}. It is
%         returned by either atti2dkt or atti2dti.
%      S0: MxNxP, the normalized baseline (a value near 1) to be applied to
%         the attenuation signal. It is returned by either atti2dkt or 
%         atti2dti.
%      gi: a Gx3 matrix with the gradient directions to be reconstructed, 
%         each row corresponding to a unit vector.
%      bi: a Gx1 vector with the b-values to be reconstructed, each one
%         corresponding to a row of gi.
%
%   OUTPUTS
%
%      atti: a MxNxPxG double array containing the attenuation signal 
%         sampled at the gradient vectors described by gi and bi.
% 
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros.
%      chunksz: the algorithm internally resorts to dti2signal, which
%         works by matrix multiplying data chunks. This parameter is 
%         directly passed to dti2signal (default: 1000).
%      maxthreads: the algorithm internally resorts to hot2signal, which is
%         multi-threaded. This is the maximum allowed number of threads, 
%         which can indeed be reduced if it exceeds the number of logical 
%         cores (default: the number of logical cores in the machine).

% -------------------------------------------------------------------------
% Check the mandatory input arguments:
if(nargin<5)
    error('At lest the dkt, the dti, the gradients table & b-values and S0 must be supplied');
end
% ---------
[M,N,P,K]        = size(dkt);
[M2,N2,P2,K2,L2] = size(dti);
if(isempty(S0))
    S0 = ones(M,N,P);
end
[M3,N3,P3]       = size(S0);
% ---------
G = size(gi,1);
assert( isequal(size(gi),[G,3] ), 'The b-vectors gi must have size Gx3, with G the 4-th dimension of atti' );
assert( isequal(size(bi),[G,1] ), 'The b-values bi must have size Gx1, with G the 4-th dimension of atti' );
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
assert(isequal([M,N,P],[M3,N3,P3]),'The FoVs of dkt and S0 do not match');
% -------------------------------------------------------------------------
% Parse the optional input arguments:
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt.chunksz = 1000;     optchk.chunksz = [true,true]; % always 1x1 double
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true]; % always 1x1 char
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% The part directly derived from the rank-2 diffusion tensor:
sig1 = dti2signal( dti, gi, 'mask', opt.mask, 'chunksz', opt.chunksz ); % M x N x P x G
sig1 = -sig1.*reshape(bi,[1,1,1,G]);
% The part directly derived from the kurtosis tensor:
sig2 = hot2signal( dkt, gi, 'mask', opt.mask, 'maxthreads',opt.maxthreads );
MD   = sum( dti(:,:,:,[1,4,6]), 4 )/3;
sig2 = sig2.*(MD.*MD).*reshape(bi.*bi,[1,1,1,G])/6;
% The compound logarithmic signal:
atti = log(S0) + sig1 + sig2;
% The final signal:
atti = exp(atti);

end
