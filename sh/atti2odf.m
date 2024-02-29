function [odf,SHodf] = atti2odf( atti, gi, gi2, varargin )
% function [odf,SHodf] = atti2odf( atti, gi, gi2, 'opt1', val1, 'opt2', val2, ... )
%
%   Takes an attenuation signal atti (DWI signal over the baseline) sampled
%   at a set of gradients gi and computes corresponding Orientation
%   Distribution Functions (ODF), either probabilistic (Jacobian weighted)
%   or non probabilistic (non-weighted) in the basis of spherical harmonics
%   by using (Tikhonov regularized) least squares over a mono-exponential 
%   model. The SH series expansion is evaluatedat a set of desired
%   directions gi2 to reconstruct the ODF.
%   This function basically implementes the ODF estimators described in the
%   following papers:
%
%     - Tristan-Vega, A., C-F. Westin, and S. Aja-Fernandez, "A new 
%     methodology for the estimation of fiber populations in the white 
%     matter of the brain with the Funk?Radon transform", NeuroImage, 
%     vol. 49, no. 2: Elsevier, pp. 1301?1315, 2010.
%
%     - Tristan-Vega, A., S. Aja-Fernandez, and C-F. Westin, "On the 
%     Blurring of the Funk-Radon Transform in Q?Ball Imaging", Medical 
%     Image Computing and Computer-Assisted Intervention?MICCAI 2009: 
%     Springer Berlin Heidelberg, pp. 415?422, 2009.
%
%     - Tristan-Vega, A., C-F. Westin, and S. Aja-Fernandez, "Estimation of
%     fiber orientation probability density functions in high angular 
%     resolution diffusion imaging", NeuroImage, vol. 47, no. 2: Elsevier, 
%     pp. 638?650, 2009.
%
%   which we ask you to cite in case you use this software for your
%   research.
%
%      atti: a MxNxPxG double array containing the S_i/S_0 (diffusion
%         gradient over non-weighted baseline) at each voxel within the
%         MxNxP image frame and for each of the G acquired image gradients.
%      gi: a Gx3 matrix with the gradients table, each row corresponding to
%         a unit vector with the direction of the acquired gradient.
%      gi2: a G2x3 matrix with the directions the ODF must be evaluated at,
%         each row corresponding to a unit vector
%
%      odf: a MxNxPxG2 array with the ODF computed at each voxel and
%         sampled at desired directions gi2
%      SHodf: a MxNxPxK array, with K=(L+1)(L+2)/2, with the coefficients
%         of the spherical harmonics representing the signal at each voxel
%         m,n,p
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%      type: a string that can be either: 'opdt', 'copdt', 'popdt'
%         (probabilistic, jacobian-weighted estimators), 'qballs',
%         'cqballs', or 'pqballs' (non probabilistic estimators) (default:
%         'opdt').
%      q0: this parmeter is used only for 'cqballs' and 'pqballs'. This is
%         the module of the variable q for the measured shell, the dual of 
%         R in the Fourier domain, defined as q0 = gamma*delta*G/(2*pi),
%         with G the module of the sensitizing gradients, gamma the 
%         gyromagnetic ratio, and delta the gradient duration. Then 
%         b = 4*pi^2*tau*q0^2, with tau the effective diffusion time
%         (default: 35 mm^(-1)).
%      L: an even integer with the maximum order of the SH to be used
%         (default: 6).
%      lambda: the Laplace-Beltrami regularization parameter for the linear
%         least squares problem (default 0.006).
%      chunksz: the LLS problem reduces to the product of the dwi signal
%         by an inverse matrix that may be pre-computed for the whole data
%         set. To improve the performance, cunksz voxels are gathered
%         together in a single matrix that is pre-multiplied by the LLS
%         inverse at each step, hence taking advantage of matlab's
%         capabilities (default: 100).
%      tl, tu: the lower and upper thresholds, respectively, defining the
%         range the dwi will lay within, so that tl should be close to 0
%         and tu should be close to 1 (default: 1.0e-5, 1-1.0e-5).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         true(M,N,P))
%      O1: due to the cropping of the infinite SH expansions, represented
%         signals have a certain error that increases if the fiber bundle
%         is far from the x-y plane. If O1>0, N directions of the raw ODF
%         (computed as dirs=icosamplesSphere(O1)) are probed at each voxel
%         to find the maximum, which is assumed to correspond to the main
%         fiber bundle at such voxel. Then, the ODF is recomputed in a new
%         basis where the x-axis is aligned with this main direction and
%         further rotated to match to the original basis (default: 0)

% Check the mandatory input argments:
if(nargin<2)
    error('At least the atti volume and the gradients table must be supplied');
end
[M,N,P,G] = size(atti);
NV = M*N*P; % Total number of voxels to be processed
if( ~ismatrix(gi) || ~ismatrix(gi2) )
    error('gi and g2must be 2-d matlab matrixes');
end
if(size(gi,1)~=G)
    error('The number of rows in gi must match the 4-th dimension of dwi');
end
if(size(gi,2)~=3)
    error('The gradients table gi must have size Gx3');
end
if(size(gi2,2)~=3)
    error('The directions gi2 must have size G2x3');
end

% Parse the optional input arguments:
opt.type = 'opdt';      optchk.type = [true,false];   % string with variable size
opt.q0 = 35;            optchk.q0 = [true,true];      % always 1x1 double
opt.L = 6;              optchk.L = [true,true];       % always 1x1 double
opt.lambda = 0.006;     optchk.lambda = [true,true];  % always 1x1 double
opt.chunksz = 100;      optchk.chunksz = [true,true]; % always 1x1 double
opt.tl = 1.0e-5;        optchk.tl = [true,true];      % always 1x1 double
opt.tu = 1-opt.tl;      optchk.tu = [true,true];      % always 1x1 double
opt.mask = true(M,N,P); optchk.mask = [true,true];    % boolean with the size of the image field
opt.O1 = 0;             optchk.O1 = [true,true];      % always 1x1 double
opt = custom_parse_inputs(opt,optchk,varargin{:});

if(opt.O1<1)
    % This is the simple case, where we need just computing the SH
    % coefficients and develop the signal:
    SHodf = atti2shodf( atti, gi, 'type', opt.type, 'q0', opt.q0, ...
        'L', opt.L, 'lambda', opt.lambda, 'chunksz', opt.chunksz, ...
        'tl', opt.tl, 'tu', opt.tu, 'mask', opt.mask );
    odf = sh2signal( SHodf, gi2, 'chunksz', opt.chunksz, ...
        'mask', opt.mask );
else
    % We need to find a proper rotation of the SH bases so that the main
    % fiber bundle becomes aligned with the x-axis
    SHodf = atti2shodf( atti, gi, 'type', opt.type, 'q0', opt.q0, ...
        'L', 2, 'lambda', opt.lambda, 'chunksz', opt.chunksz, ...
        'tl', opt.tl, 'tu', opt.tu, 'mask', opt.mask );
    dirs = icosamplesSphere(opt.O1,'O1',true,'verbose',false);
    odf = sh2signal( SHodf, dirs, 'chunksz', opt.chunksz, ...
        'mask', opt.mask );
    [~,arg] = max(odf,[],4);
    arg(~opt.mask) = 0;
    % Now, arg is MxNxP; each entry is an integer pointing to the direction
    % dirs(arg(m,n,p),:) align with the main fiber bundle at voxel m,n,p.
    atti  = reshape(atti,[NV,G]);
    arg   = reshape(arg,[NV,1]);
    odf   = zeros(NV,size(gi2,1));
    SHodf = zeros(NV,(opt.L+1)*(opt.L+2)/2);
    for j=1:size(dirs,1) % For each main direction:
        % Which voxels correspond to this main direction?
        ptr = (arg==j);
        if(~any(ptr)) % None!
            continue;
        end
        % This set of voxels are aligned with direction dirs(j,:). Design a
        % rotation that converts dirs(j,:) to [1;0;0]
        ROT = optimal_rotation_odf(dirs(j,:)');
        % Find the SH expansion under the rotated basis:
        atti_j = atti(ptr,:);
        N_j    = size(atti_j,1);
        atti_j = reshape(atti_j,[N_j,1,1,G]);
        gi_j   = (ROT*(gi'))';
        SH_j   = atti2shodf( atti_j, gi_j, ...
            'type', opt.type, 'q0', opt.q0, ...
            'L', opt.L, 'lambda', opt.lambda, 'chunksz', opt.chunksz, ...
            'tl', opt.tl, 'tu', opt.tu );
        % Develop the signal in the appropriate basis
        gi2_j  = (ROT*(gi2'))';
        odf_j  = sh2signal( SH_j, gi2_j, 'chunksz', opt.chunksz );
        % Cast the result in place:
        odf(ptr,:)   = reshape(odf_j,[N_j,size(gi2,1)]);
        SHodf(ptr,:) = reshape(SH_j, [N_j,size(SHodf,2)] );
    end
end
odf   = reshape(odf,[M,N,P,size(gi2,1)]);
SHodf = reshape(SHodf,[M,N,P,(opt.L+1)*(opt.L+2)/2]);

function ROT = optimal_rotation_odf(u)
[~,pos] = max(abs(u));
v = zeros(3,1);
v(pos) = 1;
if((v'*u)<0)
    u = -u;
end
ROT = eye(3);
w   = cross(u,v);
if(norm(w)>1.0e-3)
    st  = norm(u);
    ct  = v'*u;
    th  = atan2(st,ct);
    w   = w/sin(th);
    w   = [0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
    ROT = ROT + sin(th)*w + (1-cos(th))*(w*w);
end
