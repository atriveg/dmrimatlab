function index = mapl2index( mapl, dti, varargin )
% function index = mapl2index( mapl, dti,
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Given the coefficients of the MAP-MRI expansion, compute some related
%   scalar indices, according to the model described by  Ozarslan et al:
%    
%        Ozarslan E, Koay CG, Shepherd TM, Komlosh ME, Irfanoglu MO,
%        Pierpaoli C, Basser PJ. "Mean apparent propagator (MAP) MRI: a 
%        novel diffusion imaging method for mapping tissue microstructure". 
%        Neuroimage. 2013 Sep; 78:16-32.
%
%   Mandatory inputs:
%
%      mapl: a MxNxPxK double array containing the MAPL coefficients, where
%         K = (Nmax+2)(Nmax+4)(2*Nmax+3)/24 for some even integer Nmax>=0.
%         This should be the output to a call to atti2mapl().
%      dti: MxNxPx6, a double array with the estimated tensor model  at
%         each voxel. Note this is necessary to interpret the data in mapl,
%         since the coefficients refer to basis functions that are both
%         rotated and scaled dependint on the properties of this dti. This
%         should be the (second) output to a call to atti2mapl().
%
%   The outputs of the function are as follows:
%
%      index: a MxNxP array with the values of the desired scalar index at
%         each of the MxNxP voxels within the FOV.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   General optional parameters:
%
%      ADC0: estimated diffusivity of free water at body temperature. It is
%         used to determine the lower and upper bounds of the eigenvalues
%         of the dti and perform sanity checks (default: 3.0e-3).
%      tau: the effective diffusion time of the acquisition, in seconds, 
%         necessary to compute the indices (default: 35e-3).
%      type: the particular index to be computed, a single char, one of:
%                'o' for RTOP
%                'a' for RTAP
%                'p' for RTPP
%                'd' for PA_DTI
%                'g' for NG
%                'e' for E(0)
%                'm' for MSD
%                'q' for QIV
%                'l' for the energy of the Laplacian of E(q)
%                'u' for the u0 describing the "most similar isotropic
%                    propagator" as in eq. (48) of Ozarslan's paper
%         (default: 'o');
%
%   Other general options:
%
%      maxthreads: the algorithm is run as a multi-threaded mex. This is 
%         the maximum allowed number of threads, which can indeed be 
%         reduced if it exceeds the number of logical cores (default: the 
%         number of logical cores in the machine).
%      mask: a MxNxP array of logicals. Only those voxels where mask is
%         true are processed, the others are filled with zeros (default:
%         all trues).

% Check the mandatory input argments:
if(nargin<2)
    error('At lest the MAPL volume and the DTI volume must be supplied');
end
[M,N,P,K] = size(mapl);
assert( isequal(size(dti),[M,N,P,6]), 'The DTI volume should be M x N x P x 6, with M,N,P the same as the MAPL volume' );
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.tau = 35.0e-3;      optchk.tau = [true,true];        % always 1x1 double
opt.type = 'o';         optchk.type = [true,true];       % always 1x1 char
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Sanity checks on the type of index:
if(isempty(strfind('oapdgemqlu',opt.type)))
    error('Unknown index type: %s,',opt.type);
end
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask  = reshape(opt.mask,[M*N*P,1]);
mapl  = reshape(mapl,[M*N*P,K]);
mapl  = mapl(mask,:);
dti   = reshape(dti,[M*N*P,6]);
dti   = dti(mask,:);
index = zeros(M*N*P,1);
% -------------------------------------------------------------------------
% Call the mex function:
indexM = mapl2index_( double(mapl'), double(dti'), opt.ADC0, opt.tau, opt.type, opt.maxthreads );
% -------------------------------------------------------------------------
% Reshape back to output the result:
index(mask,:) = indexM';
index = reshape(index,[M,N,P]);
% -------------------------------------------------------------------------
end
