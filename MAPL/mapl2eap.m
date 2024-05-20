function eap = mapl2eap( mapl, dti, ui, ri, varargin )
% function eap = mapl2eap( mapl, dti, ui, ri,
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Given the coefficients of the MAP-MRI expansion, reconstruct the
%   corresponding ensemble average propagator, according to the model 
%   described by Ozarslan et al:
%    
%        Ozarslan E, Koay CG, Shepherd TM, Komlosh ME, Irfanoglu MO,
%        Pierpaoli C, Basser PJ. "Mean apparent propagator (MAP) MRI: a 
%        novel diffusion imaging method for mapping tissue microstructure". 
%        Neuroimage. 2013 Sep; 78:16-32.
%
%   This function is used to reconstruct the values of P(R) from a MAPL
%   model fitted with atti2mapl. 
%
%   Mandatory inputs:
%
%      mapl: a MxNxPxK double array containing the MAPL coefficients, where
%         K = (Nmax+2)(Nmax+4)(2*Nmax+3)/24 for some even integer Nmax>=0.
%         This should be the output of a call to atti2mapl().
%      dti: MxNxPx6, a double array with the estimated tensor model  at
%         each voxel. Note this is necessary to interpret the data in mapl,
%         since the coefficients refer to basis functions that are both
%         rotated and scaled depending on the properties of this dti. This
%         should be the (second) output to a call to atti2mapl().
%      ui: a Gx3 matrix with the directions table, each row corresponding
%         to a unit vector with the direction for which the EAP will be
%         evaluated
%      ri: a Gx1 vector with the distances to the origin corresponding to 
%         each of the G directions described above. Units are millimeters
%         (meaningful values will be obtained in the range of 0.01 mm).
%         Note that {ui,ri} describe a samplig of the R-space (EAP domain)
%         in spherical coordinates.
%
%   The outputs of the function are as follows:
%
%      eap: a MxNxPxG array with the values reconstructed for the EAP at
%         each of the MxNxP voxels within the FOV. Note G will match the
%         number of entries in the directions table gi.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   General optional parameters:
%
%      tau: the effective diffusion time of the acquisition sequence in
%         seconds (default: 70.0e-3, i.e. 70 milliseconds).
%      ADC0: estimated diffusivity of free water at body temperature. It is
%         used to determine the lower and upper bounds of the eigenvalues
%         of the dti and perform sanity checks (default: 3.0e-3).
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
if(nargin<4)
    error('At lest the MAPL volume, the DTI volume, the directions table, and the r-values must be supplied');
end
[M,N,P,K] = size(mapl);
if(~ismatrix(ui)||~ismatrix(ri))
    error('ui and ri must be 2-D matlab matrixes');
end
G = size(ui,1);
if(size(ui,2)~=3)
    error('The drections table ui must have size Gx3');
end
if(size(ri,1)~=G)
    error('The number of r-values ri must match the 1-st dimension of ui');
end
if(size(ri,2)~=1)
    error('The r-values vector must be a column vector');
end
assert( isequal(size(dti),[M,N,P,6]), 'The DTI volume should be M x N x P x 6, with M,N,P the same as the MAPL volume' );
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.tau = 70.0e-3;      optchk.tau = [true,true];        % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask = reshape(opt.mask,[M*N*P,1]);
mapl = reshape(mapl,[M*N*P,K]);
mapl = mapl(mask,:);
dti  = reshape(dti,[M*N*P,6]);
dti  = dti(mask,:);
eap  = zeros(M*N*P,G);
% -------------------------------------------------------------------------
% Call the mex function:
eapM = mapl2eap_( double(mapl'), double(dti'), ui, ri, opt.tau, opt.ADC0, opt.maxthreads );
% -------------------------------------------------------------------------
% Reshape back to output the result:
eap(mask,:) = eapM';
eap = reshape(eap,[M,N,P,G]);
% -------------------------------------------------------------------------
end
