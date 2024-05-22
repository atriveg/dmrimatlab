function odf = mapl2odf( mapl, dti, ui, varargin )
% function odf = mapl2odf( mapl, dti, ui,
%                                 'opt1', value1, 'opt2', value2, ... )
%
%   Given the coefficients of the MAP-MRI expansion, reconstruct the
%   corresponding probabilistic Orientation Distribution Function (ODF),
%   according to the model described by Ozarslan et al:
%    
%        Ozarslan E, Koay CG, Shepherd TM, Komlosh ME, Irfanoglu MO,
%        Pierpaoli C, Basser PJ. "Mean apparent propagator (MAP) MRI: a 
%        novel diffusion imaging method for mapping tissue microstructure". 
%        Neuroimage. 2013 Sep; 78:16-32.
%
%   This function is used to reconstruct the values of Phi(u) from a MAPL
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
%         to a unit vector with the direction for which the ODF will be
%         evaluated.
%
%   The outputs of the function are as follows:
%
%      odf: a MxNxPxG array with the values reconstructed for the ODF at
%         each of the MxNxP voxels within the FOV. Note G will match the
%         number of entries in the directions table ui.
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
%      contrast: the exponent of r in the radial integral that defines the
%         ODF. With contrast=2, one gets the actual jacobian of the
%         spherical coordinates, so that the ODF is the actual probability
%         density function in the orientation space, as described by
%         Tristan-Vega et al., NeuroImage (2009). With contrast=0, one gets
%         the ODF described by Tuch for Q-Balls. For any other real value
%         strictly greater than -1, one gets different contrasts (default:
%         2.0, i.e. actual probabilistic ODF).
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
if(nargin<3)
    error('At lest the MAPL volume, the DTI volume, and the directions table must be supplied');
end
[M,N,P,K] = size(mapl);
if( ~ismatrix(ui) )
    error('ui must be 2-D matlab matrixes');
end
G = size(ui,1);
if(size(ui,2)~=3)
    error('The drections table ui must have size Gx3');
end
assert( isequal(size(dti),[M,N,P,6]), 'The DTI volume should be M x N x P x 6, with M,N,P the same as the MAPL volume' );
% Make sure the ui are normalized:
ni = sqrt(sum(ui.*ui,2));
ui(:,1) = ui(:,1)./ni;
ui(:,2) = ui(:,2)./ni;
ui(:,3) = ui(:,3)./ni;
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.tau = 70.0e-3;      optchk.tau = [true,true];        % always 1x1 double
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.contrast = 2.0;     optchk.contrast = [true,true];   % always 1x1 double
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
% -------------------------------------------------------------------------
if(opt.contrast<=-1)
    error('With the selected contrast value, the radial integral is non-convergent');
end
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask = reshape(opt.mask,[M*N*P,1]);
mapl = reshape(mapl,[M*N*P,K]);
mapl = mapl(mask,:);
dti  = reshape(dti,[M*N*P,6]);
dti  = dti(mask,:);
odf  = zeros(M*N*P,G);
% -------------------------------------------------------------------------
% Call the mex function:
odfM = mapl2odf_( double(mapl'), double(dti'), ui, opt.tau, opt.ADC0, opt.contrast, opt.maxthreads );
% -------------------------------------------------------------------------
% Reshape back to output the result:
odf(mask,:) = odfM';
odf = reshape(odf,[M,N,P,G]);
% -------------------------------------------------------------------------
end
