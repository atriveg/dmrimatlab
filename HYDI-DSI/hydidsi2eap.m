function eap = hydidsi2eap( dsi, dti, lattice, Qx, Qy, Qz, ui, ri, varargin )
% function eap = hydidsi2eap( dsi, dti, lattice, Qx, Qy, Qz, ui, ri,
%                                     'opt1', value1, 'opt2', value2, ... )
%
%   Given the values of the EAP sampled in a regular Cartesian lattice,
%   as they are returned by atti2hydidsi, reconstruct the corresponding
%   ensemble average propagator at desired points by either nearest
%   neughbor, linear or cubic B-spline interpolation.
%
%   Mandatory inputs:
%
%      dsi: a MxNxPxK double array containing the sampled EAP at the
%         Cartesian grid covering (half the) R-space. This should be the
%         output of a call to atti2hydidsi().
%      dti: MxNxPx6, a double array with the estimated tensor model  at
%         each voxel. Note this is necessary to interpret the data in EAP,
%         since the EAP is described in a rotated space aligned with the
%         tensor eigenvectors. This should be the (second) output to a call
%         to atti2hydidsi().
%      lattice: 3x1, the lattice size necessary to interpret the Cartesian
%         grid where the EAP is defined. It is the third output to a call to
%         atti2hydidsi.
%      Qx, Qy, Qz: MxNxP each; they are the bandwidths estimated for the
%         EAP for each (rotated) direction, and are the 4th to 6th outputs
%         of atti2hydidsi.
%      ui: a Gx3 matrix with the directions table, each row corresponding
%         to a unit vector with the direction for which the EAP will be
%         evaluated.
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
%         number of entries in the directions table ui.
%
%   Optional arguments may be passed as name/value pairs in the regular
%   matlab style:
%
%   General optional parameters:
%
%      ADC0: estimated diffusivity of free water at body temperature. It is
%         used to determine the lower and upper bounds of the eigenvalues
%         of the dti and perform sanity checks (default: 3.0e-3).
%      interp: interpolation type, either 'n', 'l' or 's', to tell the
%         algorithm if either nearest neighbor, or linear or cubic B-spline
%         interpolation has to be used to obtain the values of the EAP
%         at continuous indices from the discrete sampling provided by
%         HYDI-DSI (default: 'l').
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
if(nargin<8)
    error('Too few input arguments');
end
[M,N,P,K] = size(dsi);
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
assert( isequal(size(dti),[M,N,P,6]), 'The DTI volume should be M x N x P x 6, with M,N,P the same as the EAP volume' );
assert( isequal([size(Qx,1),size(Qx,2),size(Qx,3)],[M,N,P]), 'The Qx volume must have size M x N x P, with M,N,P the same as the EAP volume' );
assert( isequal([size(Qy,1),size(Qy,2),size(Qy,3)],[M,N,P]), 'The Qy volume must have size M x N x P, with M,N,P the same as the EAP volume' );
assert( isequal([size(Qz,1),size(Qz,2),size(Qz,3)],[M,N,P]), 'The Qz volume must have size M x N x P, with M,N,P the same as the EAP volume' );
assert( isequal(size(lattice),[3,1]), 'Wrong lattice size. Must be 3 x 1' );
% Make sure the ui are normalized:
ni  = sqrt(sum(ui.*ui,2));
ini = ni;
ini(ni>0) = 1./ni(ni>0);
ui(:,1) = ui(:,1).*ini;
ui(:,2) = ui(:,2).*ini;
ui(:,3) = ui(:,3).*ini;
% Parse the optional input arguments:
% -------------------------------------------------------------------------
opt.ADC0 = 3.0e-3;      optchk.ADC0 = [true,true];       % always 1x1 double
opt.interp = 'l';       optchk.interp = [true,true];
opt.maxthreads = 1.0e6; optchk.maxthreads = [true,true];
opt.mask = true(M,N,P); optchk.mask = [true,true];       % boolean the same size as the image field
% -------------------------------------------------------------------------
opt = custom_parse_inputs(opt,optchk,varargin{:});
assert( any(opt.interp=='nls'), 'Interpolation type must be either ''n'', ''l'', or ''s''' );
% -------------------------------------------------------------------------
% Unroll and mask to work comfortably
mask = reshape(opt.mask,[M*N*P,1]);
dsi  = reshape(dsi,[M*N*P,K]);
dsi  = dsi(mask,:);
dti  = reshape(dti,[M*N*P,6]);
dti  = dti(mask,:);
Qx   = Qx(opt.mask);
Qy   = Qy(opt.mask);
Qz   = Qz(opt.mask);
eap  = zeros(M*N*P,G);
% -------------------------------------------------------------------------
% Call the mex function:
eapM = hydidsi2eap_( double(dsi'), double(dti'), double(lattice), ...
         double(Qx), double(Qy), double(Qz), double(ui), double(ri), ...
         opt.interp, double(opt.ADC0), double(opt.maxthreads) );
% -------------------------------------------------------------------------
% Reshape back to output the result:
eap(mask,:) = eapM';
eap = reshape(eap,[M,N,P,G]);
% -------------------------------------------------------------------------
end
