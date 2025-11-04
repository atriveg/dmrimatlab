function [atti,gi,bi] = slicervol2atti(img,bth)
% function atti = slicervol2atti(img)
% function atti = slicervol2atti(img,bth)
%
%   Parses a DWI data set created with 3-D Slicer via MatlabBridge and the
%   SlicerVolumeToMatFile module (see 3rdparty folder) and outputs an 
%   attenuation volume (i.e, the dwi volume divide by the averaged 
%   baselines), gradients table, and b-values to be used with 
%   "atti2something" functions:
%
%      img: a structure produced from 3-D Slicer with the module named
%           SlicerVolumeToMatFile.m or nrrdread.m (see 3rdparty folder)
%      bth (optional): b-values below this are assumed to correspond to
%           baseline images
%
%      atti: a MxNxPxG double array with the Si/S0 attenuation signal at
%           each voxel within the MxNxP image frame and for each i=1...G
%           gradient direction.
%      gi:  a Gx3 matrix with the gradient directions
%      bi:  a Gx1 vector withe the b-values for each gradient direction

% -------------------------------------------------------------------------
if(nargin<2)
    bth = 0.01;
end
[dwi,gi,bi] = slicervol2dwi(img);
[atti,gi,bi] = dwi2atti(dwi, gi, bi, 'bth', bth );
