function [dwi,gi,bi] = slicervol2dwi(img)
% function [dwi,gi,bi] = slicervol2dwi(img)
%
%   Parses a DWI data set read from a NRDD file with nrrdread and outputs a
%   DWI (i.e., the diffusion weighted images without normalization and
%   the baselines), gradients table (including those with zero length 
%   corresponding to the baselines), and b-values (including those zeros 
%   corresponding to the baselines) to be used with "dwi2something" 
%   functions.
%
%      img: a structure produced from 3-D Slicer with the module named
%           SlicerVolumeToMatFile.m or nrrdread.m (see 3rdparty folder)
%
%      dwi: a MxNxPx(G+B) double array with the diffusion signals and the
%           baslines at each voxel within the MxNxP image frame and for 
%           each i=1...G gradient direction and each j=1...B baseline.
%      gi:  a (G+B)x3 matrix with the gradient directions (zero-vectors
%           correspond to baselines).
%      bi:  a (G+B)x1 vector withe the b-values for each gradient direction
%           (zeros correspond to baselines).

dwi = double(permute(img.pixelData,[2,3,4,1])); % MxNxPx(G+B)
[~,~,~,G] = size(dwi);

gi = zeros(G,3);
fields = fieldnames(img.metaData);
for f=1:numel(fields)
    field = fields{f};
    if(length(field)>15)
        if(strcmp(field(1:15),'DWMRI_gradient_'))
            gidx = str2double(field(16:end));
            grad = img.metaData.(field);
            grad = sscanf(grad,'%f',Inf);
            gi(gidx+1,:) = grad(:)';
        end
    end
    if(length(field)==13)
        if(strcmp(field,'DWMRI_b_value'))
            bi = str2double(img.metaData.(field));
        end
    end
end

% Normalize according to:
%
% https://www.na-mic.org/wiki/NAMIC_Wiki:DTI:Nrrd_format#Describing_DWIs_with_different_b-values
%
% NOTE: Accroding to vtkMRMLNRRDStorageNode::ParseDiffusionInformation(),
%       each bi is computed as:
%
%          bi = b0*( ||gi|| / max_i{||gi||} )^2,
%
%       where b0 is the reference value passed with the DWMRI_b-value tag
modules = sqrt(sum(gi.*gi,2)); % Gx1
maxmod  = max(modules);
gi      = gi./modules;         % Gx3, broadcast
bi      = bi*(modules./maxmod).*(modules./maxmod);
bsidx   = (modules<0.01);      % baselines
gi(bsidx,:) = 0;
bi(bsidx,:) = 0;
