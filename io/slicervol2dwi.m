function [dwi,gi,bi] = slicervol2dwi(img)
% function dwi = slicervol2dwi(img)
%
%   Parses a DWI data set created with 3-D Slicer via MatlabBridge and the
%   SlicerVolumeToMatFile module (see 3rdparty folder) and outputs a dwi 
%   signal (i.e., the diffusion weighted images without normalization and 
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

modules = sqrt(sum(gi.*gi,2)); % Gx1
bsidx   = (modules<0.01);      % baselines

gi = gi./repmat(modules,[1,3]);
bi = bi.*modules;
gi(bsidx,:) = 0;
bi(bsidx,:) = 0;

