function [atti,gi,bi] = slicervol2atti(img)
% function atti = slicervol2atti(img)
%
%   Parses a DWI data set created with 3-D Slicer via MatlabBridge and the
%   SlicerVolumeToMatFile module (see 3rdparty folder) and outputs an 
%   attenuation volume (i.e, the dwi volume divide by the averaged 
%   baselines), gradients table, and b-values to be used with 
%   "atti2something" functions:
%
%      img: a structure produced from 3-D Slicer with the module named
%           SlicerVolumeToMatFile.m or nrrdread.m (see 3rdparty folder)
%
%      atti: a MxNxPxG double array with the Si/S0 attenuation signal at
%           each voxel within the MxNxP image frame and for each i=1...G
%           gradient direction.
%      gi:  a Gx3 matrix with the gradient directions
%      bi:  a Gx1 vector withe the b-values for each gradient direction

% -------------------------------------------------------------------------
% Avoid repeated calls to is_broadcast_available, which will always return
% the same value unless the toolbox is reconfigured:
is_broadcast_available = is_broadcast_available_test;
% -------------------------------------------------------------------------

atti = double(permute(img.pixelData,[2,3,4,1])); % MxNxPxG
[~,~,~,G] = size(atti);

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
gridx   = ~bsidx;              % gradients

baseline = mean( atti(:,:,:,bsidx), 4 ); % MxNxP
mask     = (baseline>1.0e-6);
baseline(~mask) = inf;
atti     = atti(:,:,:,gridx);            % MxNxPxGp
modules  = modules(gridx,:);  % Gpx1
gi       = gi(gridx,:);       % Gpx3
G        = size(atti,4);
bi       = bi*modules;
if(is_broadcast_available)
    atti = atti./baseline;
    gi  = gi./modules;
else
    for g=1:G
        atti(:,:,:,g) = atti(:,:,:,g)./baseline;
        gi(g,:) = gi(g,:)/modules(g);
    end
end



