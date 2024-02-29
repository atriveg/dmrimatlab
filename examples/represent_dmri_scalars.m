% represent_dmri_scalars.m
clear;

% This script illustrates a simple pipeline to read a DWI volume and
% represent ODFs at desired positions.

% Choose a valid 3D-Slicer nrrd file:
%filename = '/mnt/ARMARIO/dwiData/P6_b1300_r2_dwi.nrrd';
filename = '/Users/atriveg/Downloads/P6_b1300_r2_dwi.nrrd';

% Use 3rd party function from matlab bridge to load the volume in 3D-Slicer
% format:
vol = nrrdread(filename);

% Get metadata:
[origin,direction,space,measurement_frame] = slicervol2metadata(vol);

% Convert 3D-Slicer-like volume to a data set to be processed with the
% toolbox:
[dwi,giall,biall] = slicervol2dwi(vol); % giall and biall include baselines' gi=[0,0,0] and bi=0

% Convert to an attenuation signal:
[atti,gi,bi] = dwi2atti(dwi,giall,biall);

% Compute a mask to avoid unnecessary computations:
[~,mask] = dwi2otsuthreshold(dwi,biall);

% Create an ADC volume:
SH = atti2shadc( atti, gi, bi, 'mask', mask, 'L', 8, 'lambda', 1.0e-9 );

hf = figure(1001);
close(hf);
hf = figure(1001);
ha = axes('Parent',hf);
% ---
i1 = [11,size(SH,1)-10];
i2 = [21,size(SH,2)-20];
i3 = [round(size(SH,3)/2),round(size(SH,3)/2)];
% ---
% Try another use case: pass the background image directly. We pass the DC
% coefficient of the SH corresponding to the ADC, i.e., the mean
% diffusivity:
ha = plotdmri3d(SH,i1,i2,i3,'ha',ha, 'origin', origin, ...
    'direction', direction, 'space', space, ...
    'mframe', measurement_frame, 'bgimage', SH(:,:,:,1), 'bgsh', [], ...
    'bgalpha', 1, 'mask', mask, 'bbox', true, 'glyphs', false );
% ---
i1 = [round(size(SH,1)/2),round(size(SH,1)/2)];
i2 = [21,size(SH,2)-20];
i3 = [3,size(SH,3)-2];
% ---
ha = plotdmri3d(SH,i1,i2,i3,'ha',ha, 'origin', origin, ...
    'direction', direction, 'space', space, ...
    'mframe', measurement_frame, 'bgimage', 'fa', 'bgsh', SH, ...
    'bgalpha', 1, 'mask', mask, 'bbox', false, 'glyphs', false );
% ---
i1 = [11,size(SH,1)-10];
i2 = [round(size(SH,2)/2),round(size(SH,2)/2)];
i3 = [3,size(SH,3)-2];
% ---
ha = plotdmri3d(SH,i1,i2,i3,'ha',ha, 'origin', origin, ...
    'direction', direction, 'space', space, ...
    'mframe', measurement_frame, 'bgimage', 'color', 'bgsh', SH, ...
    'bgalpha', 1, 'mask', mask, 'bbox', false, 'glyphs', false );
% ---
axis('equal');
rotate3d('on');
xlabel('X','FontSize',16,'FontWeight','bold');
ylabel('Y','FontSize',16,'FontWeight','bold');
zlabel('Z','FontSize',16,'FontWeight','bold');


