% represent_dmri_glyphs.m
clear;

% This script illustrates a simple pipeline to read a DWI volume and
% represent ODFs at desired positions.

% Choose a valid 3D-Slicer nrrd file:
filename = '/mnt/ARMARIO/dwiData/csiro.nrrd';
%filename = '/Users/atriveg/Downloads/P6_b1300_r2_dwi.nrrd';
%filename = '/Users/atriveg/Downloads/australia.nrrd';

% Use 3rd party function from matlab bridge to load the volume in 3D-Slicer
% format:
vol = nrrdread(filename);

% Get metadata:
[origin,direction,space,measurement_frame] = slicervol2metadata(vol);

% Convert 3D-Slicer-like volume to a data set to be processed with the
% toolbox:
[dwi,giall,biall] = slicervol2dwi(vol); % giall and biall include baselines' gi=[0,0,0] and bi=0

% Estimate the noise (sigma) together with a background mask computed with
% Otsu's method:
[sigma,mask] = dwi2noisestd(dwi,biall,'estrad',[1;1;1]);
fprintf(1,'Estimated sigma: %1.2f\n',sigma);

% Clean the volumes using the denoising filters. Use the fastest approach:
fprintf(1,'Cleaning the DWI volumes to improve the PSNR. May take a while... ');
if(exist('csiro_filter.mat','file'))
    load('csiro_filter.mat');
else
    dwi = dwi2cleandwi( dwi, giall, sigma, ...
        'beta', 100.0, 'rs', [2;2;2], 'mask', mask, 'onlyUNLM', false, 'Ng', 15 );
    save('csiro_filter.mat','dwi');
end
fprintf(1,'done\n');

% Convert to an attenuation signal:
[atti,gi,bi] = dwi2atti(dwi,giall,biall);

% Compute the SH coefficients of the ADC of the volume to compute the FA as
% a background image:
SHadc = atti2shadc( atti, gi, bi, 'L', 2, 'lambda', 1.0e-3, 'mask', mask );

% Compute the OPDT to represent the glyphs:
SHodf = atti2shodf( atti, gi, 'L', 6, 'lambda', 0.006, 'type', 'copdt', 'mask', mask );

hf = figure(1001);
close(hf);
hf = figure(1001);
ha = axes('Parent',hf);
% ---
i1 = [1,size(SHodf,1)];
i2 = [1,size(SHodf,2)];
i3 = [32,32];
% ---
% Plot an axial slice without glyphs, just the FA
[ha,ho] = plotdmri3d(SHodf,i1,i2,i3,'ha',ha, 'origin', origin, ...
    'direction', direction, 'space', space, ...
    'mframe', measurement_frame, 'bgimage', 'fa', 'bgsh', SHadc, ...
    'bgalpha', 1, 'mask', mask, 'bbox', true, 'glyphs', false );
% Remove the bounding boxes, keep the text:
delete(ho.bbox(1));
delete(ho.bbox(2));
% ---
i1 = [44,60];
i2 = [44,67];
% ---
% Plot a piecce of the axial slice with the glyphs, weighted by the FA, but
% setting the transparency of the FA to 0 to avoid duplicating
[ha,ho] = plotdmri3d(SHodf,i1,i2,i3,'ha',ha, 'origin', origin, ...
    'direction', direction, 'space', space, ...
    'mframe', measurement_frame, 'bgimage', 'fa', 'bgsh', SHadc, ...
    'bgalpha', 0, 'mask', mask, 'bbox', false, 'glyphs', true, ...
    'angleres', 642, 'glyphspc', [1;1;1], 'glyphsc', 2.5, ...
    'glyphscvar', 'bg' );
light;
axis('equal');
rotate3d('on');
xlabel('X','FontSize',16,'FontWeight','bold');
ylabel('Y','FontSize',16,'FontWeight','bold');
zlabel('Z','FontSize',16,'FontWeight','bold');

title('Press ENTER to zoom in');
pause;

% Set an axial view and make a manual zoom:
CTAR = [ (size(SHodf,1)+1)/2-1; 15/4+(size(SHodf,2)+1)/2-1; (size(SHodf,3)+1)/2-1 ];
CTAR = ([direction,origin]*[CTAR;1])';
P1   = [0;0;0];
P2   = [size(SHodf,1)-1;size(SHodf,2)-1;0];
P1   = ([direction,origin]*[P1;1])';
P2   = ([direction,origin]*[P2;1])';
W    = norm(P1-P2)/sqrt(2);
CPOS = CTAR + W*[0,0,1];
set(ha,'CameraPosition',CPOS);
set(ha,'CameraTarget',CTAR);
set(ha,'CameraUpVector',[0,1,0]);
set(ha,'CameraViewAngle',11);
