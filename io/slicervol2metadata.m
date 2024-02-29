function [origin,direction,space,measurement_frame] = slicervol2metadata(img)
% function [origin,direction,space,measurement_frame] = slicervol2metadata(img)
%
%   Parses a DWI data set created with 3-D Slicer via MatlabBridge and the
%   SlicerVolumeToMatFile module (see 3rdparty folder) and outputs the meta
%   data associated with it:
%
%      img: a structure produced from 3-D Slicer with the module named
%           SlicerVolumeToMatFile.m or nrrdread.m (see 3rdparty folder)
%
%      origin: 3x1 vector of doubles
%      direction: 3x3 matrix of doubles
%            These two parameters together relate the index space within
%            the (3+1)-D array Sh and the real world, physical coordinates
%            of the imaged volume; if m\in[i1(1),i1(2)], n\in[i2(1),i2(2)],
%            and p\in[i3(1),i3(2)], the the voxel (m,n,p) corresponds to
%            the physical location:
%               [x;y;z] = [ direction, origin ]*[m;n;p;1]
%      space: a string of the form: 
%                + RAS, RAI, RPS, RPI, LAS, LAI, LPS, or LPI; or:
%                + right-anterior-superior, right-anterior-inferior, ...,
%                  left-posterior-inferior
%            describing the anatomical convention. R means the increasing
%            direction of the x axis goes from left to right; A means the
%            increasing direction of the y axis goes from the front to the
%            back; S means the increasing direction o the z axis goes
%            upwards
%      measurement_frame: the measurment frame, a 3x3 rotation matrix that
%            relates the diffusion measurements to the RAS space (NOTE:
%            regardless of the value of the space option, the measurement
%            frame is always related to the RAS space according to the NRRD
%            standard convention)

space = img.metaData.space; % something light right-anterior-superior, nothing to parse
origin = img.metaData.space_origin; % something like: (100.0,-100.0,-23.1)
direction = img.metaData.space_directions; % something like: none (2.0,0.1,0.1) (0.1,-2.0,-0.1) (0.2,-0.3,4.0)
measurement_frame = img.metaData.measurement_frame; % something like: (2.0,0.1,0.1) (0.1,-2.0,-0.1) (0.2,-0.3,4.0)

i1 = strfind(origin,'(');
i2 = strfind(origin,')');
if( length(i1)~=1 || length(i2)~=1 )
    error('Missformed origin string');
else
    origin = sscanf(origin(i1+1:i2-1),'%f,%f,%f');
end

i1 = strfind(direction,'(');
i2 = strfind(direction,')');
if( length(i1)~=3 || length(i2)~=3 )
    error('Missformed direction string');
else
    D = eye(3,3);
    for c=1:3
        D(:,c) = sscanf(direction(i1(c)+1:i2(c)-1),'%f,%f,%f');
    end
    direction = D;
end

i1 = strfind(measurement_frame,'(');
i2 = strfind(measurement_frame,')');
if( length(i1)~=3 || length(i2)~=3 )
    error('Missformed measurement frame string');
else
    D = eye(3,3);
    for c=1:3
        D(:,c) = sscanf(measurement_frame(i1(c)+1:i2(c)-1),'%f,%f,%f');
    end
    measurement_frame = D;
end
