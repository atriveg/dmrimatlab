function nii = dwinrrd2nii( nrrdfile, niifile )
% function dwinrrd2nii( nrrdfile, niifile )
% function nii = dwinrrd2nii( nrrdfile, niifile )
%
% Translates a NRRD file <nrrdfile> compatible with 3-D Slice containing a
% DWI volume to a nifti file <niifile> with extension .nii or .nii.gz.
% The gradientes table (b-values and b-vectors) are stored in two
% additional files with extensions .bval and .bvec (both are plain text
% files with an Nx1 vector and an Nx3 matrix, respectively) and with the
% same basemane as <niifile>. Optionally, an image structure <nii> is
% returned similar to that read by load_untouch_nii.

if(nargin~=2)
   error('Supply just the name of the input NRRD and the output nii files');
end

% Load NRRD file:
try
    nrrd = nrrdread(nrrdfile);
catch ME
    warning('Could not open NRRD file <%s>',nrrdfile);
    rethrow(ME);
end

% Create nii structure:
[~,name,~] = fileparts(nrrdfile);
nii = create_nii_struct(nrrd,name);

% Save nii file:
try
    save_untouch_nii(nii,niifile);
catch ME
    warning('Could not save nii file <%s>',niifile);
    rethrow(ME);
end

% Save bvec and bval files:
[path,name,ext] = fileparts(niifile);
if(strcmp(ext,'.gz'))
    [~,name,~] = fileparts(name);
end
bvecfile = sprintf('%s.bvec',name);
bvalfile = sprintf('%s.bval',name);
if(~isempty(path))
    bvecfile = sprintf('%s%s%s',path,filesep,bvecfile);
    bvalfile = sprintf('%s%s%s',path,filesep,bvalfile);
end
[~,gi,bi] = slicervol2dwi(nrrd);
save('-ascii',bvecfile,'gi');
save('-ascii',bvalfile,'bi');
end

% -------------------------------------------------------------------------
% Check:
%  https://afni.nimh.nih.gov/pub/dist/doc/nifti/nifti1_h.pdf
% for details
function nii = create_nii_struct(nrrd,filename)
% First, check which dimensions correspond to the xyz-space and which ones
% correspond to the q-space:
[X,Y,Z,G,nrrd] = check_phisical_dimensions(nrrd);
% Populate the Nifti header as needed:
% ---------------------------------
nii.hdr.hk.sizeof_hdr = 348;   % Always 348
nii.hdr.hk.data_type = '';     % UNUSED
nii.hdr.hk.db_name = '';       % UNUSED
nii.hdr.hk.extents = 0;        % UNUSED
nii.hdr.hk.session_error = 0;  % UNUSED
nii.hdr.hk.regular = 'r';      % UNUSED
nii.hdr.hk.dim_info = 0;       % NON-ESSENTIAL
% ---------------------------------
nii.hdr.dime.dim = [4,X,Y,Z,G,1,1,1];
nii.hdr.dime.intent_p1 = 0;    % UNUSED FOR DWI
nii.hdr.dime.intent_p2 = 0;    % UNUSED FOR DWI
nii.hdr.dime.intent_p3 = 0;    % UNUSED FOR DWI
nii.hdr.dime.intent_code = 0;  % UNUSED FOR DWI
% Check data precission:
[code,bpix] = datatype_code(nrrd.pixelData(1));
nii.hdr.dime.datatype = code;
nii.hdr.dime.bitpix = bpix;
nii.hdr.dime.slice_start = 0;
nii.hdr.dime.pixdim = compute_pix_dims(nrrd.ijkToLpsTransform);
% Directly taken from the HCP-MGH data (what the +4 comes from?):
nii.hdr.dime.vox_offset = nii.hdr.hk.sizeof_hdr + 4;
nii.hdr.dime.scl_slope = 1; % *1 linear scaling
nii.hdr.dime.scl_inter = 0; % +0 linear scaling
nii.hdr.dime.slice_end = 0; % not used if slice_code = 0
nii.hdr.dime.slice_code = 0;
% The NRRD header allows to define the units in a very flexible way, even
% using different units for each dimension:
%    https://teem.sourceforge.net/nrrd/format.html#units
% Nifti is not equally flexible, since it just allows using either meters,
% millimeters or unknown units. We will assume millimeters:
nii.hdr.dime.xyzt_units = 10;
nii.hdr.dime.cal_max = 0;
nii.hdr.dime.cal_min = 0;
nii.hdr.dime.slice_duration = 0;
nii.hdr.dime.toffset = 0;
nii.hdr.dime.glmax = max(nrrd.pixelData(:));
nii.hdr.dime.glmin = min(nrrd.pixelData(:));
% ---------------------------------
nii.hdr.hist.descrip = sprintf('Converted with dwinrrd2nii');
nii.hdr.hist.aux_file = '';
% Compute the homogeneous matrix to go from pixel space to physical space:
T = compute_ijk2xyz_matrix( nrrd.ijkToLpsTransform, nrrd.metaData.space );
nii.hdr.hist.qform_code = 1;
nii.hdr.hist.sform_code = 1;
% Compute quaternion form from homogeneous matrix:
q = rotmatrix2quaternion(T);
nii.hdr.hist.quatern_b = q(2);
nii.hdr.hist.quatern_c = q(3);
nii.hdr.hist.quatern_d = q(4);
nii.hdr.hist.qoffset_x = T(1,4);
nii.hdr.hist.qoffset_y = T(2,4);
nii.hdr.hist.qoffset_z = T(3,4);
% The homogeneous transform is:
nii.hdr.hist.srow_x = T(1,:);
nii.hdr.hist.srow_y = T(2,:);
nii.hdr.hist.srow_z = T(3,:);
% ---
nii.hdr.hist.intent_name = '';
nii.hdr.hist.magic = 'n+1';
% ---------------------------------
nii.filetype = 2;
nii.fileprefix = filename;
nii.machine = 'ieee-le';
nii.ext = [];
nii.img = nrrd.pixelData;
nii.untouch = 1;
end

% -------------------------------------------------------------------------
function [X,Y,Z,G,nrrd] = check_phisical_dimensions(nrrd)
% This is a cell array with the kinds for each of the dimensions:
kinds = strsplit(nrrd.metaData.kinds);
assert(length(kinds)==4,'This doesn''t look like a DWI volume');
% Look for a keyword other than 'domain' or 'space'
pp = ( strcmp('domain',kinds) | strcmp('space',kinds) );
pp = find(~pp);
assert(isscalar(pp),'Wrong ''kinds'' field in metadata');
% Keep the actual dimensions corresponding to the physical extent:
pdims = setdiff(1:4,pp);
% Permute the image buffer as needed:
nrrd.pixelData = permute( nrrd.pixelData, [pdims,pp] );
% Return the proper dimension:
[X,Y,Z,G] = size(nrrd.pixelData);
end

% -------------------------------------------------------------------------
function [code,bpix] = datatype_code(sample)
% Return a numeric code for each data type as described in:
%   https://afni.nimh.nih.gov/pub/dist/doc/nifti/nifti1_h.pdf
switch(class(sample))
    case 'char'
        code = 256;
        bpix = 8;
    case 'uchar'
        code = 2;
        bpix = 8;
    case 'uint8'
        code = 256;
        bpix = 8;
    case 'int8'
        code = 256;
        bpix = 8;
    case 'unit16'
        code = 512;
        bpix = 16;
    case 'int16'
        code = 4;
        bpix = 16;
    case 'uint32'
        code = 768;
        bpix = 32;
    case 'int32'
        code = 8;
        bpix = 32;
    case 'uint64'
        code = 1280;
        bpix = 64;
    case 'int64'
        code = 1024;
        bpix = 64;
    case 'single'
        code = 16;
        bpix = 32;
    case 'double'
        code = 64;
        bpix = 64;
    otherwise
        error('Unsupported data type');
end
end

% -------------------------------------------------------------------------
function pixdim = compute_pix_dims(ijk2lps)
ijk   = [ [0;0;0;1], [1;0;0;1], [0;1;0;1], [0;0;1;1] ];
lps   = ijk2lps*ijk;
delta = lps(1:3,2:4) - lps(1:3,1);
delta = sqrt(sum(delta.*delta,1));
% What does the pixdim of the q-space dimension stand for??? In the HCP_MGH
% volume it is set to 8.8, but there seems to be not a proper documentation
% for that:
pixdim = [1,delta,8.8,0,0,0];
end

% -------------------------------------------------------------------------
function ijk2xyz = compute_ijk2xyz_matrix( ijk2xyz, space )
% Note: nrrdread loads the NRRD file and sets the ijkToLpsTransform in
% one-based indexing, but we need it to be zero-based:
onebased2zerobased = [ [eye(3), [1;1;1] ]; [0 0 0 1] ];
ijk2xyz            = ijk2xyz*onebased2zerobased;
% Now, we have to take into account the anatomical space of the image. Note
% NRRD has a specific field "space" for it, something like:
%   right-anterior-superior (for RAS) 
% but Nifti has not. First, parse the "space" field:
if(length(space)>3)
    ds = strsplit( space, '-' );
    assert(length(ds)==3,'Wrong space field in the NRRD header');
    dx = lower(ds{1}(1));
    dy = lower(ds{2}(1));
    dz = lower(ds{3}(1));
elseif(length(space)==3)
    dx = lower(space(1));
    dy = lower(space(2));
    dz = lower(space(3));
else
    error('Wrong space field in the NRRD header');
end
% Now, dx can be either 'r' or 'l'
% Now, dy can be either 'a' or 'p'
% Now, dz can be either 's' or 'i'
% Nifti assumes RAS world coordinates, hence we need to check each row of
% the homogeneous matrix to make sure it is consistent with the space used
% in the NRRD file:
if(dx~='r') % Need to correct
    ijk2xyz(1,:) = -ijk2xyz(1,:);
end
if(dy~='a') % Need to correct
    ijk2xyz(2,:) = -ijk2xyz(2,:);
end
if(dz~='s') % Need to correct
    ijk2xyz(3,:) = -ijk2xyz(3,:);
end
end

% -------------------------------------------------------------------------
function q = rotmatrix2quaternion(T)
% --------------------------------------------
% Keep just the rotation part of T:
T = T(1:3,1:3);
% Normalize each column:
n = sqrt(sum(T.*T,1));
T = T./n;
% --------------------------------------------
m00 = T(1,1); m01 = T(1,2); m02 = T(1,3);
m10 = T(2,1); m11 = T(2,2); m12 = T(2,3);
m20 = T(3,1); m21 = T(3,2); m22 = T(3,3);
% --------------------------------------------
if( m22 < 0 )
    if ( m00 > m11 )
        t = 1 + m00 - m11 - m22;
        q = [ m21-m12, t, m01+m10, m20+m02 ];
    else
        t = 1 - m00 + m11 - m22;
        q = [ m02-m20, m01+m10, t, m12+m21  ];
    end
else
    if( m00 < -m11 )
        t = 1 - m00 - m11 + m22;
        q = [ m10-m01, m20+m02, m12+m21, t ];
    else
        t = 1 + m00 + m11 + m22;
        q = [ -t, m12-m21, m20-m02, m01-m10 ];
    end
end
% --------------------------------------------
q = q/norm(q);
if(q(1)<0)
    q = -q;
end
% --------------------------------------------
end
