function nrrd = dwinii2nrrd( niifile, nrrdfile, bvecfile, bvalfile )
% function dwinii2nrrd( niifile, nrrdfile )
% function dwinii2nrrd( niifile, nrrdfile, bvecfile )
% function dwinii2nrrd( niifile, nrrdfile, bvecfile, bvalfile )
% function nrrd = dwinii2nrrd( niifile, nrrdfile )
% function nrrd = dwinii2nrrd( niifile, nrrdfile, bvecfile )
% function nrrd = dwinii2nrrd( niifile, nrrdfile, bvecfile, bvalfile )
%
% Translates a nifti file <niifile> containing a DWI volume, with a
% gradient table described by b-vectors in <bvecfile> and b-values
% described in <bvalfile>, to a nrrd/nhdr file <nrrdfile> which can be
% loaded in 3D-Slicer. If <bvecfile> and/or <bvalfile> are ommitted, the
% function will look for their standard names. In any case the function can
% return a structure <nrrd> as loaded by nrrdread().
nrrd = [];
bth  = 25;
% --------------------------------------------
try
    nii = load_untouch_nii(niifile);
catch
    error('Could not open nii file <%s>',niifile);
end
% --------------------------------------------
if(nargin<3)
    [path,name,ext] = fileparts(niifile);
    if(strcmp(ext,'.gz'))
        name = strrep(name,'.nii','');
    end
    if(isempty(path))
        bvecfile = sprintf('%s.bvec',name);
    else
        bvecfile = sprintf('%s%s%s.bvec',path,filesep,name);
    end
    if(exist(bvecfile,'file')~=2)
        error('You didn''t provide a b-vec file and I cannot find <%s> either',bvecfile);
    end
end
try
    gi = load(bvecfile);
catch
    error('Could not load gradients file <%s>',bvecfile);
end
if( ~isnumeric(gi) || ~ismatrix(gi) )
    error('Miss-formed gradients file <%s>',bvecfile);
end
if( size(gi,2)~=3 )
    gi = gi';
end
if( size(gi,2)~=3 )
    error('Miss-formed gradients file <%s>',bvecfile);
end
% --------------------------------------------
if(nargin<4)
    [path,name,ext] = fileparts(niifile);
    if(strcmp(ext,'.gz'))
        name = strrep(name,'.nii','');
    end
    if(isempty(path))
        bvalfile = sprintf('%s.bval',name);
    else
        bvalfile = sprintf('%s%s%s.bval',path,filesep,name);
    end
    if(exist(bvalfile,'file')~=2)
        error('You didn''t provide a b-val file and I cannot find <%s> either',bvalfile);
    end
end
try
    bi = load(bvalfile);
catch
    error('Could not load b-vals file <%s>',bvalfile);
end
if( ~isnumeric(bi) || ~ismatrix(bi) )
    error('Miss-formed b-cals file <%s>',bvalfile);
end
bi = bi(:);
assert( size(bi,1)==size(gi,1), 'The size of the gradients table does not match the number of b-values' );
% --------------------------------------------
assert( size(gi,1)==size(nii.img,4), 'The size of the gradients table does not match the 4-th dimension of the nifti volume')
% --------------------------------------------
[~,name,ext] = fileparts(nrrdfile);
switch(ext)
    case '.nhdr'
        NHDR = true;
        datafile = sprintf('%s.raw.gz',name);
    case '.nrrd'
        NHDR = false;
    otherwise
        error('The only extensions allowed for the output file are .nhdr or .nrrd, got <%s>',ext);
end
% --------------------------------------------
ni = sqrt(sum(gi.*gi,2));
bs = (ni<10*eps);
ni = 1./ni;
ni(bs) = 0;
gi = gi.*ni; % All have norm 1 but baselines, which have norm 0
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
bval0 = max(bi); %
ni    = sqrt(bi/bval0);
gi    = gi.*ni;
% --------------------------------------------
[bs,ps,Ns] = auto_detect_shells(bi,bth);
shelldesc = '';
for n=1:Ns
    count     = sum(ps==n);
    if(n==1)
        shelldesc = sprintf('%1.1f x %d',bs(n),count);
    else
        shelldesc = sprintf('%s, %1.1f x %d',shelldesc,bs(n),count);
    end
end
% --------------------------------------------
data = nii.img;
[typesl,typemat] = classconvert(data);
data = cast(data,typemat);
nrrd.pixelData = permute(data,[4,1,2,3]);
% --------------------------------------------
nrrd.ijkToLpsTransform(1,:) = nii.hdr.hist.srow_x;
nrrd.ijkToLpsTransform(2,:) = nii.hdr.hist.srow_y;
nrrd.ijkToLpsTransform(3,:) = nii.hdr.hist.srow_z;
% --------------------------------------------
nrrd.metaData.type = typesl;
nrrd.metaData.dimension = '4';
nrrd.metaData.space = 'right-anterior-superior';
nrrd.metaData.sizes = sprintf('%d %d %d %d', ...
    size(data,1), size(data,2), size(data,3), size(data,4));
nrrd.metaData.space_directions = sprintf('none (%1.14f,%1.14f,%1.14f) (%1.14f,%1.14f,%1.14f) (%1.14f,%1.14f,%1.14f)', ...
    nrrd.ijkToLpsTransform(1,1), nrrd.ijkToLpsTransform(2,1), nrrd.ijkToLpsTransform(3,1),...
    nrrd.ijkToLpsTransform(1,2), nrrd.ijkToLpsTransform(2,2), nrrd.ijkToLpsTransform(3,2),...
    nrrd.ijkToLpsTransform(1,3), nrrd.ijkToLpsTransform(2,3), nrrd.ijkToLpsTransform(3,3)    );
nrrd.metaData.kinds = 'list domain domain domain';
nrrd.metaData.endian = 'little';
nrrd.metaData.encoding = 'gzip';
nrrd.metaData.space_origin = sprintf('%1.14f,%1.14f,%1.14f', ...
    nii.hdr.hist.qoffset_x, nii.hdr.hist.qoffset_y, nii.hdr.hist.qoffset_z );
if(NHDR)
    nrrd.metaData.data_file = datafile;
end
nrrd.metaData.measurement_frame = '(1.0,0.0,0.0) (0.0,1.0,0.0) (0.0,0.0,1.0)';
nrrd.metaData.DWMRI_comments = sprintf('Converted from %s',niifile);
nrrd.metaData.DWMRI_shells = shelldesc;
nrrd.metaData.DWMRI_b_value = sprintf('%1.9f',bval0);
for g=0:size(gi,1)-1
    nrrd.metaData.(sprintf('DWMRI_gradient_%s',padgnumber(g))) = sprintf('%1.6f %1.6f %1.6f', ...
        gi(g+1,1), gi(g+1,2), gi(g+1,3) );
end
nrrd.metaData.modality = 'DWMRI';
% --------------------------------------------
nrrd.metaDataFieldNames.space_directions = 'space directions';
nrrd.metaDataFieldNames.space_origin = 'space origin';
nrrd.metaDataFieldNames.measurement_frame = 'measurement frame';
if(NHDR)
    nrrd.metaDataFieldNames.data_file = 'data file';
end
nrrd.metaDataFieldNames.DWMRI_b_value = 'DWMRI_b-value';
% --------------------------------------------
nrrdwrite_patched(nrrdfile,nrrd,NHDR);
% --------------------------------------------
end

% -------------------------------------------------------------------------
function str = padgnumber(g)
if(g<10)
    pad = '000';
elseif(g<100)
    pad = '00';
elseif(g<1000)
    pad = '0';
else
    pad = '';
end
str = sprintf('%s%d',pad,g);
end
% -------------------------------------------------------------------------
function [strsl,strmat] = classconvert(vol)
strmat = class(vol);
switch(strmat)
    case 'double'
        strsl = 'double';
    case 'single'
        strsl = 'float';
    case 'uint16'
        strsl = 'unsigned short';
    case 'int16'
        strsl = 'short';
    case 'uint8'
        strsl = 'unsigned char';
    case 'int8'
        strsl = 'char';
    otherwise
        strmat = 'int16';
        strsl  = 'short';
end
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function nrrdwrite_patched(outputFilename, img, NHDR )
% Write image and metadata to a NRRD file (see http://teem.sourceforge.net/nrrd/format.html)
%   img.pixelData: pixel data array
%   img.ijkToLpsTransform: pixel (IJK) to physical (LPS, assuming 'space' is 'left-posterior-superior')
%     coordinate system transformation, the origin of the IJK coordinate system is (1,1,1) to match Matlab matrix indexing
%   img.metaData: Contains all the descriptive information in the image header. The following fields are ignored:
%     sizes: computed to match size of img.pixelData
%     type: computed to match type of img.pixelData
%     kinds: computed to match dimension of img.pixelData
%     dimension: computed to match dimension of img.pixelData
%     space_directions: ignored if img.ijkToLpsTransform is defined
%     space_origin: ignored if img.ijkToLpsTransform is defined
%   img.metaData: Contains the list of full NRRD field names for each
%     metaData field name. All fields should be listed here that have a
%     special character in their name (such as dot or space).
%   img.metaDataFieldNames: Contains full names of metadata fields that cannot be used as Matlab field names because they contains
%     special characters (space, dot, etc). Full field names are used for determining the field name to be used in the NRRD file
%     from the Matlab metadata field name.
%
% Supports writing of 3D and 4D volumes.
%
% Examples:
%
% 1. Using output from nrrdread:
%
%   img = nrrdread('testData\MRHeadRot.nrrd')
%   nrrdwrite('testOutput.nrrd', img)
%
% 2. Creating volume from scratch - minimal example
%
%   [x,y,z] = meshgrid([-10:10],[-12:15],[-8:6]);
%   img.pixelData = x/3+y/4+z/2;
%
%   nrrdwrite('testOutput.nrrd', img);
%
% 3. Creating volume from scratch
%
%   % Set pixel data
%   [x,y,z] = meshgrid([-10:10],[-12:15],[-8:6]);
%   img.pixelData = x/3+y/4+z/2;
%
%   % Define origin, spacing, axis directions by a homogeneous transformation matrix:
%   img.ijkToLpsTransform = [ 1.2 0 0 10; 0 1.2 0 12; 0 0 3.0 -22; 0 0 0 1];
%
%   % Enable compression
%   img.metaData.encoding='gzip';
%
%   nrrdwrite('testOutput.nrrd', img);
%

% Open file for writing
fid=fopen(outputFilename, 'w');
if(fid<=0)
    fprintf('Could not open file: %s\n', outputFilename);
end

standardFieldNames = { 'data file', 'type', 'dimension', 'space', 'sizes', 'space directions', 'kinds', 'endian', 'encoding', 'space origin', 'measurement frame' };
dwmriFieldNames = {'DMRI_b-value'};

fprintf(fid,'NRRD0005\n');
fprintf(fid,'# Complete NRRD file format specification at:\n');
fprintf(fid,'# http://teem.sourceforge.net/nrrd/format.html\n');

% Create/override mandatory fields

img.metaData.type = getMetaType(class(img.pixelData));
img.metaData.dimension = length(size(img.pixelData)); % ndim is not defined for int16 arrays
if ~isfield(img.metaData,'space')
    img.metaData.space = 'left-posterior-superior';
end
img.metaData.sizes=num2str(size(img.pixelData));

if isfield(img,'ijkToLpsTransform')
    % Write zero-based IJK transform (origin is at [0,0,0]) to the image header
    ijkOneBasedToLpsTransform=img.ijkToLpsTransform;
    %ijkOneBasedToIjkZeroBasedTransform=[[eye(3), [-1;-1;-1] ]; [0 0 0 1]];
    ijkOneBasedToIjkZeroBasedTransform=[[eye(3), [0;0;0] ]; [0 0 0 1]];
    ijkZeroBasedToLpsTransform=ijkOneBasedToLpsTransform*inv(ijkOneBasedToIjkZeroBasedTransform);
    axes_origin=ijkZeroBasedToLpsTransform(1:3,4);
    img.metaData.space_origin=sprintf('(%f,%f,%f)',reshape(axes_origin,1,3));
    axes_directions=ijkZeroBasedToLpsTransform(1:3,1:3);
    switch (img.metaData.dimension)
        case {3}
            img.metaData.space_directions=sprintf('(%f,%f,%f) (%f,%f,%f) (%f,%f,%f)',reshape(axes_directions,1,9));
        case {4}
            img.metaData.space_directions=sprintf('none (%f,%f,%f) (%f,%f,%f) (%f,%f,%f)',reshape(axes_directions,1,9));
        otherwise
            assert(false, 'Unsupported pixel data dimension')
    end
end

if ~isfield(img.metaData,'space_directions')
    switch (img.metaData.dimension)
        case {3}
            img.metaData.space_directions = '(1,0,0) (0,1,0) (0,0,1)';
        case {4}
            img.metaData.space_directions = 'none (1,0,0) (0,1,0) (0,0,1)';
        otherwise
            assert(false, 'Unsupported pixel data dimension')
    end
end

switch (img.metaData.dimension)
    case {3}
        img.metaData.kinds='domain domain domain';
    case {4}
        if(~isfield(img.metaData,'kinds'))
            img.metaData.kinds='list domain domain domain';
        end
        % Add a custom field to make the volume load into 3D Slicer as a MultiVolume
        isMultiVolume = true;
        if(isfield(img.metaData,'modality'))
            if(strcmp(img.metaData.modality,'DWMRI'))
                isMultiVolume = false;
            end
        end
        if(isMultiVolume)
            img = nrrdaddmetafield(img,'MultiVolume.NumberOfFrames',size(img.pixelData,4));
        end
    otherwise
        assert(false, 'Unsupported pixel data dimension')
end

if ~isfield(img.metaData,'endian')
    img.metaData.endian='little';
end

if ~isfield(img.metaData,'encoding')
    img.metaData.encoding='raw';
end

% Make sure that standard field names that contain special
% characters have their full field names defined.
for k=1:length(standardFieldNames)
    fullFieldName=standardFieldNames{k};
    fieldName=regexprep(fullFieldName,'\W','_');
    if ~strcmp(fieldName,fullFieldName)
        img.metaDataFieldNames.(fieldName)=fullFieldName;
    end
end

for k=1:length(dwmriFieldNames)
    fullFieldName=dwmriFieldNames{k};
    fieldName=regexprep(fullFieldName,'\W','_');
    fieldName=strrep(fieldName,'-','_');
    if ~strcmp(fieldName,fullFieldName)
        img.metaDataFieldNames.(fieldName)=fullFieldName;
    end
end



% Print the header data to the output file
metaDataCellArr = struct2cell(img.metaData);
fields = fieldnames(img.metaData);
for i=1:numel(fields)
    writeFieldName(fid, fields{i}, img.metaDataFieldNames, standardFieldNames);
    writeDataByType(fid,metaDataCellArr{i});
end

fprintf(fid,'\n');

if(NHDR)
    fclose(fid);
    fid = fopen(img.metaData.data_file,'w');
end

% Write pixel data
switch (img.metaData.encoding)
    case {'raw'}
        fwrite(fid, img.pixelData, class(img.pixelData));
    case {'gzip', 'gz'}
        try
            compressedPixelData = nrrd_gzip_compress(img.pixelData, class(img.pixelData));
        catch ME
            fclose('all');
            rethrow(ME);
        end
        fwrite( fid, compressedPixelData, 'uint8' );
    otherwise
        assert(false, 'Unsupported encoding')
end

fclose('all');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeFieldName(fid, fieldName, fullFieldNames, standardFieldNames)
% If full field name is listed in img.metaDataFieldNames then use that
% instead of the Matlab field name.
if isfield(fullFieldNames,fieldName)
    fullFieldName = fullFieldNames.(fieldName);
else
    fullFieldName = fieldName;
end

isStandardFieldName = ~isempty(find(strcmp(fullFieldName, standardFieldNames), 1));
if isStandardFieldName
    % Standard field names are separated by :
    fprintf(fid,'%s: ',fullFieldName);
else
    % Custom field names are separated by :=
    fprintf(fid,'%s:=',fullFieldName);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeDataByType(fid, data)
% Function that writes the header data to file based on the type of data
% params: - fid of file to write to
%         - data from header to write
if ischar(data)
    fprintf(fid,'%s\n',data);
else
    fprintf(fid,'%d ',data);
    fprintf(fid,'\n');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function metaType = getMetaType(matlabType)
% Determine the metadata type from the Matlab type
switch (matlabType)
    case {'int8'}
        metaType = 'int8';
    case {'uint8'}
        metaType = 'uint8';
    case {'int16'}
        metaType = 'int16';
    case {'uint16'}
        metaType = 'uint16';
    case {'int32'}
        metaType = 'int32';
    case {'uint32'}
        metaType = 'uint32';
    case {'int64'}
        metaType = 'int64';
    case {'uint64'}
        metaType = 'uint64';
    case {'single'}
        metaType = 'float';
    case {'double'}
        metaType = 'double';
    otherwise
        assert(false, 'Unsupported Matlab data type')
end
end

