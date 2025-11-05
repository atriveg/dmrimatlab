function nii = saveniiasin( vol, niifilename, template )
% function nii = saveniiasin( vol, niifilename, template )
%
% This function is intended to ease the writing of data volumes processed
% by the toolbox to nifti files that can be read and displayed by some
% other software. In brief:
%
%    vol: is a X x Y x Z or a X x Y x Z x W 3-D or 4-D array coming from a
%        certain processing pipeline (for example: a FA map, a denoised DWI
%        volume, a volume containing spherical harmonic coefficients...).
%
%    niifilename: is the filename of the nifti file to be written.
%
%    template: is either a structure as returned by "load_untouch_nii" or
%        the name of a nifti file (with extension .nii or .nii.gz) whose
%        anatomical information matches that of <vol>.
%    
%    nii: is the nifti structure to be written to niifilename, virtually
%        identical to the result of running:
%            >> nii = load_untouch_nii(niifilename)
%
%    EXAMPLE: Imagine you have a DWI volume in a nifti file named
%    'dwi.nii.gz'. Then you can extract its first channel to a new nifti
%    file as follows:
%
%       >> dwi = load_untouch_nii('dwi.nii.gz');
%       >> first = dwi.img(:,:,:,1);
%       >> saveniiasin( first, 'firstchannel.nii.gz', dwi );
%
%     Alterantively, the last line can also read:
%
%       >> saveniiasin( first, 'firstchannel.nii.gz', 'dwi.nii.gz' );
%
%     but it will take longer since it will have to read 'dwi.nii.gz'
%     again.

if(nargin~=3)
   error('Supply just the data volume, the file name and a template');
end

if(ischar(template))
    try
        template = load_untouch_nii(template);
    catch ME
        warning('Could not open nii file <%s>',template);
        rethrow(ME);
    end
end

% Create nii structure:
nii = replicate_nii_struct(vol, niifilename, template);

% Save nii file:
try
    save_untouch_nii(nii,niifilename);
catch ME
    warning('Could not save nii file <%s>',niifilename);
    rethrow(ME);
end

end

% -------------------------------------------------------------------------
% Check:
%  https://afni.nimh.nih.gov/pub/dist/doc/nifti/nifti1_h.pdf
% for details
function nii = replicate_nii_struct(vol, niifilename, template )
% Make sure the size of vol is compatible with the template:
nd = ndims(vol);
if( (nd<2) || (nd>7) )
    error('The input volume has an unsupported number of dimensions');
end
[X,Y,Z] = size(vol);
assert( isequal([X,Y,Z],template.hdr.dime.dim(2:4)), ...
    'The size of your volume doesn''t match the FoV of the template' );
% Most of the information can be just copied from the template:
nii = template;
% So that we just need to change the information that is actually
% different:
sz  = [ ndims(vol), size(vol) ];
sz(ndims(vol)+2:8) = 1;
nii.hdr.dime.dim   = sz;
% ---------------------------------
[code,bpix] = datatype_code(vol(1));
nii.hdr.dime.datatype = code;
nii.hdr.dime.bitpix = bpix;
nii.hdr.dime.glmax = max(vol(:));
nii.hdr.dime.glmin = min(vol(:));
% ---------------------------------
nii.hdr.hist.descrip = '';
nii.hdr.hist.aux_file = '';
% ---------------------------------
nii.fileprefix = niifilename;
nii.img = vol;
nii.untouch = 1;
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
