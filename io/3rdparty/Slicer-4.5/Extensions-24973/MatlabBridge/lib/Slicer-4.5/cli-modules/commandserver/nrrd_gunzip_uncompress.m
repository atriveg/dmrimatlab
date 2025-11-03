function pixelData = nrrd_gunzip_uncompress( compressedPixelData, DataType )

cn = find(strcmp(DataType,{'double','single','logical','char','int8','uint8',...
    'int16','uint16','int32','uint32','int64','uint64'}));

if cn == 3 || cn == 4
    DataType  = 'uint8';
end

% Create a temporary folder to work in:
tmpDir = tempname;
try
    mkdir(tmpDir);
catch
    error('Unable to create tmp folder');
end

% Create a tmp file to write zipped data:
gzipFileName = sprintf('%s%snrrdZippedData.bin.gz',tmpDir,filesep);
fid = fopen(gzipFileName,'w');
if(fid<0)
    error('Unable to write compressed data to tmp file');
end

% Save zipped byte stream:
fwrite( fid, compressedPixelData, 'uint8' );
fclose(fid);

% gunzip the zipped data to a tmp file:
try
    rawFileName = gunzip(gzipFileName);
catch
    error('Unable to gunzip compressed data from the tmp file');
end

% Read the uncompressed file to an array:
fid = fopen(rawFileName{1},'r');
if(fid<0)
    error('Unable to read uncompressed data from tmp file');
end
pixelData = fread( fid, Inf, DataType );
fclose(fid);

pixelData = cast(pixelData,DataType);

% Clean-up:
sf = check_software_platform;
if(sf==2) % Octave
    confirm_recursive_rmdir( false, 'local' );
end
rmdir(tmpDir,'s');

end

